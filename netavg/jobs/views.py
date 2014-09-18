import logging
import requests

from django.contrib.auth import authenticate
from django.contrib.auth.models import User
from django.conf import settings

from rest_framework import status, permissions, viewsets
from rest_framework.response import Response
from rest_framework.decorators import api_view, permission_classes
from rest_framework.decorators import authentication_classes
from rest_framework.authtoken.models import Token
from django_browserid.base import BrowserIDException

from .models import Job, Result
from .permissions import IsOwner
from .serializers import JobSerializer, UserSerializer, ResultSerializer
from .tasks import run_netavg_calculation

login_logger = logging.getLogger('login.audience')

class JobViewSet(viewsets.ModelViewSet):
    serializer_class = JobSerializer
    permission_classes = (permissions.IsAuthenticated, IsOwner)

    def pre_save(self, obj):
        obj.owner = self.request.user
        obj.status = Job.STATUS.submitted

    def post_save(self, obj, created=False):
        if created:
            try:
                run_netavg_calculation.delay(obj.id)
            except:
                obj.status = Job.STATUS.error
                obj.save()
                raise

    def get_queryset(self):
            """
            This view should return a list of all the jobs
            for the currently authenticated user.
            """
            user = self.request.user
            return Job.objects.filter(owner=user)


class ResultViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Result.objects.all()
    serializer_class = ResultSerializer
    permission_classes = (permissions.IsAuthenticated,)

    def get_queryset(self):
            """
            This view should return a list of all the jobs
            for the currently authenticated user.
            """
            user = self.request.user
            return Result.objects.filter(job__owner=user)


class UserViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = User.objects.all()
    serializer_class = UserSerializer
    permission_classes = (permissions.IsAdminUser,)


@api_view(['POST'])
@authentication_classes([])
@permission_classes([])
def login(request):
    assertion = request.DATA.get('assertion', None)
    if not assertion:
        return Response(
            'assertion parameter is missing',
            status.HTTP_400_BAD_REQUEST
        )
    audience = settings.API_AUDIENCE
    login_logger.error(audience)
    try:
        page = requests.post('https://verifier.login.persona.org/verify',
                             verify=True,
                             data={ "assertion": assertion,
                                    "audience": audience})
        data = page.json()
    except requests.exceptions.SSLError:
        login_logger.error("Failed login\n%s\n" % request.META)
        return Response(
            'SSL failure',
            status.HTTP_501_NOT_IMPLEMENTED,
        )
    except requests.exceptions.ConnectionError:
        login_logger.error("Failed login\n%s\n" % request.META)
        return Response(
            'Failed to contact authentication server',
            status.HTTP_502_BAD_GATEWAY,
        )

    try:
        user = authenticate(
            assertion=assertion,
            audience=audience,
        )
        if user:
            return Response({
                'email': user.email,
                'token': get_auth_token(user),
            })
        else:
            return Response(
                'Authentication failed',
                status.HTTP_500_INTERNAL_SERVER_ERROR,
            )
    except BrowserIDException as e:
        login_logger.error("Failed login\n%s\n" % request.META)
        login_logger.error("%s" % type(e.exc))
        return Response(
            'Failed to contact authentication server',
            status.HTTP_503_SERVICE_UNAVAILABLE,
        )

def get_auth_token(user):
    token, created = Token.objects.get_or_create(user=user)
    return token.key
