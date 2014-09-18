from django.conf.urls import patterns, url, include
from django.contrib import admin
from django.http import HttpResponse
admin.autodiscover()

from rest_framework.routers import DefaultRouter
from jobs.views import JobViewSet, UserViewSet, ResultViewSet


router = DefaultRouter()
router.register(r'jobs', JobViewSet, base_name='job')
router.register(r'results', ResultViewSet)
router.register(r'users', UserViewSet)

urlpatterns = patterns('',
    url(r'^', include(router.urls)),
    url(r'^login/$', 'jobs.views.login', name='login')
)

urlpatterns += patterns('backend.views',
    url(r'^admin/', include(admin.site.urls))
)

urlpatterns += patterns('',
    (r'^robots\.txt$',
     lambda r: HttpResponse("User-agent: *\nDisallow: /", mimetype="text/plain"))
)
