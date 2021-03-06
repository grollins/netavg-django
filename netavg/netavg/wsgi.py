import os
from os.path import abspath, dirname
from sys import path

SITE_ROOT = dirname(dirname(abspath(__file__)))
path.append(SITE_ROOT)

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "netavg.settings.production")
# os.environ.setdefault("DJANGO_SETTINGS_MODULE", "netavg.settings.local")

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()

