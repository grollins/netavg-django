"""Production settings and globals."""

from os import environ

from base import *

# Normally you should not import ANYTHING from Django directly
# into your settings, but ImproperlyConfigured is an exception.
from django.core.exceptions import ImproperlyConfigured

########## HOST CONFIGURATION
# See: https://docs.djangoproject.com/en/1.5/releases/1.5/#allowed-hosts-required-in-production
ALLOWED_HOSTS = ['127.0.0.1', '104.131.42.27']
########## END HOST CONFIGURATION

########## EMAIL CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#email-backend
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
########## END EMAIL CONFIGURATION

########## DATABASE CONFIGURATION
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': get_env_setting('NETAVG_DB_NAME'),
        'USER': get_env_setting('NETAVG_DB_USER'),
        'PASSWORD': get_env_setting('NETAVG_DB_PWD'),
        'HOST': get_env_setting('NETAVG_DB_HOST'),
        'PORT': get_env_setting('NETAVG_DB_PORT'),
    }
}
########## END DATABASE CONFIGURATION


########## CACHE CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#caches
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
    }
}
########## END CACHE CONFIGURATION


########## SECRET CONFIGURATION
# See: https://docs.djangoproject.com/en/dev/ref/settings/#secret-key
SECRET_KEY = get_env_setting('DJANGO_SECRET_KEY')
########## END SECRET CONFIGURATION

CORS_ORIGIN_WHITELIST = ('netavg.dillgroup.io',)
API_AUDIENCE = 'http://netavg.dillgroup.io'

# ========================
# = CELERY CONFIGURATION =
# ========================
BROKER_URL = 'redis://localhost:6379/0'
