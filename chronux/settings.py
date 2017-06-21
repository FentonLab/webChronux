"""
Django settings for chronux project.

For more information on this file, see
https://docs.djangoproject.com/en/1.7/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.7/ref/settings/
"""

import djcelery
djcelery.setup_loader()
import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
#import pyexcel.ext.xlsx  # This is required for tmp file handlers.

#import pyexcel.ext  # This is required for tmp file handlers.

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

PROJECT_BASE_FOLDER = BASE_DIR 

DATA_OUTPUT_FOLDER = BASE_DIR + "/data"

DOWNLOAD_FOLDER = BASE_DIR + "/download"

IMAGE_OUTPUT_FOLDER = BASE_DIR + "/webchronux/static/img"

INPUT_DATA_FOLDER = BASE_DIR + "/indata"

OUTPUT_DATA_FOLDER = BASE_DIR + "/outdata"

ROOT_URLCONF = 'chronux.urls'

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.7/howto/deployment/checklist/

## Make this unique, and don't share it with anybody.
SECRET_KEY = 'idmm)d5h)kwxcopo%0rirfk@^88xcki4dx#slgp5_v#0)lt7gd'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

#TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []

# Application definition
INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'chronux',
    'djcelery',
    'registration'
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.auth.middleware.SessionAuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

# table prefix name
DB_PREFIX = 'chronux'
# Local time zone for this installation. Choices can be found here:
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'chronuxdb',                      # Or path to database file if using sqlite3.
        'USER': 'root',                      # Not used with sqlite3.
        'PASSWORD': 'admin',                  # Not used with sqlite3.
        'HOST': 'localhost',                      # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '3306',                      # Set to empty string for default. Not used with sqlite3.
        'OPTIONS': {
            'init_command': 'SET default_storage_engine=INNODB',
            }
    }
}

BROKER_URL = 'amqp://guest:guest@localhost:5672//'

# Internationalization
# https://docs.djangoproject.com/en/1.7/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

STATICFILES_DIRS = [
     os.path.join(BASE_DIR, './static'),
]

STATIC_URL = '/static/'

# TEMPLATE_DIRS = (
#     BASE_DIR + '/templates/',
# )

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                #'django.template.context_processors.debug',
                #'django.template.context_processors.request',
                #'django.core.context_processors.csrf'                
                'django.contrib.auth.context_processors.auth',
                #'django.contrib.messages.context_processors.messages',
            ], 
            },
    },
]
#FILE_UPLOAD_HANDLERS = ("django_excel.ExcelMemoryFileUploadHandler",
                        #"django_excel.TemporaryExcelFileUploadHandler",
                        #"django.core.files.uploadhandler.MemoryFileUploadHandler",
                        #"django.core.files.uploadhandler.TemporaryFileUploadHandler",)

FILE_UPLOAD_HANDLERS = ("django.core.files.uploadhandler.MemoryFileUploadHandler",
                        "django.core.files.uploadhandler.TemporaryFileUploadHandler",)


LOGIN_REDIRECT_URL = "/chronux/"
