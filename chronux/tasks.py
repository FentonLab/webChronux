from __future__ import absolute_import

import csv
import numpy
from numpy import *
import scipy
import math
from scipy import *
from math import *

import sys, traceback
from chronux.models import *

import subprocess, datetime, time
from subprocess import Popen

import re, os

import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr

import sys, traceback

import pymongo
from pymongo import MongoClient

from celery import Celery
from celery import shared_task

#from tasks import add

from chronux.chronuxObjs import *
from chronux.chronuxConstants import *

#import zipfile
#import datetime

from django.conf import settings

#app = Celery('tasks', broker='amqp://guest@localhost//')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'chronux.settings')

app = Celery('chronux')

@app.task
#@shared.task
def submitLimma(filePath, basePlotPath, outputPath, factors, factorColumns, factorColumnValuesMap, factorColumnValueSet,samples, sampleNames, animals, contrasts, factorSetMap, submittedJob, projectId, userId, analysisDetailId, startingColumn, commaOrTabDelimitedFlag, allContrasts):

    contrast = []

    try:

        #print " in task - before calling R !!!! "

        startJobStatusCode = JobStatusCode.objects.all().filter(code="START")[0]
        endJobStatusCode = JobStatusCode.objects.all().filter(code="END")[0]

        submittedJob.jobStatusCode = startJobStatusCode

        submittedJob.save()

        MAT_FN = robjects.r("source('" + settings.PROJECT_BASE_FOLDER + "/limma_analysis.R'); contrasts_matrix")

        contrast = MAT_FN(filePath, basePlotPath, outputPath, robjects.StrVector(factors),  robjects.StrVector(factorColumns),robjects.ListVector(factorColumnValuesMap), robjects.StrVector(factorColumnValueSet),robjects.StrVector(samples), robjects.StrVector(sampleNames), robjects.StrVector(animals), robjects.StrVector(contrasts), robjects.ListVector(factorSetMap), str(projectId), userId, str(analysisDetailId), startingColumn, commaOrTabDelimitedFlag, allContrasts )

        submittedJob.completedTime = datetime.datetime.now()
        submittedJob.jobStatusCode = endJobStatusCode
        submittedJob.save()

        #print " &&&&&&&& in task - after calling R "

    except:

        traceback.print_exc(file=sys.stdout)

    return