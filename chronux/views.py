from django.contrib.auth.decorators import login_required

from django.contrib import messages

from django.template import RequestContext
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.http import HttpResponse

from chronux.models import *
from chronux.chronuxObjs import *
from chronux.chronuxConstants import *

import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr

import pandas as pd
from pandas import DataFrame

import chronux.tasks as ts
import seaborn as sns
from seaborn import color_palette, diverging_palette

from bokeh.plotting import *
import bokeh
from bokeh.charts import Scatter, output_file, show
from bokeh.sampledata.autompg import autompg as df
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import HoverTool

from bokeh.embed import components
from bokeh.models import Range1d

from bokeh.plotting import figure

from collections import OrderedDict

from bokeh.charts import Bar, output_file
#from django.core.context_processors import csrf
import os, sys, traceback

import zipfile
import hashlib
import zlib

import io
import re
import shutil

from django.conf import settings
import numpy
import itertools
import csv
import datetime
import os.path
from chronux.utils import * 
#import chronux.task_utils as tsk

from bokeh.charts import Bar, output_file, show
from bokeh.charts.attributes import cat, color
from bokeh.charts.operations import blend

from multiprocessing import Pool

import math
from math import *

def register(request):
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            new_user = form.save()
            return HttpResponseRedirect("/registration/")
    else:
        form = UserCreationForm()
    return render(request, "registration/register.html", {
        'form': form,
    })

#@login_required
def landing(request):
    return render(request, "chronux/landing.html", {

    })

#@login_required
#@csrf_protect
def processLanding(request):

    ''' Function process landing page submit
    Input: params.request
    Output: Delegates task to corresponding handler
    '''
    chronuxHomeButton = request.POST.get("chronuxHomeButton","1" )

    if chronuxHomeButton == "0":
        return submittedJobs ( request )
    elif chronuxHomeButton == "1" :
        return listProjects( request )
    
@login_required
def listProjects(request):

    projectObjList = []

    if request.user.is_superuser:
        projects = Project.objects.all()
    else:
        projects = Project.objects.filter(user = request.user)

    edfFileType = FileType.objects.filter(name = "edfFile")[0]

    for project in projects:

        projectObj = ProjectObj()

        projectObj.project = project

        dataFiles = Datafile.objects.filter ( project = project, fileType = edfFileType )

        projectObjList.append(projectObj)

    return render(request, 'chronux/listProjects.html', {
        "projectObjList":projectObjList,
    }, RequestContext(request))

@login_required
def addProject(request):

    try:

        projects = Project.objects.filter ( user = request.user)

    except:
        traceback.print_exc(file=sys.stdout)
        #messages.add_message(request, messages.ERROR, 'Error occurred while fetching details for sample detail id' + str(sampleDetailId) )

    return render(request, 'chronux/addProject.html', {

        "projects" : projects ,

    })

@login_required
def addDataFile(request):

    try:

        projectId = request.POST.get("projectId",0 )

        project = Project.objects.get(pk = projectId)

    except:
        traceback.print_exc(file=sys.stdout)

    return render(request, 'chronux/addDataFile.html', {

        "project" : project,

    })

@login_required
def submitAddProject(request):

    messages = []

    try:

        projectName = request.POST.get("projectName","" )
        projectDescription = request.POST.get("projectDescription","" )

        project = Project ( name = projectName, user = request.user, description = projectDescription )
        project.save()

        edfFileType = FileType.objects.filter(name = "edfFile")[0]
        dataFiles = Datafile.objects.filter ( project = project, fileType = edfFileType )

    except:
        traceback.print_exc(file=sys.stdout)

    return render(request, 'chronux/listFiles.html', {
            "dataFiles":dataFiles,
            "project":project,
        })

@login_required
def listDirectoryFiles(request):

    project = ""

    try:

        projectId = request.POST.get("projectId",0 )

        project = Project.objects.get (pk = projectId)

        edfFileType = FileType.objects.filter(name = "edfFile")[0]

        layFileType = FileType.objects.filter(name = "layFile")[0]

        # load qc files if available
        
        datafileDirectory = ''

        try:

            datafileDirectory = request.POST.get("dataFileDirectory","")
            fileList = os.listdir(datafileDirectory)
            print (fileList)            
            fileList = [x for x in fileList if x.find(".edf") != -1]
            
            print (fileList)

        except:
    
            traceback.print_exc(file=sys.stdout)

        edfFileType = FileType.objects.filter(name = "edfFile")[0]

        dataFiles = Datafile.objects.filter ( project = project, fileType = edfFileType )

    except:
        traceback.print_exc(file=sys.stdout)

    return render(request, 'chronux/listDirectoryFiles.html', {
        "project":project,
        "fileList":fileList,
        "datafileDirectory" : datafileDirectory,
    })

@login_required
def listFiles(request):

    projectId = request.POST.get("projectId",0 )

    #print " project = " + str(projectId)

    project = Project.objects.get(pk = projectId)

    edfFileType = FileType.objects.filter(name = "edfFile")[0]

    dataFiles = Datafile.objects.filter ( project = project, fileType = edfFileType )

    return render(request, 'chronux/listFiles.html', {
        "dataFiles":dataFiles,
        "project":project,
    })

@login_required
def analysisParametersSelect(request):

    projectId = request.POST.get("projectId",0 )
    
    dataFileNames = request.POST.getlist("fileName")

    #print " project = " + str(projectId)

    project = Project.objects.get(pk = projectId)

    edfFileType = FileType.objects.filter(name = "edfFile")[0]
    
    datafileDirectory = request.POST.get("datafileDirectory","" )
    
    fileList = os.listdir(datafileDirectory)

   # print (fileList)            
    fileList = [x for x in fileList if x in dataFileNames] 
    
    dataFileList = []
    
    for fileName in dataFileNames:
        filePath = datafileDirectory+ "/" + fileName
        datafiles = Datafile.objects.filter(filePath = filePath)
        if len(datafiles) == 0:
            datafile = Datafile ( filePath = datafileDirectory+ "/" + fileName, fileType = edfFileType, project = project )
            datafile.save()
        else:
            datafile = datafiles[0]
        dataFileList.append(datafile)
        
    firstFileName = fileList[0]
    
    layFileName = firstFileName[:firstFileName.index(".edf")] + ".lay"
    
    channelMap, commentsObjList = parseLayFile(datafileDirectory + "/" + layFileName )
    
    edfFileObjList = []    
    
    for fileName in fileList:
        
        layFileName = fileName[:fileName.index(".edf")] + ".lay"
        
        channelMap1, commentsObjList1 = parseLayFile(datafileDirectory + "/" + layFileName )
        
        edfFileObj = EDFFileObj () 
        
        edfFileObj.fileName = fileName
        edfFileObj.channelMap = channelMap1
        edfFileObj.commentsObjList = commentsObjList
        
        edfFileObjList.append(edfFileObj)
        
    bipolarMappingFilePath = settings.PROJECT_BASE_FOLDER + "/bipolar_mapping.txt"
    
    bipolarMappingFile = open(bipolarMappingFilePath, "r")
    
    print (" bipolarMappingFilePath = " + str(bipolarMappingFilePath))

    bipolarObjList = []    
    
    for line in bipolarMappingFile:
        line = line.replace("\n","").replace("\t","")
        if line == "":
            break
        data = line.split("-")
        
        bipolarObj = BipolarObj()
        
        bipolarObj.channel1 = data[0]
        bipolarObj.channel2 = data[1]
        
        bipolarObjList.append(bipolarObj)
    
    print (str(bipolarObjList))
        
    return render(request, 'chronux/analysisParametersSelect.html', {
        "project":project,
        "channelMap":channelMap,
        "edfFileObjList":edfFileObjList,
        "dataFileList" : dataFileList,
        "bipolarObjList": bipolarObjList,
        "eegBandMap": EEG_BANDS,
    })

@login_required
def submitAnalysis(request):

    projectId = request.POST.get("projectId",0 )
    
    project = Project.objects.get ( pk = projectId )
    datafiles = Datafile.objects.filter( project = project ) 

    spectrogramChannels = request.POST.getlist("spectrogramChannel") 
    
    samplingFrequency = float(request.POST.get('samplingFrequency', 0)) 
    removeLineNoise = request.POST.get('removeLineNoise', False)

    timeWindow = float(request.POST.get('timeWindow', '0') )

    bandWidth = float ( request.POST.get('bandWidth', '') )

    stepSize = float ( request.POST.get('stepSize', '0') )
    padding = float ( request.POST.get('padding', '0') )
    
    upperFrequency = int ( request.POST.get('upperFrequency', '0') )
    lowerFrequency = int ( request.POST.get('lowerFrequency', '0') )
    
    numTapers = int ( request.POST.get('numTapers', '0') )    
    
    coherogramChannel = request.POST.getlist('coherogramChannel')

    coherogramChanneName1 = request.POST.getlist('coherogramChanneName1')
    coherogramChanneName2 = request.POST.getlist('coherogramChanneName2')
    
    analysisObj = AnalysisObj()
        
    analysisObj.datafiles = datafiles
    analysisObj.spectrogramChannels = spectrogramChannels
    
    analysisObj.samplingFrequency = samplingFrequency
    analysisObj.removeLineNoise = removeLineNoise

    analysisObj.timeWindow = timeWindow

    analysisObj.bandWidth = bandWidth

    analysisObj.stepSize = stepSize
    analysisObj.padding = padding
    
    analysisObj.upperFrequency = upperFrequency 
    analysisObj.lowerFrequency = lowerFrequency

    print ( " analysisObj.datafiles = " + str(datafiles))
    print ( " analysisObj.spectrogramChannels = " + str(analysisObj.spectrogramChannels))
    print ( " analysisObj.samplingFrequency = " + str(analysisObj.samplingFrequency))
    print ( " analysisObj.removeLineNoise = " + str(analysisObj.removeLineNoise))
    print ( " analysisObj.timeWindow = " + str(analysisObj.timeWindow))
    print ( " analysisObj.numDataPoints = " + str(analysisObj.numDataPoints))
    print ( " analysisObj.bandWidth = " + str(analysisObj.bandWidth))
    print ( " analysisObj.stepSize = " + str(analysisObj.stepSize))
    print ( " analysisObj.padding = " + str(analysisObj.padding))
    print ( " analysisObj.upperFrequency = " + str(analysisObj.upperFrequency))

    numDataPoints = int ( math.ceil ( timeWindow * samplingFrequency ) )     
    
    analysisObj.numDataPoints = numDataPoints    
    
    timeBandWidth = timeWindow * bandWidth
    numTapersCalculated = 2 * timeBandWidth - 1
    
    if numTapers == 0:
        
        numTapersCalculated = numTapers    
        
    analysisObj.numTapers = numTapers
    
    analysisResultsObj = analyzeEDFData(analysisObj)
    
    #ts.analyzeEDFData.delay(analysisObj)   
    return render(request, 'chronux/analysisResults.html', {
        "project":project,
    })

def parseLayFile(layFileName):

    try:
        
        print (layFileName)

        f = open(layFileName, "r")
        
        impedanceMapFound = False
        chennalMapFound = False
        patientDataFound = False
        epochDataFound = False
        commentsFound = False
        startFlag = False
        endFlag = False

        commentsObjList = []
        
        channelMap = {}
        patientMap = {}
        
        commentNum = 0
        
        for index, line in enumerate(f):
          
            line = line.replace("\n","").replace("\r","").replace("\t","")
            
            if line.find("ImpedanceMap") != -1:
                impedanceMapFound = True
                continue
            
            if line.find("ChannelMap") != -1:
                impedanceMapFound = False
                channelMapFound = True
                continue            
            
            if line.find("Patient") != -1:
                channelMapFound = False
                patientDataFound = True
                continue      
             
            if line.find("Comments") != -1:
                impedanceMapFound = False
                patientDataFound = False
                channelMapFound = False
                commentsFound = True
                continue 
              
            if impedanceMapFound:
                data = line.split("=")
                print (" ####### data = " + str(data) )
                channelMap[data[0]] = data[1]
                 
            if patientDataFound:
                data = line.split("=")
                #patientMap[data[0]] = data[1]
                     
            #if commentsFound:
                #print ( str(line) ) 
                #data = line.split(",")
                #print ( str(data))
                #if line.find("START") != -1:
                    #startFlag = True
                    #endFlag = False
                    #commentsObj = CommentsObj()
                    #commentsObj.startTime = data[0]
                    #commentString = data[4]
                    #commentsObj.description = commentString[:commentString.find("START")]
                #elif line.find("END") != -1:
                    #commentsObj.endTime = data[0]
                    #print ( " ########## adding ########### ")
                    #commentsObjList.append(commentsObj)

            if commentsFound:
                
                print ( str(line) ) 
                data = line.split(",")
                print ( str(data))
                
                if line.find("loss") == -1:
                
                    commentsObj = CommentsObj()
                    commentsObj.commentNum = commentNum
                    commentNum += 1
                    commentsObj.startTime = data[0]
                    commentString = data[4]
                    commentsObj.description = commentString
                    commentsObjList.append(commentsObj)
                    
        print ( "** " + str([ x.startTime + ":: " + x.description  for x in commentsObjList ]) ) 
        print ( channelMap )
        
    except:
        traceback.print_exc(file=sys.stdout)
    return channelMap, commentsObjList 


#@login_required
#def addDataFile(request):

    #try:

        #projectId = request.POST.get("projectId",0 )

        #project = Project.objects.get(pk = projectId)

        #print (" starting = " + str(startingColumn) )

    #except:
        #traceback.print_exc(file=sys.stdout)

    #return render(request, 'chronux/addDataFile.html', {

        #"project" : project,

    #})



    

