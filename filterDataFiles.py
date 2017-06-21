import pyedflib
from scipy import *
import numpy as np
import pandas as pd
import sys
import os
import traceback
import shutil
from multitaper import *
import matplotlib.pyplot as plt
from scipy.io import savemat

samplingFrequency = 250
upperFrequency = 70
lowerFrequency = 1
timeBandWidth = 2
numTapers = 3
timeWindow = 5 # time window in seconds
STEP_SIZE = 2 # in seconds
numDataPoints =  timeWindow * samplingFrequency
stepSize = STEP_SIZE * samplingFrequency
padding = 1
paddedNumDataPoints = int ( pow (2, ceil ( math.log( numDataPoints, 2 )) + padding ))

inFilePath = "/Users/smitra/self/andre/MIT-concussion"
outfilePath = "/Users/smitra/projects/webChronux/data"
outAnalysisPath = "/Users/smitra/projects/webChronux/outdata"

BASELINE1 = "baseline_"
BASELINE2 = "baseline"

POST_SEASON = "post-season"

ANALYZE_CHANNELS = ["FP2", "O2", "O1", "FP1", "F7", "F8", "T5", "T6", "T4", "A2", "C3","C4","T3","Cz","Fz"]

individualList = []

patientObjMap = {}

class CommentsObj(object):
    def __init__(self):

        self.commentNum = 0
        self.startTime = ''
        self.description = ''

    def __unicode__(self):
        return str(self.fileName)    

class FileObj(object):
    def __init__(self):

        self.eDFFileName = ''
        self.layFileName = ''

    def __unicode__(self):
        return str(self.eDFFileName)

class PatientObj(object):
    def __init__(self):

        self.patientId = 0

        self.baselineEDFFileName = ''
        self.baselineLayFileName = ''

        self.postSeasonEDFFileName = ''
        self.postSeasonLayFileName = ''
        
        self.fileObjList = []

    def __unicode__(self):
        return str(self.baselineEDFFileName)
    
def getGridIndices():

    try:
    
        if upperFrequency and lowerFrequency:
            paddedNumDataPoints = int ( pow ( 2, ceil ( np.log2 ( numDataPoints ) ) + padding ) )
    
            frequencyResolution = float ( samplingFrequency ) / float ( paddedNumDataPoints )
            gridValues = np.arange ( 0, paddedNumDataPoints , frequencyResolution )
            gridIndices = np.where ( (gridValues >= lowerFrequency ) & (gridValues <= upperFrequency ) )
    
            upperFrequencyGrid = 0
            lowerFrequencyGrid = 0
            # print " gridIndices " + str(gridIndices)
    
            if len( gridIndices ) > 0:
    
                upperFrequencyGrid = gridIndices[0] [len(gridIndices[0]) -1]
    
                if int(gridIndices [0][0]) > 0:
    
                    lowerFrequencyGrid = int(gridIndices [0][0])

    except:
        traceback.print_exc(file=sys.stdout)
        
    #print ( " grid indices = " + str(gridIndices))
    
    return int(lowerFrequencyGrid), int(upperFrequencyGrid) , gridIndices
    
def analyzeData(patientId, fileType, filePath):

    try:

        print ( " in analyze data for patient = " + str(patientId) + " path = " + str(filePath))
        
        patientDir = outAnalysisPath + "/" + str(patientId)

        if not os.path.isdir(patientDir):

            os.mkdir(patientDir )          
        
        f = pyedflib.EdfReader(inFilePath + "/" + filePath)
        n = f.signals_in_file
      
        signal_labels = f.getSignalLabels()
      
        #print ( " num signals = " + str(n) )
        #print ( " signal labels = " + str(signal_labels) )
        
        beginWin = 0
        endWin = 0
      
        numTapers = 2 * timeBandWidth -1
      
        [tapers, eigenValues] = dpss_windows(int(numDataPoints), float(timeBandWidth), int(numTapers) )
      
        spectrumChannelSumData = [0] * ( upperFrequency - lowerFrequency + 1 )
      
        numTapers = len(tapers)
      
        lowerFrequencyGrid, upperFrequencyGrid, gridIndices = getGridIndices ()
        
        channelDataList = []
      
        for channelIndex in range(n):
            
            if signal_labels[channelIndex] not in ANALYZE_CHANNELS:
                
                continue
            
            print ( " channel is = " + str(signal_labels[channelIndex]) ) 

            spectrogramData = []
        
            channelData = f.readSignal(channelIndex)
        
            channelDataList.append(list(channelData))

            print ( " !!!!!! shape is = " + str(np.shape(channelDataList)) )  

            #channelDataList.reshape()
            
            savemat("outmat_" + filePath + ".mat", mdict={'arr': channelDataList} )
    
            #numWindows = round ( len(channelData) / numDataPoints ) + 1
        
            numWindows = int ( ( len ( channelData ) - numDataPoints + 1) / ( stepSize  ) )
        
            #print ( " numWindows = " + str(numWindows) )
        
            for windowNum in range ( numWindows ) :
        
                beginWin = windowNum * numDataPoints * STEP_SIZE
                endWin = beginWin + numDataPoints
        
                #print ( " beginWin = " + str(beginWin) + " endWin = " + str(endWin))
        
                windowData = channelData [ beginWin : endWin]
        
                #print ( " windowData = " + str(windowData) )
        
                if len(windowData) == 0:
        
                    break
        
                for taperIndex, taper in enumerate ( tapers ) :
          
                    taperData = [float(a*b) for a,b in zip(windowData,taper)]
                    
                    fftData = fft(taperData,paddedNumDataPoints)
          
                    spectrumChannelData = np.array([log(abs(x*conj(x))) for x in fftData])
          
                    #print ( " padded num = " + str(paddedNumDataPoints) + " spectrum len = " + str(len(spectrumChannelData)))
          
                    spectrumChannelData = list(spectrumChannelData[gridIndices])
          
                    spectrumChannelData = (1 / float(samplingFrequency) ) * np.array(spectrumChannelData)
                    spectrumChannelData = spectrumChannelData[lowerFrequency:upperFrequency+1]
        
                    #print (spectrumChannelData)
          
                    #print (" for taper = " + str(taperIndex ) + "  spectrumChannelData = " + str(spectrumChannelData) )
                    spectrumChannelSumData = spectrumChannelSumData + array(spectrumChannelData)
        
                spectrumChannelAvgData = np.array( spectrumChannelSumData ) / numTapers
        
                spectrogramData.append(spectrumChannelAvgData)
        
                #print (spectrogramData)
        
                #print (" for window = " + str(windowNum) + " spectrogramData = " + str(spectrogramData) )
        
                #break
        
            np.savetxt(patientDir + "/" + str(fileType) + " channel_spectrogram_data" + str(signal_labels[channelIndex]) + ".txt", spectrogramData )
            print (spectrogramData)
        
            #np.savetxt("outdata/" + str(patientId) + " channel_spectrogram_data" + str(channelIndex) + ".txt", spectrogramData )
            spectrumPSD = [float(sum(col))/len(col) for col in zip(*spectrogramData)]
            spectrumPSD = np.array(spectrumPSD)/100

            #np.savetxt("outdata/channel_spectrum_PSD_data" + str(channelIndex) + ".txt", spectrumPSD )
        
            plt.plot(spectrumPSD.transpose())
            plt.savefig(patientDir + "/" + str(fileType) + "_channel_spectrum_psd_" + str(signal_labels[channelIndex]) + ".png" )
        
            plt.clf()
        
            plt.imshow(np.array(spectrogramData).transpose())
            plt.savefig(patientDir + "/" + str(fileType) + "_channel_spectrogram_" + str(signal_labels[channelIndex]) + ".png" )
        
            plt.clf()
            
       
    except:
        traceback.print_exc(file=sys.stdout)
        return
            

def moveFiles():

    try:
        
        fileList = os.listdir(inFilePath)
        
        if os.path.isdir(outfilePath):
            shutil.rmtree(outfilePath)
            
        if os.path.isdir(outAnalysisPath):
            shutil.rmtree(outAnalysisPath)

        os.mkdir(outfilePath)
        os.mkdir(outAnalysisPath)

        for infileName in fileList:
            
            if ( infileName.find(BASELINE1) == -1 and infileName.find(BASELINE2) == -1 ) or infileName.find(".lay") == -1:
                continue

           # print(infileName)

            channelMap = {}
            commentsObjList = []

            channelMap, commentsObjList = parseLayFile(inFilePath + "/" + infileName )
            if len(commentsObjList) == 0: # skip if not comments
                continue
        
            if infileName.find(BASELINE1) != -1:      
        
                patientId =  infileName[9:9 + infileName[9:].find("_")]
                
            else:
        
                patientId =  infileName[8:8 + infileName[8:].find("_")]
    
            #print (" patientId = " + str(patientId) + " infileName = " + str(infileName))

            if patientId not in patientObjMap:
                
                patientObj = PatientObj()
                
                patientObj.baselineLayFileName = infileName
                
                parsedBaseLineFileName = infileName[:infileName.find(".lay")]
                
                #print ()
                
                patientObj.baselineEDFFileName = parsedBaseLineFileName + ".edf"
                
                # 257894_-post-season_20161201_133919.edf
                
                # baseline_257894_20160930_130415.lay
                
                patientFileNameStart = str(patientId)
                
                otherFiles = [x for x in fileList if x.startswith (str(patientId)) ]

                postSeasonFiles = [x for x in otherFiles if x.find("post-season") != -1]

                nonPostSeasonFiles = [x for x in otherFiles if x.find("post-season") == -1]

                #print ( " post season files = " + str(postSeasonFiles) ) 
                #print ( " non post season files = " + str(nonPostSeasonFiles) ) 

                if len(postSeasonFiles) > 0:
                    
                    for postSeasonFileName in postSeasonFiles:
                    
                        if postSeasonFileName.find(".lay") != -1: 
                    
                            patientObj.postSeasonLayFileName = postSeasonFileName
                            
                            parsedPostSeasonFileName = postSeasonFileName[:postSeasonFileName.find(".lay")]
                            
                            #print (" postSeasonFileName " + str(postSeasonFileName) + " parsedPostSeasonFileName = " + str(parsedPostSeasonFileName))
                            
                            patientObj.postSeasonEDFFileName = parsedPostSeasonFileName + ".edf"
                            
                            break
                
                patientObjMap[patientId] = patientObj
                #print (" adding " + str(patientId) + " baseline " + patientObj.baselineEDFFileName +                        " post season " + str(patientObj.postSeasonEDFFileName))
                
                #break

        for patientId, patientObj in patientObjMap.items():
            
            patientDir = outfilePath + "/" + str(patientId)

            print ( " in = " + patientObj.postSeasonEDFFileName.strip() + " out = " + patientObj.postSeasonEDFFileName.strip() )
        
            #if patientObj.postSeasonEDFFileName != '' or patientObj.postSeasonLayFileName != '':
                #continue
        
            if os.path.isdir(inFilePath + "/" + patientObj.postSeasonEDFFileName) :
                #print ( " 1111 "  +  patientObj.postSeasonEDFFileName) 
                continue
                
            if os.path.isdir(inFilePath + "/" + patientObj.postSeasonLayFileName):
                #print ( " 2222 "  +  patientObj.postSeasonEDFFileName) 
                continue
            #print ( " ADDING " ) 
        
            os.mkdir(patientDir )        

            shutil.copyfile (inFilePath + "/" + patientObj.baselineEDFFileName, patientDir + "/" + patientObj.baselineEDFFileName)
            shutil.copyfile (inFilePath + "/" + patientObj.baselineLayFileName, patientDir + "/" + patientObj.baselineLayFileName)
        
            shutil.copyfile (inFilePath + "/" + patientObj.postSeasonEDFFileName, patientDir + "/" + patientObj.postSeasonEDFFileName)
            shutil.copyfile (inFilePath + "/" + patientObj.postSeasonLayFileName, patientDir + "/" + patientObj.postSeasonLayFileName)
            
            analyzeData(patientId, "baseline" , patientObj.baselineEDFFileName)
            #break
            #analyzeData(patientId, "postSeason", patientObj.postSeasonEDFFileName)
        
        print (patientObjMap.keys())

    except:
        traceback.print_exc(file=sys.stdout)
    return


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
                #print (" ####### data = " + str(data) )
                channelMap[data[0]] = data[1]

            if commentsFound:

                #print ( str(line) )
                data = line.split(",")
                #print ( str(data))

                if line.find("loss") == -1:

                    commentsObj = CommentsObj()
                    commentsObj.commentNum = commentNum
                    commentNum += 1
                    commentsObj.startTime = data[0]
                    if len (data) >= 5:
                        commentString = data[4]
                        commentsObj.description = commentString
                        commentsObjList.append(commentsObj)

       # print ( "** " + str([ x.startTime + ":: " + x.description  for x in commentsObjList ]) )
       # print ( channelMap )

    except:
        traceback.print_exc(file=sys.stdout)
    return channelMap, commentsObjList

moveFiles()
print (" end 999 ")
