import pyedflib
from scipy import *
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import traceback

channelMap = {}

#The first number is the time of the event in seconds since the beginning of the trial.  The second is the duration, which we always define as 0, because we don't have on/off control timing for user annotations.  The user needs to simulate on/off by creating an event like your "Eyes Open" one, and then another like the "Eyes Closed", each of which is treated like a single event of duration 0.  For signal loss situations (where the PC stops receiving from the microEEG), the duration will be non-zero, and will be equal to the time interval that the signal was lost for.  The other two numbers, 0 and 100, never change.  They have some role in the Insight software that I don't know about, but in all of our software, they are always 0 and 100.

class CommentsObj(object):
    def __init__(self):

        self.startTime = 0
        self.stopTime = 0
        self.description = ''

    def __unicode__(self):
        return str(self.project.name)

def parseLayFile(layFileName):

    try:

        f = open(layFileName, "r")
        
        impedanceMapFound = False
        patientDataFound = False
        epochDataFound = False
        commentsFound = False
        startFlag = False
        endFlag = False

        commentsObjList = []
        
        channelMap = {}
        patientMap = {}
        
        for index, line in enumerate(f):
          
            line = line.replace("\n","").replace("\r","").replace("\t","")
            if line.find("ImpedanceMap") != -1:
                impedanceMapFound = True
                continue
            
            if line.find("Patient") != -1:
                impedanceMapFound = False
                patientDataFound = True
                continue      
             
            if line.find("Comments") != -1:
                impedanceMapFound = False
                patientDataFound = False
                commentsFound = True
                continue 
              
            if impedanceMapFound:
                data = line.split("=")
                channelMap[data[0]] = data[1]
                 
            if patientDataFound:
                data = line.split("=")
                patientMap[data[0]] = data[1]
                     
            if commentsFound:
                print ( str(line) ) 
                data = line.split(",")
                print ( str(data))
                if line.find("START") != -1:
                    startFlag = True
                    endFlag = False
                    commentsObj = CommentsObj()
                    commentsObj.startTime = data[0]
                    commentString = data[4]
                    commentsObj.description = commentString[:commentString.find("START")]
                elif line.find("END") != -1:
                    commentsObj.endTime = data[0]
                    print ( " ########## adding ########### ")
                    commentsObjList.append(commentsObj)
        print ( "** " + str([ x.startTime + ":: " + x.endTime  for x in commentsObjList ]) ) 
        print ( channelMap )
        
    except:
        traceback.print_exc(file=sys.stdout)
    return

parseLayFile('test.lay')
print (" end 999 ")
