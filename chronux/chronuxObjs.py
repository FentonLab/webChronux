class ProjectObj(object):
    def __init__(self):

        self.project = ''
        self.numContrastMatrixFiles = 0

    def __unicode__(self):
        return str(self.project.name)
    
class EDFFileObj(object):
    def __init__(self):

        self.fileName = ''
        self.channelMap = ''
        self.commentsObjList = []

    def __unicode__(self):
        return str(self.fileName)    

class CommentsObj(object):
    def __init__(self):

        self.commentNum = 0
        self.startTime = ''
        self.description = ''

    def __unicode__(self):
        return str(self.fileName)    

class AnalysisObj(object):
    def __init__(self):

        self.filePaths = []
        self.spectrogramChannels = []
        
        self.samplingFrequency = 0
        self.removeLineNoise = False
    
        self.timeWindow = 0
        self.numDataPoints = 0
    
        self.bandWidth = 0
    
        self.stepSize = 0
        self.padding = 0
        
        self.upperFrequency = 0
        self.lowerFrequency = 0
        
        self.numTapers = 0
        
    def __unicode__(self):
        return str(self.fileName)    

class BipolarObj(object):
    def __init__(self):

        self.channel1 = ''
        self.channel2 = ''

    def __unicode__(self):
        return str(self.channel1) + "-" + str(self.channel2)  