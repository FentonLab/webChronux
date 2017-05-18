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
