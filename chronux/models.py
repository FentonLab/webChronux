from django.db import models
from django.contrib.auth.models import User
    
class Project(models.Model):
    name = models.CharField(max_length=256)
    user = models.ForeignKey(User, blank = True, null=True)    
    description = models.CharField(max_length=512, blank = True, null = True)
    def __str__(self):
        return self.name

class FileType(models.Model):
    name = models.CharField(max_length=256)
    description = models.CharField(max_length=512)
    def __str__(self):
        return self.name

class JobStatusCode(models.Model):
    code = models.CharField(max_length=10)
    description = models.CharField(max_length=255, null=True)
    def __unicode__(self):
        return self.code

class SubmittedJobType(models.Model):
    name = models.CharField ( max_length=10)  
    description = models.CharField ( max_length=255)   
    def __unicode__(self):
        return str(self.description)     

class LeadList(models.Model):
    name = models.CharField ( max_length=50) 
    description = models.CharField ( max_length=255, null = True, blank = True) 
    def __unicode__(self):
        return self.name    
    
class LeadListChannel(models.Model):
    leadList = models.ForeignKey(LeadList)

    channelId = models.CharField(max_length=10) 
    baseChannelId = models.CharField(max_length=10) 
    channelOrder = models.IntegerField(null = True, blank = True)

    xCoordinate = models.FloatField (null = True, blank = True) 
    yCoordinate = models.FloatField (null = True, blank = True)     
    zCoordinate = models.FloatField (null = True, blank = True) 

    isInterior = models.NullBooleanField()
    def __unicode__(self):
        return self.channelId
    
class BipolarMontage(models.Model):
    channelId1 = models.ForeignKey(LeadListChannel, related_name="channelId1")
    channelId2 = models.ForeignKey(LeadListChannel, related_name="channelId2")
    def __unicode__(self):
        return str(self.channelId1 ) + str(self.channelId2 )    

class Datafile(models.Model):
    name = models.CharField(max_length=256)
    description = models.CharField(max_length=512)
    uploadedBy = models.ForeignKey(User, blank = True, null = True )
    uploadedDate = models.DateTimeField(blank = True, null = True )
    filePath = models.CharField(max_length=512)    
    fileType = models.ForeignKey(FileType, blank = True, null = True )
    project = models.ForeignKey ( Project ) 
    leadList = models.ForeignKey(LeadList)
    def __str__(self):
        return self.name

class Epoch(models.Model):
    epochValue1 = models.FloatField()
    epochValue2 = models.FloatField()
    epochValue3 = models.FloatField()
    epochValue4 = models.FloatField()
    epochType = models.CharField(max_length = 256)
    datafile = models.ForeignKey(Datafile)
    def __str__(self):
        return self.epochType
    
class AnalysisDetail(models.Model):
    
    name = models.CharField ( max_length=255, blank = True, null=True)

    project = models.ForeignKey(Project) 
    description = models.CharField ( max_length=255)

    timeWindow = models.FloatField(blank = True, null=True)
    frequencyBandWidth = models.FloatField(blank = True, null=True)
    numTapers = models.IntegerField(blank = True, null=True)
    stepSize = models.FloatField(blank = True, null=True)
    padding = models.IntegerField(blank = True, null=True)

    upperFrequency = models.FloatField(blank = True, null=True)
    lowerFrequency = models.FloatField(blank = True, null=True)

    analysisOutputPath = models.CharField(max_length=512, blank = True, null=True)

    def __unicode__(self):
        return str(self.name)

class AnalysisResultFile(models.Model):
 
    name = models.CharField ( max_length=128)  
    filePath = models.CharField ( max_length=255)
    analysisDetail = models.ForeignKey(AnalysisDetail) 
    resultFileName = models.CharField ( max_length=128)      
    
    def __unicode__(self):
        return str(self.description) 

class AnalysisPlot(models.Model):
 
    name = models.CharField ( max_length=128)  
    plotPath = models.CharField ( max_length=255)
    analysisDetail = models.ForeignKey(AnalysisDetail) 
    plotFileName = models.CharField ( max_length=128)  
    
    def __unicode__(self):
        return str(self.description) 
    
class SubmittedJob(models.Model):
    name = models.CharField(max_length=512)
    description = models.CharField(max_length=512, null=True, blank = True)
    submittedBy = models.ForeignKey(User) 
    submittedOn = models.DateTimeField( null=True, blank = True)

    submittedJobType = models.ForeignKey(SubmittedJobType)
    
    jobStatusCode = models.ForeignKey(JobStatusCode)
    analysisDetail = models.ForeignKey(AnalysisDetail, null=True, blank = True)
    downloadDataFileLink = models.CharField(max_length=512, null=True, blank = True)
    completedTime = models.DateTimeField(null=True, blank=True)
    def __unicode__(self):
        return self.name    
