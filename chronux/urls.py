from django.conf.urls import url, include
from django.contrib import admin
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from chronux.views import * 

admin.autodiscover()

urlpatterns = [

    url(r'^admin/', include(admin.site.urls)),
    
    url(r'^$', landing, name = 'landing'),

    url(r'^chronux/$', landing, name = 'landing'),
    url(r'^chronux/processLanding/$', processLanding, name = 'processLanding'),  

    url(r'^chronux/listProjects/$', listProjects, name = 'listProjects'),  
    url(r'^chronux/addProject/$', addProject, name = 'addProject'),  
    url(r'^chronux/submitAddProject/$', submitAddProject, name = 'submitAddProject'),  

    url(r'^chronux/listFiles/$', listFiles, name = 'listFiles'),  
    url(r'^chronux/addDataFile/$', addDataFile, name = 'addDataFile'),  
    url(r'^chronux/listDirectoryFiles/$', listDirectoryFiles, name = 'listDirectoryFiles'),  
    url(r'^chronux/analysisParametersSelect/$', analysisParametersSelect, name = 'analysisParametersSelect'),  
    #url(r'^chronux/displayFileDetails/$', displayFileDetails, name = 'displayFileDetails'),  
    url(r'^chronux/submitAnalysis/$', submitAnalysis, name = 'submitAnalysis'),  
   
    #url(r'^chronux/listSubmittedJobs/$', listSubmittedJobs, name = 'listSubmittedJobs'),  
    
    #url(r'^chronux/addProject/$', addProject, name = 'addProject'),  
    #url(r'^chronux/addDataFile/$', addDataFile, name = 'addDataFile'),  

    #url(r'^chronux/submitAddProject/$', submitAddProject, name = 'submitAddProject'),  
    #url(r'^chronux/submitAddDataFile/$', submitAddDataFile, name = 'submitAddDataFile'),  

    #url(r'^chronux/displayFileDetails/$', displayFileDetails, name = 'displayFileDetails'),  
    #url(r'^chronux/analyzeFileSelect/$', analyzeFileSelect, name = 'analyzeFileSelect'),  

    #url(r'^chronux/analyzeFileSelectFactors/$', analyzeFileSelectFactors, name = 'analyzeFileSelectFactors'), 

    #url(r'^chronux/analyzeFileShowContrastMatrix/$', analyzeFileShowContrastMatrix, name = 'analyzeFileShowContrastMatrix'), 
    
    #url(r'^chronux/analyzeFileSubmit/$', analyzeFileSubmit , name = 'analyzeFileSubmit'),     
    
    #url(r'^chronux/analyzeFileSelectColumns/$', analyzeFileSelectColumns, name = 'analyzeFileSelectColumns'),   
    
    #url(r'^chronux/analyzeProject/$', analyzeProject, name = 'analyzeProject'),       

    ##url(r'^chronux/listAnalysisDetails/$', name = 'listAnalysisDetails'),      
    #url(r'^chronux/displayAnalysisDetail/$', displayAnalysisDetail, name = 'displayAnalysisDetail'),

    #url(r'^chronux/submittedJobDetail/$', submittedJobDetail, name = 'submittedJobDetail'),  
    
    #url(r'^chronux/sampleDetail/$', sampleDetail, name = 'sampleDetail'),     

    #url(r'^chronux/downloadDEGData/$', downloadDEGData, name = 'downloadDEGData'),   

    #url(r'^chronux/downloadData/$', downloadData, name = 'downloadData'),   

    #url(r'^chronux/downloadImage/$', downloadImage, name = 'downloadImage'),
    
    #url(r'^chronux/listAnalyses/$', listAnalyses, name = 'listAnalyses'),   

    #url(r'^chronux/deletePhenotypeFile/$', deletePhenotypeFile, name = 'deletePhenotypeFile'),   
    
    #url(r'^chronux/deleteProjectFiles/$', deleteProjectFiles, name = 'deleteProjectFiles'),   
    
    #url(r'^chronux/deleteProject/$', deleteProject, name = 'deleteProject'),  
    
    #url(r'^chronux/deleteAnalysisDetail/$', deleteAnalysisDetail, name = 'deleteAnalysisDetail'),  
    
    #url(r'^chronux/displayNewProjects/$', displayNewProjects, name = 'displayNewProjects'),   

    #url(r'^chronux/submitAddNewProject/$', submitAddNewProject, name = 'submitAddNewProject'),   

    #url(r'^chronux/showPlots/$', showPlots, name = 'showPlots'), 
    
    #url(r'^chronux/scatterPlot/$', scatterPlot, name = 'scatterPlot'), 

    url(r'^chronux/register/', register, name = 'register'),
    url(r'^accounts/', include('registration.urls'))

]

urlpatterns += staticfiles_urlpatterns()

#print str(urlpatterns)
