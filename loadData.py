from chronux.models import *
import os, sys, traceback
import pandas as pd

#QC_FILE_BASE_PATH = "/static/qc/"

jobStatusCode = JobStatusCode(code= "QUEUE", description = "In Queue")
jobStatusCode.save()

jobStatusCode = JobStatusCode(code= "START", description = "Started")
jobStatusCode.save()

jobStatusCode = JobStatusCode(code= "END", description = "Ended")
jobStatusCode.save()

jobStatusCode = JobStatusCode(code= "TERM", description = "Terminated")
jobStatusCode.save()

submittedJobType = SubmittedJobType(name = "Download", description = "Create download Zip File")
submittedJobType.save()

submittedJobType = SubmittedJobType(name = "Upload", description = "Upload files")
submittedJobType.save()

submittedJobType = SubmittedJobType(name = "Analysis", description = "Analysis of data")
submittedJobType.save()

dataFileType = FileType(name = "layFile", description = "ContrastMatrix")
dataFileType.save()

dataFileType = FileType(name = "edfFile", description = "Zip of EDF files")
dataFileType.save()
