import re, sys, traceback, math, numpy , scipy 
from chronuxplus.models import *
from numpy import *
from scipy import *
# import pdb
# pdb.set_trace()

def loadShellDetails():
    shell = Shell (name = "Shell 1", radius = 1.0, conductivity = 1, sequence = 1)
    shell.save()
    shell = Shell (name = "Shell 2", radius = 0.9467, conductivity = 0.0125, sequence = 2)
    shell.save()
    shell = Shell (name = "Shell 3", radius = 0.8667, conductivity = 3.0, sequence = 3)
    shell.save()
    shell = Shell (name = "Shell 4", radius = 0.84, conductivity = 1.0, sequence = 4)
    shell.save()
shell = Shell (name = "Shell 1", radius = 1.0, conductivity = 1, sequence = 1)
shell.save()
shell = Shell (name = "Shell 2", radius = 0.9467, conductivity = 0.0125, sequence = 2)
shell.save()
shell = Shell (name = "Shell 3", radius = 0.8667, conductivity = 3.0, sequence = 3)
shell.save()
shell = Shell (name = "Shell 4", radius = 0.84, conductivity = 1.0, sequence = 4)
