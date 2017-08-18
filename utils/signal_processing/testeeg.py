import numpy as np
from math import *
x = np.linspace(0,200,5000);

s1 = np.sin( x * pi / 180 )

s2 = np.sin( x * pi / 18 )

s2 = np.array ( [x*y for x,y in zip(s1,s2) ] ) 

s3 = s1 + s2 

eeg = s3

print (eeg)
