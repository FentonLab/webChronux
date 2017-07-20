x = 0:30:90;

s1 = sin( x * pi / 180 )

s2 = sin( x * pi / 18 ) % s2 = 10* w1

s2 = s1.*s2 

%print ( s2)
s3 = s1 + s2 

%print ( s2)

%To evaluate the performance of your code.

%You can even add 60 Hz noise to test the hypothesis that the difference in the python and may lab examples is due effects of 60Hz:

%S3 +sin(60*x)
