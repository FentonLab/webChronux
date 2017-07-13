lcut = 9;
hcut = 13;
attendB = 40;
attenHz = 2;
nyq = 250/2;
Fstop1 = (lcut - attenHz)  / nyq;  
Fpass1 = lcut  / nyq; 
Fpass2 = hcut / nyq;
Fstop2 = (hcut + attenHz) / nyq;
Astop1 = attendB;
Apass  = 1;
Astop2 = attendB;
h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                 Fpass2, Fstop2, Astop1, Apass, Astop2);
Hd = design(h, 'kaiserwin');

%d = fdesign.lowpass('Fp,Fst,Ap,Ast',Fpass1,Fstop1,1,40);

%Hd = design(h, 'kaiserwin');
fvtool(Hd)
%bAlpha = Hd.Numerator;
%disp( bAlpha);