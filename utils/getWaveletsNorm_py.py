import numpy as np

eegFS = 250
wFactor = 8
bands = {}

def getWaveletsNorm( bands, wFactor, eegFS ):

    bands ["alpha"] = np.linspace(4,7)
    bands ["beta"]  = np.linspace(9,12)
    bands ["low_gamma"] = np.linspace(15,50,5)
    bands ["high_gamma"] = np.linspace(70,90,5)

    wavelets = {} # wavelets in each band    
    #cell(1,length(bands));
    
    for  key, value in bands.items():
        
        f = bands[bandIndex]
        sigmaF = f/wFactor    #practical setting for spectral bandwidth (Tallon-Baudry)
        sigmaT = 1/(sigmaF*pi*2) #wavelet duration
        t = np.linspace ( -4*sigmaT*eegFS , 4*sigmaT*eegFS )
        t = t/eegFS
        #length of t has to be even
        if len(t) % 2 != 0 :
            t = t[:-1]

        S1 = np.exp(    (-1*(np.pow(t,2))) / (2 * ( np.pow(sigmaT,2) ) )    ) # gaussian curve
        S2 = np.exp(2*1j*math.pi*f*t) # sinewave
        
        A = np.pow( (sigmaT * sqrt(pi)) , -0.5 ) #normalization for total power = 1
        psi = A*S1.dot(S2)
        
        #psi = real(psi); %to make output real
        s = sum(abs(real(psi)))/2 #divide by total energy so convolution results in one
        psi = psi / s
        
        wavelets[key] = psi;

    return wavelets
