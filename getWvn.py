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

        S1 = np.exp(    (-1*(t**)) / (2*(sigmaT**2))    ) # gaussian curve
        S2 = np.exp(2*1i*pi*f*t) # sinewave
        
        A = (sigmaT * sqrt(pi))**(-0.5) #normalization for total power = 1
        psi = A*S1.dot(S2)
        
        #psi = real(psi); %to make output real
        s = sum(abs(real(psi)))/2 #divide by total energy so convolution results in one
        psi = psi / s
        
        wavelets[key] = psi;

    return wavelets

wavelets = getWaveletsNorm( bands, wFactor, eegFS );

#filter
for filtI in range( length(wavelets) )

        psi = wavelets{filtI};

        %convolution of wavelet and signal
        c = conv(eeg,psi);
        %fix start and end
        N = round((length(psi)-1)/2);
        c = c(N:length(c)-N);
        if length(c) > size(eeg,2)
            c = c(1:size(eeg,2));
        end
        power = (abs(c)).^2;

        powers(filtI) = mean(power)


%create Morlet Wavelets
%normalizes wavelet for the total energy so it creates same output as FIR
%filter
function wavelets = getWaveletsNorm( bands, wFactor, eegFS )

wavelets = cell(1,length(bands));

for bI = 1:length(bands)
    f = bands(bI);
    sigmaF = f/wFactor;    %practical setting for spectral bandwidth (Tallon-Baudry)
    sigmaT = 1/(sigmaF*pi*2); %wavelet duration
    t = -4*sigmaT*eegFS:4*sigmaT*eegFS;
    t = t/eegFS;
    %length of t has to be even
    if rem(length(t),2) == 0
        t = t(1:end-1);
    end
    S1 = exp((-1*(t.^2)) / (2*sigmaT^2)); %gaussian curve
    S2 = exp(2*1i*pi*f*t); %sinewave
    A = (sigmaT * sqrt(pi))^(-0.5); %normalization for total power = 1
    psi = A*S1.*S2;
    
    %psi = real(psi); %to make output real
    s = sum(abs(real(psi)))/2; %divide by total energy so convolution results in one
    psi = psi / s;
    
    wavelets{bI} = psi;
end
