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


