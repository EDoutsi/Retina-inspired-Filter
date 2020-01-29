%%%%%%%%%%%%%%%%%%%%%%   BUILD  THE RIF FILTER     %%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% INPUT:  1. sc, ss -> the standard deviations of the center and surround Gaussian filters, respectively, 
%         2. x -> the linespace, 
%         3. tsamp -> the sampling rate of the time vector,
%         4. tauC, tauS and tauG -> constant parameters related to the temporal filters, 
%         5. wC and wS -> two weights 
%         6. T -> the size of the time window 
%         7. tmax -> a time parameter which is used for illustration purpose
% -----------------------------------------------------------------------
% OUTPUT: 1. Filter -> The RIF filter in space domain
%         2. fftFilter -> the RIF filter in frequency domain
% -----------------------------------------------------------------------
% SET OF BIO-PLAUSIBLE PARAMETERS 
%       sc  = 0.5; ss  = 3*sc;
%       max = 10; min  = -max; samp = 3*max+1; 
%       tauC = 20.*10^-3; tauS = 4.*10^-3; tauG = 5.*10^-3;
%       MaxIter = 10000;
%       wC = 1; wS = 1;                 % weight surround
%       T  = 150; tmax = T; tsamp = 1;  % number od time instances
%       t1 = 1:tsamp:tmax;              % temporal vector
%       x  = linspace(min,max,samp);    % spatial vector
% -----------------------------------------------------------------------
function [Filter,fftFilter]=...
    RIFFilter(sc,ss,x,tsamp,tauC,tauS,tauG,wC,wS,T,tmax)
    %% -- 2D GUASSIAN FILTERS
    [GausC, FGausC] = GaussianKernel(x,sc); % in tests 11
    [GausS, FGausS] = GaussianKernel(x,ss); % in tests 11
    [Gx,Gy] = size(GausC);
    [fGx , fGy] = freqspace([Gx Gy]);
    %% -- 1D TEMPORAL FILTERS
    R_C = ComputingRc(T, tmax, wC, tauC, tauG);
    R_S = ComputingRs(T, tmax, wS, tauC, tauG, tauS);
    tmax2 = 3*T;
    t2 = 1:tmax2;
    R_C_2T = ComputingRc(T, tmax2, wC, tauC, tauG);
    R_S_2T = ComputingRs(T, tmax2, wS, tauC, tauG, tauS);
    %% -- 2D RIF Filter 
    for j = 1:tsamp:T
        Filter(j,:)  = wC * GausC(:) .* R_C(j) - wS * GausS(:) .* R_S(j);       % space
        fftFilter(j,:)  = wC * FGausC(:) .* R_C(j) - wS * FGausS(:) .* R_S(j);  % frequency
    end
    [FMx,FMy] = size(Filter);
end
