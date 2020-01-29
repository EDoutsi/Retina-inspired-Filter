%%%%%%%%%%%%%%%%%%   COMPONENTS OF THE RIF FILTER     %%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% INPUT --> sc and ss are the standard deviations of the center and surround 
%           Gaussian filters respectively, x is the linespace, tsamp is the 
%           sampling rate of the time vector, tauC, tauS and tauG are constant 
%           parameters related to the temporal filters, wC and wS are two weights, 
%           T is the size of the time window and tmax
% -----------------------------------------------------------------------
% OUTPUT --> GausC, GausS, FGausC, FGausS, R_C and R_S
% -----------------------------------------------------------------------
function [Filter,fftFilter]=...
    RIF_Kernel(sc,ss,x,tsamp,tauC,tauS,tauG,wC,wS,T,tmax)
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
