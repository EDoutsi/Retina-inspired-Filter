%%%%%%%%%%%%%%%%%%   COMPONENTS OF THE RIF FILTER     %%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% FUNCTIONS: 1. Build the 2D-Gaussian filters GausC and GausS
%            2. Build the temporal filters R_C and R_S 
% -----------------------------------------------------------------------
% OUTPUT --> GausC, GausS, FGausC, FGausS, R_C and R_S
% -----------------------------------------------------------------------
function [Filter,fftFilter]=...
    RIF_Kernel(sc,ss,max,min,samp,tsamp,tauC,tauS,tauG,wC,wS,T,tmax,t1)
    %% -- 2D GUASSIAN FILTERS
    [GausC, FGausC] = GaussianKernel(min,max,samp,sc); % in tests 11
    [GausS, FGausS] = GaussianKernel(min,max,samp,ss); % in tests 11
    [Gx,Gy] = size(GausC);
    [fGx , fGy] = freqspace([Gx Gy]);
    x = linspace(min,max,samp);
    %% -- 1D TEMPORAL FILTERS
    R_C = ComputingRc(T, tmax, wC, tauC, tauG);
    R_S = ComputingRs(T, tmax, wS, tauC, tauG, tauS);
    tmax2 = 3*T;
    t2 = 1:tmax2;
    R_C_2T = ComputingRc(T, tmax2, wC, tauC, tauG);
    R_S_2T = ComputingRs(T, tmax2, wS, tauC, tauG, tauS);
    %% -- 2D RIF Filter 
    for j = 1:tsamp:T
        Filter(j,:)  = wC * GausC(:) .* R_C(j) - wS * GausS(:) .* R_S(j);    % space
        fftFilter(j,:)  = wC * FGausC(:) .* R_C(j) - wS * FGausS(:) .* R_S(j);  % frequency
    end
end
