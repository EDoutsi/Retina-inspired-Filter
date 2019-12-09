%%%%%%%%%%%%%%%%%%   COMPONENTS OF THE RIF FILTER     %%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% FUNCTIONS: 1. Build the 2D-Gaussian filters GausC and GausS
%            2. Build the temporal filters R_C and R_S 
% -----------------------------------------------------------------------
% OUTPUT --> GausC, GausS, FGausC, FGausS, R_C and R_S
% -----------------------------------------------------------------------
function [GausC,FGausC,GausS,FGausS,R_C,R_S]=...
    Build_Filters(sc,ss,max,min,samp,tauC,tauS,tauG,wC,wS,T,tmax,t1)
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

    %% -- PLOT THE FILTERS 
    figure;
    subplot(2,3,1);
    mesh(x,x, GausC);
    title('Center Gaussian - space');
    subplot(2,3,2);
    mesh(fGx, fGy, abs(FGausC));
    title('Center Gaussian - frequency');
    subplot(2,3,3);
    plot(t1,R_C,'b',t1,R_S,'r');
    xlabel('time');
    subplot(2,3,4);
    mesh(x,x, GausS);
    title('Surround Gaussian - space');
    subplot(2,3,5);
    mesh(fGx, fGy,abs(FGausS));
    title('Surround Gaussian - frequency')
    subplot(2,3,6);
    plot(t2,R_C_2T,'b',t2,R_S_2T,'r');
    xlabel('time');axis tight;
end