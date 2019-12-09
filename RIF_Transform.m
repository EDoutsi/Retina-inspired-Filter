%% MAIN FILE WHICH CONNECTS ALL THE FUNCTIONS
% E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079INSTRUCTIONS:
%
% 1. Build the retinal-inspired frame.
%       1.1 We need first to build the Gaussian Filters (GaussianKernel.m)
%       1.2 The temporal filters (ComputingRc.m and ComputingRs.m)
%       1.3 Build the retinal-inspired non-separable spatiotemporal filter (retinalFilter.m)

close all; clc;

%% SET OF BIO-PLAUSIBLE PARAMETERS 
sc  = 0.5;       % standard deviation sc = 0.5 according to Wohrer et al. 
ss  = 3*sc;      % bio-plausible ss = 3*sc according to Wohrer et al.
max = 10;        % 4 in the experiments
min  = -max;  
samp = 50*max+1; % 100*max+1 in the experiments
tauC = 20.*10^-3;
tauS = 4.*10^-3;
tauG = 5.*10^-3;
wC = 1;          % weight  center/ constant valq  ue wC=0.75
wS = 1;          % weight surround
T  = 150;        % bio-plausible processing time
tmax = T;        % i) t = T or ii) t > T  (the difference is in the convergence)
tsamp = 1;       % number od time instances
t1 = 1:tsamp:tmax;              % temporal vector
x  = linspace(min,max,samp);    % spatial vector

%% BUILD THE SPATIAL AND TEMPORAL FILTERS
[Filter,fftFilter] = RIF_Kernel(sc,ss,max,min,samp,tsamp,tauC,tauS,tauG,wC,wS,T,tmax,t1);
[Gx,Gy] = size(GausC);                 
[fGx , fGy] = freqspace([Gx Gy]);  

%% WEIGHTED DoG CASES
[GausC,FGausC] = GaussianKernel(min,max,samp,1);
[GausS,FGausS] = GaussianKernel(min,max,samp,1.8);
b = 1;
a = [0.1 0.3 0.4 0.7 0.9 1.1 2 4];
% g = (b*(sc^4))/(a*(ss^4));    
figure;
for i = 1:size(a,2)
    DoG = a(i)*GausC - b*GausS;
    fDoG = a(i)*FGausC - b*FGausS; 
    subplot(2,8,i); 
    plot(x,DoG(floor(Gx/2),:));
    subplot(2,8,i+8);
    plot(x/2,abs((fDoG(floor(Gx/2),:))))
end