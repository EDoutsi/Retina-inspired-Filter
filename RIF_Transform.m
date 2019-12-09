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
sc  = 0.5; ss  = 3*sc;      % bio-plausible parameters 
max = 10; min  = -max; samp = 50*max+1; % 100*max+1 in the experiments
tauC = 20.*10^-3; tauS = 4.*10^-3; tauG = 5.*10^-3;
wC = 1; wS = 1;             % weight surround
T  = 150;  tmax = T;  tsamp = 1;    % number odd time instances
t1 = 1:tsamp:tmax;                  % temporal vector
x  = linspace(min,max,samp);        % spatial vector

%% BUILD THE SPATIAL AND TEMPORAL FILTERS
[Filter,fftFilter] = RIF_Kernel(sc,ss,max,min,samp,tsamp,tauC,tauS,tauG,wC,wS,T,tmax,t1);
