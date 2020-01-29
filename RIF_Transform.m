% E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079INSTRUCTIONS:
%
% Generate the Retinal-Inspired Filter (RIF). 
     

close all; clc;

%% SET OF BIO-PLAUSIBLE PARAMETERS 
sc  = 0.5; ss  = 3*sc;      % bio-plausible parameters 
maxx = 10; minx  = -maxx; samp = 50*maxx+1; % 100*max+1 in the experiments
tauC = 20.*10^-3; tauS = 4.*10^-3; tauG = 5.*10^-3;
wC = 1; wS = 1;             % weight surround
T  = 150;  tmax = T;  tsamp = 1;    % number odd time instances
t1 = 1:tsamp:tmax;                  % temporal vector
x  = linspace(minx,maxx,samp);        % spatial vector

%% BUILD THE SPATIAL AND TEMPORAL FILTERS
[Filter,fftFilter] = RIF_Kernel(sc,ss,x,tsamp,tauC,tauS,tauG,wC,wS,T,tmax,t1);
