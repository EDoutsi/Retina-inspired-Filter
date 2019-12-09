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
[GausC,FGausC,GausS,FGausS,R_C,R_S] = Build_Filters(sc,ss,max,min,samp,tauC,tauS,tauG,wC,wS,T,tmax,t1);
[Gx,Gy] = size(GausC);                 
[fGx , fGy] = freqspace([Gx Gy]);  

%% BUILD THE RETINAL-INSPIRED FILTER IN SPACE AND FREQUENCY

[Filter,fftFilter] = retinalFilter(GausC, FGausC, GausS, FGausS, R_C, R_S, wC, wS, tsamp); 
[FMx,FMy] = size(Filter);

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

%% PLOT THE RETINA-INSPIRED FILTER IN SPACE
for i = 1:FMx
        F = Filter(i,:);
        filter = reshape(F,Fx,Fy);
        fftfilter = fftFilter(i,:);
        F = reshape(fftfilter,Fx,Fy);
        filter    = Matrix_reshape(Filter,i,Gx,Gy);
        fftfilter = Matrix_reshape(fftFilter,i,Gx,Gy);
        maxF    = filter(floor(Gx/2)+1,:);
        figure(21);
        hold on;
        clr = lines(FMy); 
        lineWidth = (i-1)*0.1+0.1;
        cut_maxF(i,:) = maxF;
        maxSpec = fftfilter(floor(Gx/2)+1,:);
        cut_spec(i,:) = maxSpec;
        subplot(1,2,1);
        plot(2*pi*fGx,abs(maxSpec),'LineWidth',lineWidth,'Color',clr(i,:));
        title('transversal Cut of the Spectrum')
        hold on;
        subplot(1,2,2);
        plot(x,maxF,'LineWidth',lineWidth,'Color',clr(i,:));
        title('transversal Cut of the Filter in space')
end
hold off;
str = cellstr( num2str((1:tsamp)','time bin %d') );
h = legend('$t_{1}$=1msec','$t_{2}$=30msec','$t_{3}$=60msec','$t_{4}$=90msec','$t_{5}$=120msec','Location','Best' );
set(h, 'Interpreter', 'latex')

figure(23);
mesh(x, t1, cut_maxF);
xlabel('r');
ylabel('time');
colorbar('location','southoutside');
figure(24);
subplot(1,4,1);
mesh(x, t1, cut_maxF);
xlabel('x');
ylabel('time');
subplot(1,4,2);
imagesc(x, t1, cut_maxF);
xlabel('x');
ylabel('time');
subplot(1,4,3);
mesh(2*pi*fGx,t1,abs(cut_spec));axis tight;
xlabel('omega');
ylabel('time');
subplot(1,4,4);
imagesc(2*pi*fGx,t1,abs(cut_spec));
xlabel('omega');
ylabel('time');

%% PLOT THE BANDWIDTH OF THE RETINA-INSPIRED FILTER  
[fmin,fmax1, fmax2,Spectrum_min1,Spectrum_max1,Spectrum_max2,band1,band2,band3,f_bdwth,omegaTotal] = study_RC_and_RS(wC,wS,R_C,R_S,ss,sc);
wmin    = 2*pi*fmin;
wmax1   = 2*pi*fmax1;
wmax2   = 2*pi*fmax2;
w_bdwth1 = 2*pi*f_bdwth(1,:);
w_bdwth2 = 2*pi*f_bdwth(2,:);

figure(42);
plot(wmin,t1,'b',w_bdwth1,t1,'m-',-w_bdwth1,t1,'m-',w_bdwth2,t1,'g-',-w_bdwth2,t1,'g-',omegaTotal(3,:),t1,'m-.',-omegaTotal(3,:),t1,'m-.',omegaTotal(2,:),t1,'g-.',-omegaTotal(2,:),t1,'g-.');
title('Approximation');
h = legend('$\omega=0$', '$\omega_{H}^{exact}$','$-\omega_{H}^{exact}$','$\omega_{L}^{exact}$','$-\omega_{L}^{exact}$','$\omega_{H}^{approx}$','$-\omega_{H}^{approx}$','$-\omega_{L}^{approx}$','$-\omega_{L}^{approx}$','Location','Best');
xlabel('$\omega$');
ylabel('time');
set(h, 'Interpreter', 'latex')
figure(43);
plot(wmin,t1,'b',omegaTotal(3,:),t1,'m-.',-omegaTotal(3,:),t1,'m-.',omegaTotal(2,:),t1,'g--',-omegaTotal(2,:),t1,'g--');
h = legend('$\omega=0$', '$\omega_{H}$','$-\omega_{H}$','$\omega_{L}$','$-\omega_{L}$','$-\omega_{L}$','$-\omega_{L}$','Location','Best');
title('Exact Solution');
