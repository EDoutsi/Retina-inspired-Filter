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
        Filter(j,:)  = wC * GausC(:) .* R_C(j) - wS * GausS(:) .* R_S(j);    % space
        fftFilter(j,:)  = wC * FGausC(:) .* R_C(j) - wS * FGausS(:) .* R_S(j);  % frequency
    end
    [FMx,FMy] = size(Filter);
     %% -- PLOT THE GAUSSIAN AND EXPONENTIAL FILTERS 
% %     figure;
% %     subplot(2,3,1);
% %     mesh(x,x, GausC);
% %     title('Center Gaussian - space');
% %     subplot(2,3,2);
% %     mesh(fGx, fGy, abs(FGausC));
% %     title('Center Gaussian - frequency');
% %     subplot(2,3,3);
% %     plot(t1,R_C,'b',t1,R_S,'r');
% %     xlabel('time');
% %     subplot(2,3,4);
% %     mesh(x,x, GausS);
% %     title('Surround Gaussian - space');
% %     subplot(2,3,5);
% %     mesh(fGx, fGy,abs(FGausS));
% %     title('Surround Gaussian - frequency')
% %     subplot(2,3,6);
% %     plot(t2,R_C_2T,'b',t2,R_S_2T,'r');
% %     xlabel('time');axis tight;
    
    %% PLOT THE RETINA-INSPIRED FILTER IN SPACE
% %     for i = 1:FMx
% %             filter     = Matrix_reshape(Filter,i,Gx,Gy);
% %             fftfilter = Matrix_reshape(fftFilter,i,Gx,Gy);
% %             maxF    = filter(floor(Gx/2)+1,:);
% %             figure(21);
% %             hold on;
% %             clr = lines(FMy); 
% %             lineWidth = (i-1)*0.1+0.1;
% %             cut_maxF(i,:) = maxF;
% %             maxSpec = fftfilter(floor(Gx/2)+1,:);
% %             cut_spec(i,:) = maxSpec;
% %             subplot(1,2,1);
% %             plot(2*pi*fGx,abs(maxSpec),'LineWidth',lineWidth,'Color',clr(i,:));axis tight;
% %             title('transversal Cut of the Spectrum')
% %             hold on;
% %             subplot(1,2,2);
% %             plot(x,maxF,'LineWidth',lineWidth,'Color',clr(i,:));axis tight;
% %             title('transversal Cut of the Filter in space')
% %     end
% %     hold off;
% %     str = cellstr( num2str((1:tsamp)','time bin %d') );
% %     h = legend('$t_{1}$=1msec','$t_{2}$=30msec','$t_{3}$=60msec','$t_{4}$=90msec','$t_{5}$=120msec','Location','Best' );
% %     set(h, 'Interpreter', 'latex')
% % 
% %     figure(23);
% %     mesh(x, t1, cut_maxF);
% %     xlabel('r');
% %     ylabel('time');
% %     colorbar('location','southoutside');
% %     figure(24);
% %     subplot(1,4,1);
% %     mesh(x, t1, cut_maxF);
% %     xlabel('x');
% %     ylabel('time');
% %     subplot(1,4,2);
% %     imagesc(x, t1, cut_maxF);
% %     xlabel('x');
% %     ylabel('time');
% %     subplot(1,4,3);
% %     mesh(2*pi*fGx,t1,abs(cut_spec));axis tight;
% %     xlabel('omega');
% %     ylabel('time');
% %     subplot(1,4,4);
% %     imagesc(2*pi*fGx,t1,abs(cut_spec));
% %     xlabel('omega');
% %     ylabel('time');

    %% PLOT THE BANDWIDTH OF THE RETINA-INSPIRED FILTER  
% %     [fmin,fmax1, fmax2,Spectrum_min1,Spectrum_max1,Spectrum_max2,band1,band2,band3,f_bdwth,omegaTotal] = study_RC_and_RS(wC,wS,R_C,R_S,ss,sc);
% %     wmin    = 2*pi*fmin;
% %     wmax1   = 2*pi*fmax1;
% %     wmax2   = 2*pi*fmax2;
% %     w_bdwth1 = 2*pi*f_bdwth(1,:);
% %     w_bdwth2 = 2*pi*f_bdwth(2,:);
% % 
% %     figure(42);
% %     plot(wmin,t1,'b',w_bdwth1,t1,'m-',-w_bdwth1,t1,'m-',w_bdwth2,t1,'g-',-w_bdwth2,t1,'g-',omegaTotal(3,:),t1,'m-.',-omegaTotal(3,:),t1,'m-.',omegaTotal(2,:),t1,'g-.',-omegaTotal(2,:),t1,'g-.');
% %     title('Approximation');
% %     h = legend('$\omega=0$', '$\omega_{H}^{exact}$','$-\omega_{H}^{exact}$','$\omega_{L}^{exact}$','$-\omega_{L}^{exact}$','$\omega_{H}^{approx}$','$-\omega_{H}^{approx}$','$-\omega_{L}^{approx}$','$-\omega_{L}^{approx}$','Location','Best');
% %     xlabel('$\omega$');
% %     ylabel('time');
% %     set(h, 'Interpreter', 'latex')
% %     figure(43);
% %     plot(wmin,t1,'b',omegaTotal(3,:),t1,'m-.',-omegaTotal(3,:),t1,'m-.',omegaTotal(2,:),t1,'g--',-omegaTotal(2,:),t1,'g--');
% %     h = legend('$\omega=0$', '$\omega_{H}$','$-\omega_{H}$','$\omega_{L}$','$-\omega_{L}$','$-\omega_{L}$','$-\omega_{L}$','Location','Best');
% %     title('Exact Solution');

end