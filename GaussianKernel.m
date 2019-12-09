%%%%%%%%%%%%%%%%%%%%%    BUILD THE GAUSSIAN KERNEL    %%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% INPUTS-->  1. min: the min value of the linspace
%            2. max: The max value of the linspace
%            3. samp: the number of points between min and max
%            4. s: is the standard deviation of the Gaussian which defines the spread
%
% OUTPUT-->  1. Gaus2D: the 2D Gaussian filter
%            2. FGaus2D: the spectrum of the filter

%% Gaussian Filter 
function [GausC2D,FGaus2D] = GaussianKernel(min,max,samp,s)
    x = linspace(min,max,samp);
    [X,Y] = meshgrid(x,x);
    %% SPACE
    GausC2D = (1/(2*pi*(s^2)))*exp( (-(Y.^2)-(X.^2))/(2*(s^2)) );
    % GausC2D = GausC2D/sum(GausC2D(:)); % normalize the filter
    [Gx,Gy] = size(GausC2D);
    [fGx,fGy] = freqspace([Gx Gy]);
    %% FREQUENCY
    [fX,fY] = meshgrid(fGx,fGy);
    FGaus2D =1/(2*pi)* (exp((-(2*pi*fX*s).^2 - (2*pi*fY*s).^2)/2 ));
    % % figure(9);
    % % subplot(1,2,1);
    % % mesh(x,x,(GausC2D));title('Gaussian'); set(gca,'Fontsize',12);
    % % subplot(1,2,2);
    % % mesh(fGx,fGx,spectrum2D);title('FFT Gaussian'); set(gca,'Fontsize',12);
    % % figure(10);
    % % subplot(1,2,1);
    % % mesh(fGx,fGx,abs(fftshift(fft2(GausC2D))));
    % % subplot(1,2,2);
    % % mesh(fGx,fGy,spectrum2D);
end