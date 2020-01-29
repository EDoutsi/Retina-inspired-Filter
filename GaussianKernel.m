%%%%%%%%%%%%%%%%%%%%%    BUILD THE GAUSSIAN KERNEL    %%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% INPUTS-->  1. x: linespace
%            2. s: is the standard deviation of the Gaussian which defines the spread
% -----------------------------------------------------------------------
% OUTPUT-->  1. Gaus2D: the 2D Gaussian filter
%            2. FGaus2D: the spectrum of the filter
% -----------------------------------------------------------------------
function [GausC2D,fftGaus2D] = GaussianKernel(x,s)
    [X,Y] = meshgrid(x,x);
    %% SPACE
    GausC2D = (1/(2*pi*(s^2)))*exp( (-(Y.^2)-(X.^2))/(2*(s^2)) );
    [Gx,Gy] = size(GausC2D);
    [fGx,fGy] = freqspace([Gx Gy]);
    %% FREQUENCY
    [fX,fY] = meshgrid(fGx,fGy);
    FGaus2D = 1/(2*pi)* (exp((-(2*pi*fX*s).^2 - (2*pi*fY*s).^2)/2 ));
    Gaus2D_nor = GausC2D./sum(GausC2D(:));
    fftGaus2D = fftshift(fft2(Gaus2D_nor));
end
