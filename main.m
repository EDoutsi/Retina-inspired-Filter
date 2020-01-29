%% MAIN FILE WHICH CONNECTS ALL THE FUNCTIONS
% E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079INSTRUCTIONS:

close all; clc;

% LOAD THE INPUT IMAGE
I = double(rgb2gray(imread('testimage1.tiff'))); 

%% APPLY THE RIF TRANSFORM
[RIF_I,fftRIF_I,fftRIF_Filter] = RIF_Transform(I);

%% APPLY THE INVERSE RIF TRANSFORM
[Iout,MSE] = Inverse_Filtering(I,fftRIF_Filter,fftRIF_I,10000);
