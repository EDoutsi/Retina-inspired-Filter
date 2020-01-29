%% MAIN FILE WHICH CONNECTS ALL THE FUNCTIONS
% E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079INSTRUCTIONS:

close all; clc;


for p = 1:12
    if p == 1
            class = 'beach';sizep = 2;
    elseif p == 2
            class = 'dogs';sizep = 5;
    elseif p == 3
           class = 'sleeper';sizep = 7;
    elseif p == 4
           class = 'beer_bottle';sizep = 9;
    elseif p == 5
            class = 'retriever';sizep = 9;
    elseif p == 6
            class = 'Stradivarius';sizep = 9;
    elseif p == 7
            class = 'beer_glass';sizep = 9;
    elseif p == 8
            class = 'maillot';sizep = 10;
    elseif p == 9
            class = 'golden_retriever';sizep = 10;
    elseif p == 10
            class = 'social_insect';sizep = 10;
    elseif p == 11
            class = 'steam_iron';sizep = 10;
    elseif p == 12
            class = 'yardstick';sizep = 10;
    end
    for s = 1:sizep;
        %% OPEN FILE
        dir = '/Users/effrosynidoytsi/Documents/Codes/ImageNet_Images/';
        dir = strcat(dir,class,'/image',num2str(s),'.tiff');
        dir1 = '/Users/effrosynidoytsi/Documents/Codes/ImageNet_Images/retriever/image9.tiff';
        %% READ INPUT IMAGE
        I = imread(dir1); 
        I = rgb2gray(I);
        I = double(I);
        %% APPLY THE RIF TRANSFORM
        tic;
        [obs_Matrix,fftobs_Matrix,fftF_Matrix] = RIF_Transform(I);
        %% INVERSE RIF FILTER
        [Iout,MSE,cost_Vector,numbIter] = Inverse_Filtering(I,fftF_Matrix,fftobs_Matrix,MaxIter);
        a = toc
        figure;imagesc(Iout);colormap(gray);
    end
end
%% WEIGHTED DoG CASES
% % [GausC,FGausC] = GaussianKernel(min,max,samp,1);
% % [GausS,FGausS] = GaussianKernel(min,max,samp,1.8);
% % b = 1;
% % a = [0.1 0.3 0.4 0.7 0.9 1.1 2 4];
% % % g = (b*(sc^4))/(a*(ss^4));    
% % figure;
% % for i = 1:size(a,2)
% %     DoG = a(i)*GausC - b*GausS;
% %     fDoG = a(i)*FGausC - b*FGausS; 
% %     subplot(2,8,i); 
% %     plot(x,DoG(floor(Gx/2),:));
% %     subplot(2,8,i+8);
% %     plot(x/2,abs((fDoG(floor(Gx/2),:))))
% % end