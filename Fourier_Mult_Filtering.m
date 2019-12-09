%%%%%%%%%%%%%%%%%%     FITLERING IN FOURIER DOMAIN     %%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
% This is a function which does instead of a convolution in space domain, 
% a mulitplication in frequency domain.
% INPUT --> 1. Input image I
%           2. Input filter F
%    !!!! We test the size of the two matrices. If they have the same size
%    then we don't need to add any zero padding. As a result, we use the
%    fft2 (fast fourier transform for 2D) in the default way !!!!
%                          ------BUT-------
%    !!!! In case the two matrices have different size we need to do the
%    zero-padding (padx = Ix+ Fx-1; pady = Iy+Fy-1;) !!!!                      
% ------------------------------------------------------------------------
% OUTPUT --> The filter version in SPACE (obs) and FREQUENCY (fftobs)
%            domain.
% ------------------------------------------------------------------------

function  [obsFreq,obsSpac,fftobsSpectrum,fftobsSpace,fftFSpectrum,fftI] = Fourier_Mult_Filtering(I,fftFSpectrum,F, noise)
        [Ix,Iy] = size(I);
        [Fx,Fy] = size(F);
        if (Ix ~= Fx && Iy ~= Fy) 
                %% COMPUTE THE CONVOLUTION BETWEEN THE FILTER AND THE IMAGE
                % We actually use the mulitplication in fourier domain instead of the 
                % convolution in space. We compute the fft of the input image I 
                % and the filter F. Finally, we calculate the (fftobs) and we inverse
                % this to extract the observations in space. Finally, we save the
                % obs and the fftobs results in matrices.
                h = F;
                ffth = fftFSpectrum;
                s = (size(fftFSpectrum)-1)/2;
                padx = (Ix-Fx)/2;
                pady = (Iy-Fy)/2;
                fftHs = zeros(size(I,1),size(I,2));
                fftHs(padx+1:Fx+padx,pady+1:Fy+pady) = ffth;
                fftHs = fftshift(fftHs);
                fftFSpectrum = fftHs;
                Hs = zeros(size(I,1),size(I,2));
                Hs(padx+1:Fx+padx,pady+1:Fy+pady) = h;
                Hs = fftshift(Hs);
                F = fftshift(Hs);
                fftI = fft2(I);
                fftFSpace = fft2(F);
% %                 figure(40);
% %                 subplot(1,2,1);
% %                 mesh(abs(fftF));
% %                 subplot(1,2,2);
% %                 mesh(abs(FFTofF));
                %% CONVOLUTION in SPACE = MULTIPLICATION in FOURIER                          
                fftobsSpectrum = fftI.*fftFSpectrum;            % using the Gaussian in frequency 
                fftobsSpace      = fftI.* fftFSpace;                 % using the Gaussian in space
            
                obsFreq  = real(ifft2(fftobsSpectrum));         % using the Gaussian in frequency  
                obsSpac = real(fftshift(ifft2(fftobsSpace))); % using the Gaussian in space
                %% NOISE
                obsFreq = obsFreq + noise*randn(Ix,Iy);     
                fftobs = fft2(obsFreq);

% %         F = F0;
% %         padx = Ix+Fx-1;
% %         pady = Iy+Fy-1;
% %         difx = Fx/2;
% %         dify = Fy/2;
% %         fftI = fft2(I,padx,pady);
% %         fftF = fft2(F,padx,pady);
% %         fftobs = fftI.*fftF;
% %         obs = real(ifft2(fftobs));
    else
            fftI = fft2(I);
% %             fftF = fft2(F0);
            fftobs = fftI.*fftFSpectrum;
            obsFreq = real(ifft2(fftobs));
    end
% %     figure(2);
% %     subplot(2,3,1);
% %     imagesc(I);colormap(gray);
% %     title('Original Signal');
% %     subplot(2,3,2);
% %     mesh(F0);
% %     title('Filter in space');
% %     subplot(2,3,3);
% %     imagesc(obs);colormap(gray);
% %     title('Observations');
% %     subplot(2,3,4);
% %     imagesc(log(abs(fftshift(fftI))));colormap(gray);
% %     title('FFT of the signal');
% %     subplot(2,3,5);
% %     imagesc(log(abs(fftshift(fftF))));
% %     title('FFT of the filter');
% %     subplot(2,3,6);
% %     imagesc(log(abs(fftshift(fftobs))));colormap(gray);
% %     title('FFT of the observation');
% %     
    
end