%%%%%%%%%%%%%%%%%%     FITLERING IN FOURIER DOMAIN     %%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
% This is a function which does instead of a convolution in space domain, 
% a mulitplication in frequency domain.
% INPUT --> 1. Image "I"
%           2. RIF Filter "fftFilter" in frequency domain     
%           3. RIF Filter "Filter" in space domain
%           4. "noise" (default 0)
% ------------------------------------------------------------------------
% ALGORITHM --> We test the size of the two matrices "I" and "Filter". If they have the same
% size then we don't need to add any zero padding. As a result, we use the 
% fft2 (fast fourier transform for 2D) in the default way. In case the two 
% matrices have different size we need to do the zero-padding                 
% ------------------------------------------------------------------------
% OUTPUT --> 1. The filtered images "obs_Matrix" in space domain
%            2. The filtered images "fftobs_Matrix" in frequency domain
%            3. The zero-padded filter "fftF_Matrix" in frequency domain
% ------------------------------------------------------------------------

function  [obs_Matrix,fftobs_Matrix,fftF_Matrix] = Filtering(I,fftFilter, Filter, noise)
        [Ix,Iy] = size(I);
        [Fx,Fy] = size(Filter);
        for i = 1:Fx                    
                F = Matrix_reshape(Filter,i,sqrt(Fy),sqrt(Fy));                          % extract 1-1 the decomposition layers without noise
                fftFSpectrum = Matrix_reshape(fftFilter,i,sqrt(Fy),sqrt(Fy));            % extract 1-1 the decomposition layers without noise
                if (Iy ~= Fy)
                    %% COMPUTE THE CONVOLUTION BETWEEN THE FILTER AND THE IMAGE
                    % We actually use the mulitplication in fourier domain instead of the 
                    % convolution in space. We compute the fft of the input image I 
                    % and the filter F. Finally, we calculate the (fftobs) and we inverse
                    % this to extract the observations in space. Finally, we save the
                    % obs and the fftobs results in matrices.
                    h = F;
                    ffth = fftFSpectrum;
                    s = (size(fftFSpectrum)-1)/2;
                    padx = round((Ix-sqrt(Fy))/2); 
                    pady = round((Iy-sqrt(Fy))/2);
                    fftHs = zeros(size(I,1),size(I,2));
                    fftHs(padx+1:sqrt(Fy)+padx,pady+1:sqrt(Fy)+pady) = ffth;
                    fftHs = fftshift(fftHs);
                    fftFSpectrum = fftHs;
                    Hs = zeros(size(I,1),size(I,2));
                    Hs(padx+1:sqrt(Fy)+padx,pady+1:sqrt(Fy)+pady) = h;
                    F = fftshift(Hs);
                    fftI = fft2(I);
                    fftF = fft2(F);

                   %% CONVOLUTION in SPACE = MULTIPLICATION in FOURIER                          
                    fftWhite = fft2(noise*randn(Ix,Iy));
                    fftobsSpectrum = fftI.*fftFSpectrum;            % using the Gaussian in frequency 
                    fftobsSpace      = fftI.* fftF + fftWhite;                 % using the Gaussian in space
                    obs = real((ifft2(fftobsSpace))); % using the Gaussian in space
% %                     figure(1);
% %                     imagesc(obs);colormap(gray);
                    %% NOISE
% %                     obsFreq = obsFreq + noise*randn(Ix,Iy);     
% %                     fftobs = fft2(obsFreq);
                else
                    fftI = fft2(I);
                    fftobsSpace = fftI.*fftF;
                end
                obs_Matrix(i,:) = obs(:);
                fftobs_Matrix(i,:) = fftobsSpace(:);
                fftF_Matrix(i,:) = fftF(:);
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