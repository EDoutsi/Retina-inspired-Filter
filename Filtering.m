%%%%%%%%%%%%%%%%%%     FITLERING IN FOURIER DOMAIN     %%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
% This is a function which does instead of a convolution in space domain, 
% a mulitplication in frequency domain.
% INPUT --> 1. Image "I"
%           2. RIF Filter "fftFilter" in frequency domain     
%           3. RIF Filter "Filter" in space domain
%           4. "noise" (default 0)
% ------------------------------------------------------------------------
% ALGORITHM --> We test the size of the two matrices "I" and "Filter". 
% If they have the same size then we don't need to add any zero padding. 
% As a result, we use the fft2 (fast fourier transform for 2D) in the default 
% way. In case the two matrices have different size we need to do the zero-padding.                 
% ------------------------------------------------------------------------
% OUTPUT --> 1. The filtered images "obs_Matrix" in space domain
%            2. The filtered images "fftobs_Matrix" in frequency domain
%            3. The zero-padded filter "fftF_Matrix" in frequency domain
% ------------------------------------------------------------------------
function  [obs_Matrix,fftobs_Matrix,fftF_Matrix] = Filtering(I,fftFilter, Filter, noise)
        [Ix,Iy] = size(I);
        [Fx,Fy] = size(Filter);
        for i = 1:Fx                    
                F = Matrix_reshape(Filter,i,sqrt(Fy),sqrt(Fy));                  
                fftFSpectrum = Matrix_reshape(fftFilter,i,sqrt(Fy),sqrt(Fy));    
                if (Iy ~= Fy)
                    %% COMPUTE THE CONVOLUTION BETWEEN THE FILTER AND THE IMAGE
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
                    fftWhite = fft2(noie*randn(Ix,Iy));
                    fftobsSpectrum = fftI.*fftFSpectrum;            % using the Gaussian in frequency 
                    fftobsSpace      = fftI.* fftF + fftWhite;                 % using the Gaussian in space
                    obs = real((ifft2(fftobsSpace))); % using the Gaussian in space
                else
                    fftI = fft2(I);
                    fftobsSpace = fftI.*fftF;
                end
                obs_Matrix(i,:) = obs(:);
                fftobs_Matrix(i,:) = fftobsSpace(:);
                fftF_Matrix(i,:) = fftF(:);
        end
end
