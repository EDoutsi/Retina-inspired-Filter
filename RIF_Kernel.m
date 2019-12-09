%%%%%%%%%%%    NON-SEPARABLE SPATIOTEMPORAL FILTER     %%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% This function is responsible to build a bio-inspired filter which is 
% non-separable in space and time. That means we have several decomposition
% layers, all of them have the same size (so there is no pyramidal
% structure), and at the same time each one of them has a different kind of
% a Difference of Gaussian (DoG) filter. 
% -----------------------------------------------------------------------
% FUNCTIONS: 1. Read an image (256x256) or (512x512)
%            2. Construct two 2D-Gaussian filters GausC and GausS
%            3. Construct the temporal filters R_C and R_S 
%            4. Use the R_C and R_S as weights of DoG = GausC - GausS
%            5. The number of samples defines the number of the
%            decomposition layers we are going to consider in the
%            reconstruction. samples could be 2,4,10 or T
% -----------------------------------------------------------------------
% OUTPUT --> Filter
% -----------------------------------------------------------------------


function [F_Matrix,fftF_Matrix] = RIF_Kernel(GausC,FGausC,GausS,FGausS,R_C,R_S,wC,wS,samples)
    T = numel(R_C);
    [Gx,Gy] = size(GausC);
    iter = 0;
    for j = 1:samples:T
        iter = iter+1;
        F      =  GausC .* R_C(j) - wS*  GausS .* R_S(j);    % space
        fftF  =  FGausC .* R_C(j) - wS* FGausS .* R_S(j);  % frequency
        F_Matrix(iter,:) = F(:);
        fftF_Matrix(iter,:) = fftF(:);
        
    end
end