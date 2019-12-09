%%%%%%%%%%%%%%%%%%%%%%%%%%    RIF_Kernel     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% This function builds the Retina-inspired filter in space and frequency
% -----------------------------------------------------------------------
function [F_Matrix,fftF_Matrix] = retinalFilter(GausC,FGausC,GausS,FGausS,R_C,R_S,wC,wS,samples)
    T = numel(R_C);
    [Gx,Gy] = size(GausC);
    iter = 0;
    for j = 1:samples:T
        iter = iter+1;
        F      =  GausC .* R_C(j) - wS*  GausS .* R_S(j);  % space
        fftF  =  FGausC .* R_C(j) - wS* FGausS .* R_S(j);  % frequency
        F_Matrix(iter,:) = F(:);
        fftF_Matrix(iter,:) = fftF(:);
    end
end
