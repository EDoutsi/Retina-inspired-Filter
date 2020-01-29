function [reconstruction,distortion,cost_Vector,numbIter]= RIF_Inverse_Transform(I,fftF_Matrix,fftQobs_Matrix,MaxIter)
            %% GRADIENT DESCEND
            [Ix,Iy] = size(I);
            [FMx,FMy] = size(fftF_Matrix);
            Fx = sqrt(FMy);
            Fy = Fx;
            fftAsum = 0;
            dsum = 0;
            for i = 1:FMx
                fftA = Matrix_reshape(fftF_Matrix, i, Ix, Iy);            % Extract each on of the filters layers
                fftAT = transpose(fftA);                                  % Compute the conjugate of the filter
                fftATA =  fftAT.*fftA;                                    % Compute the positive definite matrix
                fftAsum = fftAsum + fftATA;                               % Sum all the filter layers
                fftobs = Matrix_reshape(fftQobs_Matrix, i, Ix, Iy);       % Extract each on of the observation layers
                d = real(ifft2(fftAT.*fftobs)) ; 
                dsum = dsum  + d;                                         % Sum all the decomposition layers
            end
            xn = zeros(Ix,Iy);
            xn = xn;
            rn = dsum;
            pn = rn;
            for i = 1:MaxIter
                % CGD algortihm
                r1 = sum(rn.*rn); 
                fftpn = fft2(pn);
                temp = real(ifft2(fftAsum.*fftpn));
                r0 = sum(pn.*temp);
                an = r1/r0;               % step
                xn = xn + an*pn;          % solution 
                rn = rn - an*temp;        % residual 
                r2 = sum(rn.*rn);
                bn = r2/r1;               % correction of the search direction 
                pn = rn + bn.*pn;         % search direction 
                % error 
                err = norm(I-xn,'fro')^2/numel(I);
                err_vector(i) = err;
                % cost function
                fftx = fft2(xn);
                h = real(ifft2(fftAsum.*fftx));
                cost = norm(dsum-h,'fro')^2/numel(I);
                cost_Vector(i) = cost;
                cost_fun = 0.5*((sum(xn(:).*h(:)))) - (sum(xn(:).*dsum(:)));
                cost_Function(i) = cost_fun;
                numbIter = i;
                % criterion 
                if i>1
                    if  err_vector(i-1) - err_vector(i)>10e-10
                        if cost < 10e-20
                            disp('Break the gradient descent');
                            break
                        elseif err < 10e-20
                            disp('Break the gradient descent');
                            break
                        end
                    else
                            disp('Break the gradient descent');
                            break
                    end
                end
            end
            distortion = err_vector(end);
            reconstruction = xn;
end
