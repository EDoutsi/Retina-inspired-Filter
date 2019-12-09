
%%%%%%%%%%%%%%%%% BUILD THE CENTER TEMPORAL FILTER  %%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% Description: Here we compute the temporal center filter R_C which is the
% weight for for the surround spatial Gaussian filter. This definition is 
% according to the mathematical proof and the results we obtained in [1].
%
%[1] E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079
%
% INPUTS --> 1. T    : is the length of the temporal vector
%            2. wC   : is a constant weight parameter
%            3. tauC : is a constant that defines slope exponential function
%            4. tauG : is a constant of the gamma filter
%            
% OUTPUTS --> The temporal 1D surround filter R_C 

function [ R_C ] = ComputingRc(T, tmax, wC , tauC, tauG )
        n = 5;
        A = 1./tauG ; 
        B = (tauC - tauG)./(tauC.*tauG);
        Conv1 = zeros(1,tmax);  
        for i = 0:tmax-1
             sum1 = zeros(1, n+1);
             sum11 = zeros(1, n+1);
             sum12 = zeros(1, n+1);
             if (i <= T)
                for k = 0:n
                    sum1(k+1) = - (factorial(n)./ (factorial(n-k).* (tauG.^(n+1)) .* (A.^(k+1)))).* ((i.*10^-3).^(n-k) .* exp(-(i.*10^-3)/tauG));
                end
                Conv1(i+1) = sum(sum1) + factorial(n);
             else
                for k = 0:n
                    sum11(k+1) = - (factorial(n)./ (factorial(n-k).* (tauG.^(n+1)) .* (A.^(k+1)))).* ((i.*10^-3).^(n-k) .* exp(-(i.*10^-3)/tauG));
                    sum12(k+1) = - (factorial(n)./ (factorial(n-k).* (tauG.^(n+1)) .* (A.^(k+1)))).* (((i-T).*10^-3).^(n-k) .* exp(-((i-T).*10^-3)/tauG));
                end
                Conv11 = sum(sum11) + factorial(n);
                Conv12 = sum(sum12) + factorial(n);
                Conv1(i+1) = Conv11 - Conv12;
             end
        end      

        Conv2 = zeros(1,tmax);
        for i = 0:tmax-1
              sum2 = zeros(1, n+1);
              sum21 = zeros(1, n+1);
              sum22 = zeros(1, n+1);
              if (i <= T)
                    for k = 0:n
                        m = n-k;
                        sum1 = zeros(1,m+1);
                        for l = 0:m;
                            sum1(l+1) = (factorial(n)./( factorial(m-l).* (B.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauC)).* ((i.*10^-3).^(m-l)) .* exp(-A.*(i.*10^-3));
                        end
                        sum2(k+1) = sum(sum1) - ((factorial(n))./((B.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauC));
                    end
                    Conv2(i+1)  = sum(sum2) + (factorial(n).*(1-exp(-(i.*10^-3)./tauC)))./(B.^(n+1).*(tauG.^(n+1)));
              else
                  for k = 0:n
                        m = n-k;
                        sum11 = zeros(1,m+1);
                        sum12 = zeros(1,m+1);
                        for l = 0:m;
                            sum11(l+1) = (factorial(n)./( factorial(m-l).* (B.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauC)).* ((i.*10^-3).^(m-l)) .* exp(-A.*(i.*10^-3));
                            sum12(l+1) = (factorial(n)./( factorial(m-l).* (B.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauC)).* (((i-T).*10^-3).^(m-l)) .* exp(-A.*((i-T).*10^-3));
                        end
                        sum21(k+1) = sum(sum11) - ((factorial(n))./((B.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauC));
                        sum22(k+1) = sum(sum12) - ((factorial(n))./((B.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauC));
                  end
                  Conv21  = sum(sum21) + (factorial(n).*(1-exp(-((i).*10^-3)./tauC)))./(B.^(n+1).*(tauG.^(n+1)));
                  Conv22  = sum(sum22) + (factorial(n).*(1-exp(-((i-T).*10^-3)./tauC)))./(B.^(n+1).*(tauG.^(n+1)));
                  Conv2(i+1) = Conv21 - Conv22 ;
              end    
        end
        R_C = Conv1 - wC.*Conv2 ;
        % % figure;
        % % plot(R_C);
        % % title('center filter');   
end
