
%%%%%%%%%%%%%%%%% BUILD THE SURROUND TEMPORAL FILTER  %%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------------
% Description: Here we compute the temporal surround filter R_S which is the
% weight for for the surround spatial Gaussian filter. This definition is 
% according to the mathematical proof and the results we obtained in [1].
%
%[1] E. Doutsi, L. Fillatre, M. Antonini and J. Gaulmin, "Retina-Inspired Filter," 
% in IEEE Transactions on Image Processing, vol. 27, no. 7, pp. 3484-3499, July 2018. 
% doi: 10.1109/TIP.2018.2812079
%
% INPUTS --> 1. T    : is the length of the temporal vector
%            2. wC   : is a constant parameter
%            3. tauC : is a constant that defines slope center exponential function
%            4. tauG : is the parameter which defined the gamma filter
%            5. tauS : is a constant that defines slope surround exponential function
%
% OUTPUTS --> The temporal 1D surround filter R_S 

function [ R_S ] = ComputingRs( T, tmax, wC , tauC, tauG, tauS )
        n = 5;
        A = 1./tauG ; 
        B = (tauC - tauG)./(tauC.*tauG);
        g = (tauS - tauG)./(tauG.*tauS); 
        phi = (tauS - tauC)./(tauC.*tauS);
        At = zeros(1,tmax);
        for i = 0:tmax-1
               sum2 = zeros(1, n+1);
               sum21 = zeros(1, n+1);
               sum22 = zeros(1, n+1);
               if(i <= T)
                    for k = 0:n
                        m = n-k;
                        sum1 = zeros(1,m+1);
                        for l = 0:m;
                            sum1(l+1) = (factorial(n)./( factorial(m-l).* (g.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauS)).* ((i.*10^-3).^(m-l)) .* exp(-A.*(i.*10^-3));
                        end
                        sum2(k+1) = sum(sum1) - (factorial(n)./((g.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauS));
                    end
                    At(i+1)  = sum(sum2) + (factorial(n)./(g.^(n+1).*(tauG.^(n+1)))) .*(1 - exp(-(i.*10^-3)./tauS));
               else
                   for k = 0:n
                        m = n-k;
                        sum11 = zeros(1,m+1);
                        sum12 = zeros(1,m+1);
                        for l = 0:m;
                            sum11(l+1) = (factorial(n)./( factorial(m-l).* (g.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauS)).* ((i.*10^-3).^(m-l)) .* exp(-A.*(i.*10^-3));
                            sum12(l+1) = (factorial(n)./( factorial(m-l).* (g.^(k+1)) .* (A.^(l+1)).*(tauG.^(n+1)).*tauS)).* (((i-T).*10^-3).^(m-l)) .* exp(-A.*((i-T).*10^-3));
                        end
                        sum21(k+1) = sum(sum11) - (factorial(n)./((g.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauS));
                        sum22(k+1) = sum(sum12) - (factorial(n)./((g.^(k+1)).*(A.^(m+1)).*(tauG.^(n+1)).*tauS));
                    end
                    At1  = sum(sum21) + (factorial(n)./(g.^(n+1).*(tauG.^(n+1)))) .*(1 - exp(-(i.*10^-3)./tauS));
                    At2  = sum(sum22) + (factorial(n)./(g.^(n+1).*(tauG.^(n+1)))) .*(1 - exp(-((i-T).*10^-3)./tauS));
                    At(i+1) = At1 - At2;
               end 
        end

        Bt = zeros (1, tmax);
        for i = 0:tmax-1
              sum3 = zeros (1 , n+1);
              sum31 = zeros (1 , n+1);
              sum32 = zeros (1 , n+1);
              if (i <= T)
                  for k= 0:n
                        m = n-k;
                        sum2 = zeros(1, m+1);
                        for l = 0:m
                            p = m-l;
                            sum1 = zeros(1, p+1);
                            for r = 0:p
                                sum1(r+1) = (- factorial(n) ./ (factorial(p-r) .*(B.^(k+1)).* (g.^(l+1)) .* (A.^(r+1)) .* (tauG.^(n+1)) .* tauC .* tauS)) .* ((i.*10^-3).^(p-r)).* exp(-A.*(i.*10^-3));
                            end
                            sum2(l+1) = sum(sum1) + (factorial(n) ./ ((B.^(k+1)).* (g.^(l+1)) .* (A.^(p+1)).* (tauG.^(n+1)) .* tauC .*tauS));
                        end
                        sum3(k+1) = sum(sum2) - (factorial(n) ./ ((B.^(k+1)).* (g.^(m+1)).* (tauG.^(n+1)) .* tauC)) .*(1 - exp(-(i.*10^-3)./tauS));
                  end
                  Bt(i+1) = sum(sum3) + (factorial(n)./((B.^(n+1)).* phi.* (tauG.^(n+1)) .* tauC)).*(1 - exp(-(i.*10^-3)./tauS)) - (factorial(n)./((B.^(n+1)).* phi .* (tauG.^(n+1)) .*tauS)).*(1 - exp(-(i.*10^-3)./tauC));
              else
                  for k= 0:n
                        m = n-k;
                        sum21 = zeros(1, m+1);
                        sum22 = zeros(1, m+1);
                        for l = 0:m
                            p = m-l;
                            sum11 = zeros(1, p+1);
                            sum12 = zeros(1, p+1);
                            for r = 0:p
                                sum11(r+1) = (- factorial(n) ./ (factorial(p-r) .*(B.^(k+1)).* (g.^(l+1)) .* (A.^(r+1)) .* (tauG.^(n+1)) .* tauC .* tauS)) .* ((i.*10^-3).^(p-r)).* exp(-A.*(i.*10^-3));
                                sum12(r+1) = (- factorial(n) ./ (factorial(p-r) .*(B.^(k+1)).* (g.^(l+1)) .* (A.^(r+1)) .* (tauG.^(n+1)) .* tauC .* tauS)) .* (((i-T).*10^-3).^(p-r)).* exp(-A.*((i-T).*10^-3));
                            end
                            sum21(l+1) = sum(sum11) + (factorial(n) ./ ((B.^(k+1)).* (g.^(l+1)) .* (A.^(p+1)).* (tauG.^(n+1)) .* tauC .*tauS));
                            sum22(l+1) = sum(sum12) + (factorial(n) ./ ((B.^(k+1)).* (g.^(l+1)) .* (A.^(p+1)).* (tauG.^(n+1)) .* tauC .*tauS));
                        end
                        sum31(k+1) = sum(sum21) - (factorial(n) ./ ((B.^(k+1)).* (g.^(m+1)).* (tauG.^(n+1)) .* tauC)) .*(1 - exp(-(i.*10^-3)./tauS));
                        sum32(k+1) = sum(sum22) - (factorial(n) ./ ((B.^(k+1)).* (g.^(m+1)).* (tauG.^(n+1)) .* tauC)) .*(1 - exp(-((i-T).*10^-3)./tauS));
                  end
                  Bt1 = sum(sum31) + (factorial(n)./((B.^(n+1)).* phi.* (tauG.^(n+1)) .* tauC)).*(1 - exp(-(i.*10^-3)./tauS)) - (factorial(n)./((B.^(n+1)).* phi .* (tauG.^(n+1)) .*tauS)).*(1 - exp(-(i.*10^-3)./tauC));
                  Bt2 = sum(sum32) + (factorial(n)./((B.^(n+1)).* phi.* (tauG.^(n+1)) .* tauC)).*(1 - exp(-((i-T).*10^-3)./tauS)) - (factorial(n)./((B.^(n+1)).* phi .* (tauG.^(n+1)) .*tauS)).*(1 - exp(-((i-T).*10^-3)./tauC));
                  Bt(i+1) = Bt1 - Bt2; 
             end
        end
        R_S = At - wC.*Bt ;
        % % figure;
        % % plot(R_S);
        % % title('Surround Filter');
end
