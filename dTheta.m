%% Derivative of Thetaf w.r.t. psi

function y=dTheta(psi)
global alpha thetas thetar n m  aa bb h1m h2m sigma1 sigma2 w model

switch model
    
    %% Van Genuchten
    case 0
        if (psi<=0)
            %unsaturated medium
            y = alpha*n*m*(thetas-thetar)/((1+abs(alpha*psi)^n)^(m+1)) *abs(alpha*psi)^(n-1);
        else
            %saturated medium
            %y = 9.81*( aa + thetas*bb );
            y = 0.0;
        end
        
        %% Romano
    case 1
        if (psi<0)
            gamma1 = exp( -( log(psi/h1m)/(sigma1*sqrt(2)) )^2 );
            gamma2 = exp( -( log(psi/h2m)/(sigma2*sqrt(2)) )^2 );
            %unsaturated medium
            y = 1/(sqrt(2*pi)*psi/h1m)*(thetas-thetar)*( w/sigma1*gamma1 + (1-w)/sigma2*gamma2);
        else
            %saturated medium
            y = 0;
        end
        
end
