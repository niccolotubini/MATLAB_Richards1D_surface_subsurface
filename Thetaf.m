%% SWRC

function y=Thetaf(psi)
global alpha thetas thetar n m aa bb h1m h2m sigma1 sigma2 w model

switch model
    
    %% Van Genuchten
    case 0
        if (psi<=0)
            y = thetar + (thetas-thetar)/( (1+abs(alpha*psi)^n)^m );
        else
            %y = thetas + 9.81*( aa + thetas*bb )*psi;
            y = thetas;
    
        end
    
    %% Romano    
    case 1
        if (psi<0)
            y = thetar + (thetas-thetar)*( w/2*erfc( log(psi/h1m)/(sigma1*sqrt(2)) ) + (1-w)/2*erfc( log(psi/h2m)/(sigma2*sqrt(2)) ) );
        else
            y = thetas;
        end
        
end

