%% Hydrayulic conductivity: Van Genucthen parametrization following Mualem's assumpution
%% vg
function K=kappa(psi)
global  thetas thetar m Ks model sigma1 sigma2 h1m h2m w

switch model
    
    %% Van Genuchten
    case 0
        temp = (Thetaf(psi)-thetar)/(thetas-thetar);
        
        if(temp<1)
            K = Ks*temp^(1/2)*(1-(1-temp^(1/m))^m)^2;
        else
            K = Ks;
        end
        
        %% Romano
    case 1
        if(psi<0)
            sat = (Thetaf(psi)-thetar)/(thetas-thetar); %saturation is between 0 1
            aa = (sigma1^2+log(psi/h1m))/(sigma1*sqrt(2));
            bb = (sigma2^2+log(h1m/h2m*psi/h1m))/(sigma2*sqrt(2));
            r = h1m/h2m*(1-w)/w*exp(0.5*(sigma1^2-sigma2^2));
            K = Ks*sqrt(sat)*( 0.5*erfc(aa)/(1+r) +  0.5*erfc(bb)/(1+r) )^2;
        else
            K = Ks;
        end
end