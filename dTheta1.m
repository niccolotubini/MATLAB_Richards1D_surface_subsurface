%% Jordan decomposition: dTheta1

function y=dTheta1(psi)
global psic psic1 psic2 psic3 model

switch model
    
    %% Van Genuchten
    case 0
        if (psi<=psic)
            %left of critical value, take the original derivative
            y = dTheta(psi);
        else
            %on the right of the critical value, keep the maximum derivative
            y = dTheta(psic);
        end
        
        %% Romano
    case 1
        if (psi<=psic1)
            y = dTheta(psi);
        elseif(psic1<psi && psi<psic3)
            y = dTheta(psic1);
        elseif(psic3<=psi && psi<=psic2)
            y = dTheta(psi) + dTheta(psic1);
            %    y = dTheta(psi) + dTheta(psic1) - dTheta(psic3);
        else
            y = dTheta(psic2) + dTheta(psic1);
            %     y = dTheta(psic2) + dTheta(psic1) - dTheta(psic3);
        end
end

