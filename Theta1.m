%% Jordan decomposition: Theta1

function y = Theta1(psi)
global psic psic1 psic2 psic3 model

switch model
    
    %% Van Genuchten
    case 0
        if (psi<=psic)
            y = Thetaf(psi);
        else
            y = Thetaf(psic) + dTheta(psic)*(psi-psic);
        end
        
        %% Romano
    case 1
        if (psi<=psic1)
            y = Thetaf(psi);
        elseif (psic1<psi && psi<=psic3)
            y = Thetaf(psic1) + dTheta(psic1)*(psi-psic1);
        elseif (psic3<psi && psi<=psic2)
            y = Thetaf(psi) - Thetaf(psic3) + Thetaf(psic1) +dTheta(psic1)*(psi-psic1);
            %     y = Thetaf(psi) - Thetaf(psic3) + Thetaf(psic1) + dTheta(psic1)*(2*psi-psic1-psic3) - dTheta(psic3)*(psi-psic3);
        else
            y = Thetaf(psic1) + Thetaf(psic2) - Thetaf(psic3) + dTheta(psic2)*(psi-psic2) + dTheta(psic1)*(psi-psic1);
            %     y = Thetaf(psic1) + Thetaf(psic2) - Thetaf(psic3) + dTheta(psic1)*(2*psic2-psic1-psic3) - dTheta(psic3)*(psic2-psic3) + (dTheta(psic2)+dTheta(psic1)-dTheta(psic3))*(psi-psic2);
        end
        
end
