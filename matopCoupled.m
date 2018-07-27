%% matop computes the product between the coefficient matrix and the unknown vector for Richards equation coupled with surface flow

function Apsi = matopCoupled(psi)
global di dx  dt K IMAX  KL  
Apsi = di.*psi;
for i=1:IMAX+1
    
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
%             Kp = max(K(i),K(i+1));
%             Km = max(K(i),KL);
            Apsi(i) =  Apsi(i)+ dt/dx*(2*Km+Kp)*psi(i) - dt/dx*Kp*psi(i+1);
            upper(i) =- dt/dx*Kp;
            main(i) = + dt/dx*(2*Km+Kp);
            lower(i) =0;
        elseif(i==IMAX+1)
            Km = 0.5*(K(i) + K(i-1));
%             Km = max(K(i),K(i-1));
            % Neumann boundary condition: fluxTop(t) = Rain(t)
            Apsi(i) =  Apsi(i) + 2*dt/dx*(Km)*psi(i) - 2*dt/dx*Km*psi(i-1);
            upper(i) =0;
            main(i) = + 2*dt/dx*(Km);
            lower(i) = - 2*dt/dx*Km;
            %Apsi(i) =  Apsi(i);
        elseif(i==IMAX)
            Kp = 0.5*(K(i) + K(i+1));
            Km = 0.5*(K(i) + K(i-1));
%             Kp = max(K(i),K(i+1));
%             Km = max(K(i),K(i-1));
            Apsi(i) =  Apsi(i) - 2*dt/dx*Kp*psi(i+1) + dt/dx*(2*Kp+Km)*psi(i) - dt/dx*Km*psi(i-1);
            upper(i) =- 2*dt/dx*Kp;
            main(i) = + dt/dx*(Km+2*Kp);
            lower(i) = - dt/dx*Km;
        else
            Kp = 0.5*(K(i) + K(i+1));
            Km = 0.5*(K(i) + K(i-1));
            Apsi(i) =  Apsi(i) - dt/dx*Kp*psi(i+1) + dt/dx*(Km+Kp)*psi(i) - dt/dx*Km*psi(i-1);
            upper(i) = - dt/dx*Kp;
            main(i) = + dt/dx*(Km+Kp);
            lower(i) = - dt/dx*Km;
            
        end
end
