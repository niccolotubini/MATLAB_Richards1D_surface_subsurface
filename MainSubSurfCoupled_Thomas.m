%% BTCS scheme for the 1D Richards equation using the nested Newton method of
% Casulli & Zanolli.
% coupling surface-subsurface flow
% F. Gugole, 'A fast semi-implicit 3D algorithm for the solution of coupled
%  free-surface and variably saturated sub-surface flows', 2017


%% SWRC models implemented
% 0 Van Genuchten
% 1 Romano

format long
clear all
close all
clc


global  alpha thetas thetar n m Ks psic K KL di IMAX dx dt h1m h2m sigma1 sigma2 w psic1 psic2 psic3 aa bb model
model = 0;
aa=0;%10^(-7);
bb=0;%4.4*10^(-10);
%Phisical model parameters in SI units
%day     = %24*3600;
Ks      = 0.00003697;%0.062/day;    %[meter/second]
thetas  = 0.5;%0.41;         %[-] saturated water content
thetar  = 0.07;%0.095;        %[-] residuel water content
n       = 1.16;%1.31;         %[-] parameter n
m       = 1-1/n;        %[-] parameter m
alpha   = 5.88;%1.9;          %[m^(-1)]
psic    = -1/alpha*((n-1)/n)^(1/n);  %critical value of psi where the maximum of the derivative is located

h1m     = -1.25;
h2m     = -0.04;
sigma1  = 1;
sigma2  = 0.5;
w       = 0.4;
psic1   = -0.4598493014643029;
psic3   = -0.19355295121650085;
psic2   = -0.031152031322856197;

%Domain
zL = 0;                 % bottom
zR = 2;                 % surface
IMAX = 300;             % number of control volumes
dx = (zR-zL)/IMAX;      % mesh spacing
z = linspace(zL+dx/2,zR-dx/2,IMAX);
z = [z,2];
tend =3e8;              % set the final simulation time
time = 0;               % initial time

% initialize variables
theta = zeros(1,IMAX);
volume = zeros(1,IMAX+1);
volumeNew = zeros(1,IMAX+1);
psi = zeros(1,IMAX+1);
a = zeros(1,IMAX+1);
b = zeros(1,IMAX+1);
c = zeros(1,IMAX+1);
rhs = zeros(1,IMAX+1);
K = zeros(1,IMAX+1);
f  = zeros(1,IMAX+1);
fk = zeros(1,IMAX+1);

%set the initial condition for psi both in soil and surface
for i=1:IMAX+1
    if(i==IMAX+1)
        psi(i) = 0.3;%psi(IMAX)-dx/2;   % psi at surface: if >0 ponding, otherwise not
    else
        psi(i) = -z(i);     % hydrostatic pressure
    end
end
psi0 = psi;
% time cycle
NMAX=90;
for nit=1:20
    disp(sprintf('time iteration:%d', nit ));
    dt = 10;   %BTCS is unconditionally stable
    if(time+dt>tend)
        dt=tend-time;
    end
    
    % Right boundary condition
    %     if time< 15000
    %         Rain = 0.0006/300; % rainfall rate m/s
    %     else
    %         Rain=0;
    %     end
    Rain = 0/1000/300;
    % Left boundary condition
    psiL = 0.0;
    KL = kappa(psiL);
    
    for i=1:IMAX+1
        if(i==IMAX+1)
            volume(i) = Hf(psi(i));
            K(i) = kappa(psi(i)); % hydaraulic conductivity at soil surface
        else
            volume(i) = Thetaf(psi(i))*dx;
            theta(i) = Thetaf(psi(i));
            K(i) = kappa(psi(i));
        end
    end
    
    % plot at time level n
    subplot(2,1,1)
    plot(z,psi,'o')
    ylim([-2,0.9])
    ylabel('psi [L]')
    subplot(2,1,2)
    plot(z(1:IMAX),theta,'o')
    ylabel('volume [L]')
    xlabel('z [m]')
    title(sprintf('Current time = %f',time))
    drawnow
    
    if(time>=tend)
        break
    end
    
    %compute right hand side and the linear part of the system
    for i=1:IMAX+1
        if(i==1)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+KL);
            %Kp = max(K(i),K(i+1));
            %Km = max(K(i),KL);
            a(i) = 0;                   % lower diagonal -> moved to the rhs multiplied by (-2)*psiL
            b(i) = +(2*Km+Kp)*dt/dx;  % diagonal(not rhs!)
            c(i) = -Kp*dt/dx;         % upper diagonal
            rhs(i) = volume(i) + dt*(Kp-Km) + 2*Km*dt/dx*psiL;  % right hand side <- the "b" of lecture notes
        elseif(i==IMAX+1)
            Kp = 0;
            Km = 0.5*(K(i)+K(i-1));
            %Kp = 0;
            %Km = max(K(i),K(i-1));
            %rainfall rate (Rain) boundary condition
            a(i) = - dt/dx*Km*2;
            b(i) = + dt/dx*(2*Km+2*Kp);
            c(i) = - 2*dt/dx*Kp;
            rhs(i) = volume(i) + dt*(Rain-Km);
        % questo if serve perche` a differenza del codice java in cui ho le
        % distanze dei centroidi qui ho solo un dx costante
        elseif(i==IMAX)
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            a(i) = - dt/dx*Km;
            b(i) = + dt/dx*(Km+2*Kp);
            c(i) = - 2*dt/dx*Kp;
            rhs(i) = volume(i) + dt*(Kp-Km);    
        else
            Kp = 0.5*(K(i)+K(i+1));
            Km = 0.5*(K(i)+K(i-1));
            a(i) = -Km*dt/dx;         % lower diagonal
            b(i) = +(Km+Kp)*dt/dx;    % diagonal(not rhs!)
            c(i) = -Kp*dt/dx;         % upper diagonal
            rhs(i) = volume(i) + dt*(Kp-Km);  % right hand side <- the "b" of lecture notes
        end
    end
    
    
    %% Newton method
    tol = 10e-12;
    
    % initial guess
    psi(IMAX+1) = max(psi(IMAX+1),0.1);
    psi(1:IMAX) = min(psi(1:IMAX),psic);%psic1

    for iNewton=1:100 %outer Newton iterations
        
        di = zeros(1,IMAX+1);
       
        for i=1:IMAX+1
            if(i==1)
                f(i) = Thetaf(psi(i))*dx + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
            elseif(i==IMAX+1)
                f(i) = Hf(psi(i)) + a(i)*psi(i-1) + b(i)*psi(i) - rhs(i);
            else
                f(i) = Thetaf(psi(i))*dx + a(i)*psi(i-1) + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
            end
        end
        %f(:)
        outres = sqrt(sum(f.*f)); %outer residual
        disp(sprintf('    Outer iteration %d, outres=%e',iNewton, outres));
        if(outres<tol)
            break %tolerance has been reached
        end
        
        psik = psi;  % save the value at the current outer iteration
        %psi(IMAX+1) = max(psi(IMAX+1),1);
        %psi(1:IMAX) = max(psi(1:IMAX),psic2);%psic2 % initial guess for the inner iteration

        for inner=1:100 % inner Newton step
            
            di = zeros(1,IMAX+1);
            
            for i=1:IMAX+1
                if(i==1)
                    fk(i) = Theta1(psi(i))*dx - ( Theta2(psik(i))*dx + dTheta2(psik(i))*dx*(psi(i)-psik(i)) ) + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
                    di(i) = dTheta1(psi(i))*dx - dTheta2(psik(i))*dx;
                elseif(i==IMAX+1)
                    fk(i) = H1(psi(i)) - ( H2(psik(i)) + dH2(psik(i))*(psi(i)-psik(i)) ) + a(i)*psi(i-1) + b(i)*psi(i) - rhs(i);
                    di(i) = dH1(psi(i)) - dH2(psik(i));
                else
                    fk(i) = Theta1(psi(i))*dx - ( Theta2(psik(i))*dx + dTheta2(psik(i))*dx*(psi(i)-psik(i)) ) + a(i)*psi(i-1) + b(i)*psi(i) + c(i)*psi(i+1) - rhs(i);
                    di(i) = dTheta1(psi(i))*dx - dTheta2(psik(i))*dx;
                end
            end

            inres=sqrt(sum(fk.*fk));
            disp(sprintf('        Inner iteration %d, inres= %e', inner,inres));
            if(inres<tol)
                break
            end

            dpsi = Thomas(a,b+di,c,fk);
            psi(:) = psi(:) - dpsi(:);  %update psi at the inner iteration
            
        end
        
        
    end
    PSI(nit,:)=psi;
    for i=1:IMAX+1
        if(i==IMAX+1)
            volumeNew(i) = Hf(psi(i));
        else
            volumeNew(i) = Thetaf(psi(i))*dx;
        end
    end
    for i=1:IMAX+1
        if i==0
            kM=0.5*(KL+K(i));
            v(i)=-kM* (psi(1)-psiL)/(dx/2) -kM;
        elseif i==IMAX+1
            kP=0.5*(K(i)+K(i));
            v(i)=-kP* (psi(i)-psi(i-1))/(dx/2) -kP;
        else
            kP=0.5*(K(i)+K(i+1));
            v(i)=-kP* (psi(i+1)-psi(i))/dx -kP;
        end
    end
    psi(IMAX+1);
    psi(IMAX);
    %errorMass(nit)= sum(volumeNew) - sum(volume) - dt*( -0.5*(K(IMAX+1)+K(IMAX+1))*( (psi(IMAX+1)-psi(IMAX) )/(dx/2) +1 ) ...
    %                                                                         + 0.5*( KL+K(1) )*( (psi(1)-psiL)/(dx/2) +1 ) );
    errorMass1(nit)= sum(volumeNew) - sum(volume) - dt*( +Rain - 0.5*( KL+K(1) )*( (psi(1)-psiL)/(dx/2) +1 ) );
    %disp(sprintf('    Error mass:%e', errorMass(nit) ));
    disp(sprintf('    Error mass1:%e', errorMass1(nit) ));
    time = time+dt;     %advance time
end








