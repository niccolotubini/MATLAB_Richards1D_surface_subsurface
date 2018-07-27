global  alpha thetas thetar n m Ks psic K KL di IMAX dx dt h1m h2m sigma1 sigma2 w psic1 psic2 psic3 aa bb

h1m     = -1.25;
h2m     = -0.04;
sigma1  = 1;
sigma2  = 0.5;
w       = 0.4;
psic1   = -0.4598493014643029;
psic3   = -0.19355295121650085;
psic2   = -0.031152031322856197;
thetas  = 0.5;%0.41;         %[-] saturated water content
thetar  = 0.02;%0.095;
alpha =1.47;
n=1.7;
m=1-1/n;
aa = 10^(-7);
bb = 4.4*10^(-10);
psi = linspace(-10,-0.001,1000);

for i=1:1000
    theta(i) = Thetaf(psi(i));
    theta1(i) = Theta1(psi(i));
    theta2(i) = Theta2(psi(i));
    dtheta(i) = dTheta(psi(i));
    dtheta1(i) = dTheta1(psi(i));
    dtheta2(i) = dTheta2(psi(i));
    
end
%plot(psi,dtheta-(dtheta1-dtheta2),'-')
plot(psi,dtheta1,'-')

% plot(psi,theta-(theta1-theta2),'-')
%plot(psi,theta,'-')
for i=1:999
    diff1(i) = dtheta1(i+1)-dtheta1(i);
    diff2(i) = dtheta2(i+1)-dtheta2(i);
end