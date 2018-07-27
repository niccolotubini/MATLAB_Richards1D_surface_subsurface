global  alpha thetas thetar n m Ks psic K KL di IMAX dx dt h1m h2m sigma1 sigma2 w psic1 psic2 psic3

psic1   = h1m;
psic2   = -0.193552951216436;
psic3   = h2m;
h1m     = -0.125;
h2m     = -0.04;
sigma1  = 1;
sigma2  = 0.5;
w       = 0.4;
thetas  = 0.5;%0.41;         %[-] saturated water content
thetar  = 0.07;%0.095;
psi = linspace(-10,2,1000);

for i=1:1000
    h(i) = Hf(psi(i));
    h1(i) = H1(psi(i));
    h2(i) = H2(psi(i));
    dh(i) = dH(psi(i));
    dh1(i) = dH1(psi(i));
    dh2(i) = dH2(psi(i));
    
end
plot(psi,h2,'-')
%plot(psi,dtheta1,'-')
