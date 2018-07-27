psi = linspace(-10,2,1000);

for i=1:1000
    h(i) = Hf(psi(i));
    h1(i) = H1(psi(i));
    dh(i) = dH(psi(i));
    dh1(i) = dH1(psi(i));
    
end
plot(psi,h,'-')
