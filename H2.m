%% Jordan decomposition: Hf = H1-H2

function y = H2(psi)
if(psi<=0)
    y = 0;
else
    y = H1(psi)-Hf(psi);
end