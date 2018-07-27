%matrix-free conjugate method for the solution of A*x=b with A symmetric positive definite
%Input: (we don't need the matrix as input)
%       b = right hand side
%output: 
%       x = solution vector
%the user has to provide a function "matop(x)" wich computes the matrix-vector product A*x
% CGop � lo stesso sia per 1D che 2D quello che cambia � la funzione matop

function x=CGop(b)
x = b;                %initial guess
N = numel(b);         %get the problem size 
r = b-matopCoupled(x);     %residual
p = r;                %initial search direction
alpha = sum(sum(r.*r));  %square of r (scalar product) %or r'*r; 
for k=1:4*N           %2* or 4*N to have more precision
    if (sqrt(alpha)<1e-14)   %nominal machine operation is 1e-16 so we use something smaller
        %norm of the residual is sufficently small
        break
    end
    v = matopCoupled(p);      
    lambda = alpha/sum(sum(p.*v));
    x = x + lambda*p;   
    r = r - lambda*v;
    alphak = alpha;
    alpha = sum(sum(r.*r));
    p = r + alpha/alphak*p;   
end


