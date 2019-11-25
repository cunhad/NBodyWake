%[fk,ier]=nufft1d1(length(y),y,u,-1,1e-5,8); fk*length(y)
function [out, ier]=ToeplitzResamp(f,y,x,n, EPS)
% The function resamples a trigonometric polynomial using fast inversion of
% Toeplitz matrices. Suppose A is a non uniform Fourier matrix, we wish to
% find x that minimizes ||Ax-b||. Instead of solving it directly using LS,
% We do the following.
% Minimizing ||Ax-b|| is equal to solving A'*A*x=A'*b.
% Note that A'*A is a square symmetric Toeplitz matrix. To find the
% Toeplitz vector it is enough to compute A'*A(:,1), this is equivalent to
% applying non equally space Fourier transform to the column A(:,1). This
% column is a complex exponent with fixed frequency: exp(1i*y*k(1))
%f, y, x - column vectors
persistent yy;
persistent D1 D2 D3 D4;
global preptime;

f=f(:);
y=y(:);
x=x(:);
%eps=1e-15;
k = [-n/2:n/2-1];
if length(yy)==length(y)
    err=norm(y-yy);
else
    err=1;
end
%fprintf('%d, %d \n',length(y), length(x))
if length(yy)~=length(y)
    h=tic;
    
    % The line below compute A'*A(:,1). A(:,1) is a vector of complex
    % exponent with frequency k(1) sampled at points y. The multiplication
    % is done using the nufft1d1.
    [c,ier1]=nufft1d1(length(y),y,exp(1i*y*k(1)),-1,EPS,n);
    c=c*length(y);
    % Now c is our A'*A(:,1)
    [m1,m2,m3,m4]=topinv(c,c); % Compute the inverse using Gohberg-Semencul factorization
    [D1,D2,D3,D4]=topprepinv(m1,m2,m3,m4);
    % Here we inverted A'*A, which is a symmetric Toeplitz matrix represented by c.
    
    preptime=preptime+toc(h);
end
   % Here we compute A'*b, with f uses as b. This is done again with the
   % nufft since A is a non equally spaced Fourier matrix.
   [v,ier2]=nufft1d1(length(y),y,f,-1,EPS,n);
   v = v*length(y);
   % Now v is our A'*f
   yy=y;
   alpha=topinvmul(D1,D2,D3,D4,v); % Multiply v with the inverse of the Toeplitz.
   % alpha is the solution for the LS problem, alpha=inv(A'*A)*A'*b

   % Now compute the values of the trigonometric polynomial, whose
   % coefficients are given by the vector alpha at points x. These values
   % are given by A*x, which is computed with the nufft1d2 function, again,
   % since A is a non uniform Fourier matrix.
   [out, ier3] = nufft1d2(length(x),x,1,EPS,n,alpha);
%ier=ier1+ier2+ier3;

end
