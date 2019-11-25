% function x=adjF3d(y)
%
% Adjoint of the operator F for 3d case.
%
% Parameters:
%      y        vector of length n+1 (n even)
%      x        vector of length n
%
% Oren  11/12/2005
function x=adjF3d(y)

% Check that the vector y is of length n+1 with n even.
n = length(y)-1;
if mod(n,2)~=0
    error('y must be of length n+1 (n even)');
end
m=3*n+1;

y2 = zeros(1,m);
y2(1:3:m) = y;
x = m.*icfft(y2);
x = x(n+1:2*n);
