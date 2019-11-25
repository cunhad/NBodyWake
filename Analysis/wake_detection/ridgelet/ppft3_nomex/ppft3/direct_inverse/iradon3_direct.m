% function Y = iRadond3d(res,EPS)
%
% function iRadon3d restores the initial image from its 3d-radon transformaton. In fact it only turns Radon transformation
% to ppfft and then applyies ippfd3d to it.
%
% Oren 11/12/2005
function Y = iradon3_direct(res,EPS);

temp = cfftd(res,2);
Y=ippft3_direct(temp,EPS);