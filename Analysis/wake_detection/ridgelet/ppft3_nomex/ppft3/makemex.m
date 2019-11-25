% Compile the required mex files.
%
% Yoel Shkolnisky 20/05/2013

cd conv_inverse
mex -O nufft3dauxmx.cpp
cd ..