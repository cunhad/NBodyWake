function [  ] = big_arrays( N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


array=zeros(N,1);

end

% big_arrays(3000000000)
% Error using zeros
% Requested 3000000000x1 (22.4GB) array exceeds maximum array size preference. Creation of arrays greater than this limit may take a long time and cause
% MATLAB to become unresponsive. See array size limit or preference panel for more information.

% big_arrays(2000000000)
% Error using zeros
% Out of memory. Type HELP MEMORY for your options.
% 
% Error in big_arrays (line 6)
% array=zeros(N,1);
%  

%no cluster com mem=64G deixa fazer ate N=70000000000 (70 bilhoes):
% >> a=zeros(70000000000,1);
% Error using zeros
% Requested 70000000000x1 (521.5GB) array exceeds maximum array size preference. Creation of arrays greater than this limit may take a long time and cause MATLAB to become unresponsive. See <a href="matlab:
% helpview([docroot '/matlab/helptargets.map'], 'matlab_env_workspace_prefs')">array size limit</a> or preference panel for more information.
