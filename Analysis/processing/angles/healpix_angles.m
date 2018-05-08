function [  ] = healpix_angles( N_side )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab','/home/asus');

% addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab');
tic;

addpath('/home/asus/Programs/s2let/src/main/matlab');



% %prelocate memory not recomended, because if one uses a function we endup
% alocating twice
% 
% thetas=zeros(N_side);
% 
% phis=zeros(N_side);


[thetas, phis] = s2let_hpx_sampling_ring(N_side);

% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_t.cvs'),thetas,'single');
% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_p.cvs'),phis,'single');

% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_t.cvs'),thetas);
% csvwrite(strcat('/home/asus/Documents/angles',num2str(N_side),'_p.cvs'),phis);


toc;
end

            