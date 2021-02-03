function [outputArg1,outputArg2] = Test_healpixRange()

NSIDE=8;

aperture_angle=pi/4;

% angles_hpx(1,:)  = dlmread(strcat('../../../../python/angles',num2str(NSIDE),'_t.cvs'));
% angles_hpx(2,:) = dlmread(strcat('../../../../python/angles',num2str(NSIDE),'_p.cvs'));

%molweide projection
    

addpath(genpath('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/Analysis/wake_detection/slices_curv_2d'));

[thetas, phis] = s2let_hpx_sampling_ring(NSIDE);

phi_range_low=-aperture_angle+pi/2;
phi_range_high=aperture_angle+pi/2;

signal=double(thetas>=phi_range_low&thetas<=phi_range_high);

% signal=zeros(size(thetas))


% signal=zeros(size(thetas));
% signal=rand(size(thetas));

[x, y] = ssht_mollweide(thetas, phis,0,0);

fig=figure;
gridDelaunay = delaunay(x,y);
h = trisurf(gridDelaunay,x,y,signal(:)*0.0,signal(:));
set(h, 'LineStyle', 'none')
axis equal
axis off
campos([0 0 1])
camup([0 1 0])

end

