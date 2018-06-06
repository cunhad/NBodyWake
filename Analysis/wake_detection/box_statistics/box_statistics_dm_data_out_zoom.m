function [  ] = box_statistics_dm_data_out_zoom( root,root_box_in,root_plot_out,root_snan_out,spec,aux_path,aux_path_box_in,aux_path_plot_out,aux_path_snan_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part_in,part_out,num_cores,level_window,dwbasis,angle_range,n_angles_1d,n_max)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%(example)  [  ] = box_statistics_dm_data_out_zoom('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast/','/home/asus/Dropbox/extras/storage/graham/small_res/box_plot/','/home/asus/Dropbox/extras/storage/graham/small_res/box_snan/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','','','10.000xv0.dat',2,1,[0,0,0],4,1,1,4,1,'sym6',0.01,20,1);

addpath(genpath('../../processing/'));


cd('../../preprocessing');

[~,redshift_list,~,~,nc,np,~,~,~,Gmu,~] = preprocessing_info(root,spec,aux_path );

bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];


[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);


path_data=strcat(strcat(root_box_in,spec,aux_path),'data/',aux_path_box_in,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/');

fileID = fopen(strcat(path_data,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_in),'_NSIDE',num2str(NSIDE),'.bin'));
display(strcat(path_data,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_in),'_NSIDE',num2str(NSIDE),'.bin'))

out_filtered_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
fclose(fileID);

% display(out_filtered_proj1d_angles(1,:))

[thetas, phis] = s2let_hpx_sampling_ring(NSIDE);

out_filtered_proj1d_angles_sort=flip(unique(sort(out_filtered_proj1d_angles(1,:),'descend')));

peaks_maxima=zeros(n_max);

for max_id=1:n_max
    peak=out_filtered_proj1d_angles_sort(max_id);
    display(peak)
    f_index_max=find(out_filtered_proj1d_angles==peak,1);
    phis_max=phis(f_index_max);
    thetas_max=thetas(f_index_max);
    display(phis_max)
    display(thetas_max)
    
    [xq,yq] = meshgrid(-angle_range/2:angle_range/n_angles_1d:angle_range/2, -angle_range/2:angle_range/n_angles_1d:angle_range/2);
    peaks_fine_filtered_proj1d_angles=zeros(size(xq));  %matrix with the peaks to be computed
    thetas_sq=zeros(size(xq));
    phis_sq=zeros(size(xq));
    [phis_q,thetas_q,~] = cart2sph(xq(:),yq(:),1);  %planar approximation
%     thetas_q=thetas_q-pi/2                          %elevation to theta
%     clearvars xq yq
    [rx(1,:),rx(2,:),rx(3,:)] = sph2cart(phis_q,thetas_q,1);           %points in the sphere
    
    %rotation to centralize on peak
    
    Ry = [cos(thetas_max) 0 -sin(thetas_max); 0 1 0; sin(thetas_max) 0 cos(thetas_max)];
    Rz = [cos(phis_max) -sin(phis_max) 0; sin(phis_max) cos(phis_max) 0; 0 0 1];
    
    rx=Rz*rx;
    rx=Ry*rx;
    
    
    [phis_q,thetas_q,~] = cart2sph(rx(1,:),rx(2,:),rx(3,:));
    thetas_q=thetas_q-pi/2;                      %elevation to theta

    
    proj1d_angles=zeros(length(bins)-1,length(phis_q));     %stores the plain 1d projections
    
    
    for part_id = 1  :   part_in
        
        [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,~, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part_in,part_id);
        
        Pos=mod(Pos,nc);
        
        Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
        Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
        Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
        
        
        lim_pre=(1/(lenght_factor))*nc;
        
        Pos(:,abs(Pos(1,:))>lim_pre|abs(Pos(2,:))>lim_pre|abs(Pos(3,:))>lim_pre)=[];
        
        histogr_1d_angles=zeros(length(bins)-1,length(phis_q));      %auxiliary variable that will stores all 1d projs
        
        
        parfor (angl_ind=1:length(phis_q),num_cores)
            
            theta=thetas_q(angl_ind);
            phi=phis_q(angl_ind);
            
            nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
            dz=nz*Pos;
            histogr_1d_angles(:,angl_ind)=histcounts(dz,bins);
            
            
        end
        
        proj1d_angles=proj1d_angles+histogr_1d_angles;
        
    end
    
    %this part only computes a normalization factor, since diferent
    %slices on the projection of the particles inside the cube have diferent
    %areas. The normalization factor is this area.
    
    [v f] = createCube; v = (v-[0.5 0.5 0.5])*nc ;
    
    
    parfor (angl_ind=1:length(phis_q),num_cores)
        theta=thetas_q(angl_ind);
        phi=phis_q(angl_ind);
        nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
        bin_sz=(length(bins)-1);
        area=zeros(bin_sz,1);
        for dist_1d=1:bin_sz
            plane = createPlane(nz*((bins(1,dist_1d+1)+bins(1,dist_1d))/2), nz);
            plane_overlap = intersectPlaneMesh(plane, v, f);
            area(dist_1d,1)=polygonArea3d(plane_overlap);
        end
        areas(:,angl_ind)=area(:,1);
    end
    
    %do the normalization
    
    proj1d_angles=proj1d_angles./areas;
    
    
    %here the wavelet filter takes place (on both fileds, the particle number and the density contrast)
    
    filtered_proj1d_angles=zeros(length(bins)-1,length(phis_q));
    level=floor(log2(length(bins)-1)); %level of the wavelet transformation
    
    parfor i=1:length(phis_q)
        
        proj1d=proj1d_angles(:,i);
        [dwt_proj1d,levels] = wavedec(proj1d,level,dwbasis);
        
        D=zeros(length(bins)-1,1);
        for lev_win = 1:length(level_window)
            lvwin=level_window(lev_win);
            D(:,1)=D(:,1)+wrcoef('d',dwt_proj1d,levels,dwbasis,lvwin);
        end
        filtered_proj1d_angles(:,i) = D(:,1);
        
    end
    
    %in this part the peak of each 1d projection is taked, Now we are
    %performing these operations on the filtered 1d projections
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
    max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
    
    for i=1:length(phis_q)
        
        peaks_fine_filtered_proj1d_angles(i)=max_amplitude_filtered_proj1d_angles(i);
        thetas_sq(i)=thetas_q(i);
        phis_sq(i)=phis_q(i);
    end
    
    %     fig0=figure('Visible', 'off');
    fig0=figure;
%     mesh(thetas_sq,phis_sq,peaks_fine_filtered_proj1d_angles);
        mesh(xq,yq,peaks_fine_filtered_proj1d_angles);
    tot_plot_path_out=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
    mkdir(tot_plot_path_out);
    mkdir(tot_plot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/'));
    saveas(fig0,strcat(tot_plot_path_out,'level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_peak_zoom_z',num2str(z),'_plot.png'));
%     saveas(fig0,strcat(tot_plot_path_out,'level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_peak_zoom_z',num2str(z),'_plot'));
    
    peaks_maxima(max_id)=max(peaks_fine_filtered_proj1d_angles(:));
    display(peaks_maxima)
    
end

snan_tot_path_out=strcat(root_snan_out,spec,aux_path,'snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
mkdir(snan_tot_path_out);
mkdir(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'));
snan=[max(peaks_maxima) find(peaks_maxima==max(peaks_maxima),1)];
dlmwrite(strcat(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_snan_box_zoom_z',num2str(z),'_data.txt'),snan,'delimiter','\t');    

cd('../wake_detection/box_statistics');


end

function [thetas, phis] = s2let_hpx_sampling_ring(nside)

% s2let_hpx_sampling_ring 
% Compute the Healpix sampling nodes
% Default usage :
%
%   [thetas, phis] = s2let_hpx_sampling_ring(nside)
%
% nside is the input Healpix resolution,
% thetas contains the colatitudes of the nodes,
% phis contains the longitudes of the nodes.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

%     npix = 12 * nside^2;
%     thetas = zeros(npix,1);
%     phis = zeros(npix,1);
%     for i = 1:npix
%        [t,p] = healpix_pix2ang_ring(i, nside);
%        thetas(i) = t;
%        phis(i) = p;
%     end

    npix = 12 * nside^2;
    ipix = 1:npix;
    nl2 = 2 * nside;
    ncap = 2  * nside * (nside - 1);  % points in each polar cap, =0 for nside =1

    thetas = zeros(npix,1);
    phis = zeros(npix,1);
    
    ind = (ipix <= ncap); % North polar cap
        hip   = ipix(ind) .* 0.5;
        fihip = fix(hip);
        iring = fix( sqrt( hip - sqrt(fihip) ) ) + 1 ;% counted from North pole
        iphi  = ipix(ind) - 2*iring.*(iring - 1) ;

        thetas(ind) = acos( 1.0 - iring.^2.0 ./ (3.0*nside^2.0) );
        phis(ind)   = (iphi - 0.5) .* pi ./ (2.0.*iring);

    ind = (ipix <= nl2*(5*nside+1) & ipix > ncap); % Equatorial region

       ip    = ipix(ind) - ncap - 1;
       nl4   = 4*nside;
       iring = floor( ip / nl4 ) + nside; % counted from North pole
       iphi  = mod(ip,nl4) + 1;

       fodd  = 0.5 * (1 + mod(iring+nside,2)) ; % 1 if iring+nside is odd, 1/2 otherwise
       thetas(ind) = acos( (nl2 - iring) ./ (1.5*nside) );
       phis(ind)   = (iphi - fodd) * pi /(2.0*nside);

    ind = (ipix > nl2*(5*nside+1) & ipix > ncap); % South polar cap

       ip    = npix - ipix(ind) + 1;
       hip   = ip*0.5;
       fihip = round(hip);
       iring = floor( sqrt( hip - sqrt(fihip) ) ) + 1 ; % counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring.*(iring-1));

       thetas(ind) = acos( -1.0 + iring.^2 ./ (3.0*nside^2) );
       phis(ind)   = (iphi - 0.5) * pi ./ (2.0*iring);

end

function [theta, phi] = healpix_pix2ang_ring(ipix, nside)

    npix = 12 * nside^2;
    nl2 = 2 * nside;
    ncap = 2  * nside * (nside - 1);  % points in each polar cap, =0 for nside =1

    if (ipix <= ncap)
        hip   = ipix * 0.5;
        fihip = fix(hip);
        iring = fix( sqrt( hip - sqrt(fihip) ) ) + 1 ;% counted from North pole
        iphi  = ipix - 2*iring*(iring - 1) ;

        theta = acos( 1.0 - iring^2.0 / (3.0*nside^2.0) );
        phi   = (iphi - 0.5) * pi/(2.0*iring);

    elseif (ipix <= nl2*(5*nside+1)) % Equatorial region ------

       ip    = ipix - ncap - 1;
       nl4   = 4*nside;
       iring = floor( ip / nl4 ) + nside; % counted from North pole
       iphi  = mod(ip,nl4) + 1;

       fodd  = 0.5 * (1 + mod(iring+nside,2)) ; % 1 if iring+nside is odd, 1/2 otherwise
       theta = acos( (nl2 - iring) / (1.5*nside) );
       phi   = (iphi - fodd) * pi /(2.0*nside);

    else % South Polar cap -----------------------------------

       ip    = npix - ipix + 1;
       hip   = ip*0.5;
       fihip = round(hip);
       iring = floor( sqrt( hip - sqrt(fihip) ) ) + 1 ; % counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

       theta = acos( -1.0 + iring^2 / (3.0*nside^2) );
       phi   = (iphi - 0.5) * pi/(2.0*iring);

    end

end

% 
% 
% function [x, y] = ssht_mollweide(thetas, phis,thetas_shift, phis_shift)
% % ssht_mollweide - Compute Mollweide projection
% %
% % Compute Mollweide projection of spherical coordinates.
% %
% % Usage is given by
% %
% %   [x,y] = ssht_mollweide(thetas, phis)
% %
% % where thetas and phis are spherical coordinates and x and y are the
% % projected Mollweide coordinates.
% %
% % Author: Jason McEwen (www.jasonmcewen.org)
% 
% MAX_ITERATIONS = 1e5;
% TOL = 1e-10;
% 
% 
% % Convert theta to longitude.
% thetas = pi/2 - thetas;
% phis = phis - pi;
% 
% % 
% [rx(1,:),rx(2,:),rx(3,:)] = sph2cart(phis,thetas,1);
% 
% theta=pi/2+thetas_shift;
% phi=-phis_shift;
% 
% Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
% 
% rx=Rz*rx;
% rx=Ry*rx;
% 
% [phis,thetas,~] = cart2sph(rx(1,:),rx(2,:),rx(3,:));
% 
% 
% t = thetas;
% for it = 1:MAX_ITERATIONS
% 
%    dt = (t + sin(t) - pi.*sin(thetas)) ./ (1 + cos(t));
%    t = t - dt;
%    
%    if(max(abs(dt)) < TOL)
%       break;
%    end
%    
% end
% t = t/2;
% x = 2 .* sqrt(2) ./ pi .* phis .* cos(t);
% y = sqrt(2) .* sin(t);
% end
% 
