function [ out_filtered_proj1d_angles ] = box_statistics_dm_analysis_dwt_peaks_s2let( root,root_box_in,root_plot_out,root_snan_out,spec,aux_path,aux_path_box_in,aux_path_plot_out,aux_path_snan_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part_in,part_out,num_cores,level_window,dwbasis,angle_range,n_angles_1d,n_max)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%(example)  [  ] = box_statistics_dm_analysis_dwt_peaks('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast/','/home/asus/Dropbox/extras/storage/graham/small_res/box_plot/','/home/asus/Dropbox/extras/storage/graham/small_res/box_snan/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','','','10.000xv0.dat',2,1,[0,0,0],4,1,1,4,1,'sym6',0.01,20,1);
%(example)  [  ] = box_statistics_dm_analysis_dwt_peaks('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast/','/home/asus/Dropbox/extras/storage/graham/small_res/box_plot/','/home/asus/Dropbox/extras/storage/graham/small_res/box_snan/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','','','10.000xv0.dat',2,1,[0,0,0],64,1,1,4,1,'sym6',0.5,20,3); 
%(example)  [  ] = box_statistics_dm_analysis_dwt_peaks('/home/asus/Dropbox/extras/storage/graham/', '/home/asus/Dropbox/extras/storage/graham/box_stat_cubic_fast_ap/','/home/asus/Dropbox/extras/storage/graham/box_plot_/','/home/asus/Dropbox/extras/storage/graham/box_snan_/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample2001/','','','','10.000xv0.dat',2,1,[0,0,0],512,16,1,4,1,'sym6',0.5,0,3);

%(example)  [  ] = box_statistics_dm_analysis_dwt_peaks('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast_ap_cic/','/home/asus/Dropbox/extras/storage/graham/small_res/box_plot_cic/','/home/asus/Dropbox/extras/storage/graham/small_res/box_snan_cic/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','','','10.000xv0.dat',2,1,[0,0,0],64,1,1,4,1,'sym6',0.5,20,3); 

%(example)  [  ] = box_statistics_dm_analysis_dwt_peaks_s2let('/home/asus/Dropbox/extras/storage/guillimin/', '/home/asus/Dropbox/extras/storage/guillimin/box_stat_cubic_fast_ap_cic/','/home/asus/Dropbox/extras/storage/guillimin/box_plot_cic_test/','/home/asus/Dropbox/extras/storage/guillimin/box_snan_cic_test/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0003/','','','','10.000xv0.dat',2,1,[0,0,0],512,16,4,1,1,'sym6',0.5,20,3);

% 
% myCluster = parcluster('local');
% myCluster.NumWorkers=num_cores;
% saveProfile(myCluster);
% 
% p = parpool(num_cores);

addpath(genpath('../../processing/'));

% addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab','/home/asus');

addpath('/home/asus/Programs/s2let/src/main/matlab');


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

sz = size(out_filtered_proj1d_angles(1,:));
nsideguessed = sqrt(max(sz)/12);
    L = 2*nsideguessed;
    B=2;
    J_min=0;
    
    [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(out_filtered_proj1d_angles(1,:),'B',B,'L',L,'J_min',J_min);
    
    % Plot
J = s2let_jmax(L, B);


for j = J_min:J-1
    f_wav{j+1,1}(1,:)=0;
end
% 
% for j = J:J
%     f_wav{j+1,1}(1,:)=0;
% end

out_filtered_proj1d_angles(1,:)=s2let_transform_axisym_synthesis_hpx(f_wav, f_scal, 'B',B,'L',L,'J_min',J_min);
% out_filtered_proj1d_angles(1,:)=f_wav{J-J_min+1}+f_wav{J-J_min};



% out_filtered_proj1d_angles_sort=flip(unique(sort(out_filtered_proj1d_angles(1,:),'descend')));
[f_index_max] = find_peaks(thetas, phis,out_filtered_proj1d_angles(1,:),NSIDE,n_max);

out_filtered_proj1d_angles_sort=out_filtered_proj1d_angles(1,f_index_max);
% phis_max=phis(f_index_max(1));
% thetas_max=thetas(f_index_max(1));



tot_plot_path_out=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/');
mkdir(tot_plot_path_out);

tot_plot_path_out_zoom=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
mkdir(tot_plot_path_out_zoom);

snan_tot_path_out=strcat(root_snan_out,spec,aux_path,'snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/',strcat('level_window',mat2str(level_window(:)),'/peak/'));
snan_tot_path_out_zoom=strcat(root_snan_out,spec,aux_path,'snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/',strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'),'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');

mkdir(snan_tot_path_out);
mkdir(snan_tot_path_out_zoom);

% mkdir(snan_tot_path_out_zoom,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'));
% mkdir(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak/'));


peaks_maxima=zeros(n_max);

for max_id=1:n_max
    peak=out_filtered_proj1d_angles_sort(max_id);
    display(peak)
%     f_index_max=find(out_filtered_proj1d_angles(1,:)==peak,1);
%     phis_max=phis(f_index_max);
%     thetas_max=thetas(f_index_max);
    
%       f_index_max=find(out_filtered_proj1d_angles(1,:)==peak,1);
     phis_max=phis(f_index_max(max_id));
     thetas_max=thetas(f_index_max(max_id));
    display(phis_max)
    display(thetas_max)
    
    %plot the spherical projection
    
    [x, y] = ssht_mollweide(thetas, phis,0,0);        
    
%     fig1=figure('Visible', 'off');
    fig1=figure;
    set(gcf, 'Position', [0 0 1200 600]);
    ax1 = axes('Position',[0.05 0.13 0.7 0.7]);
    
    gridDelaunay = delaunay(x,y);
    h = trisurf(gridDelaunay,x,y,out_filtered_proj1d_angles(1,:)*0.0,out_filtered_proj1d_angles(1,:));
    colorbar('southoutside');
    
%     display(out_filtered_proj1d_angles(1,:));
    
    set(h, 'LineStyle', 'none')
    
    hold on;
    
    axis equal
    axis off
    campos([0 0 1])
    camup([0 1 0])
    
    [x_max, y_max] = ssht_mollweide(thetas_max, phis_max,0,0);       
    scatter(x_max, y_max,255,'r');
    
    descr = {strcat('max = ',num2str(peak));
    strcat('std = ',num2str(std(out_filtered_proj1d_angles(1,:))));
    strcat('stn = ',num2str(peak/std(out_filtered_proj1d_angles(1,:))))};
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    txt=text(0.75,0.5,descr);
    set(txt,'Parent',ax1,'interpreter', 'latex');
    
    hold off;
    

%     mkdir(tot_plot_path_out,strcat('/peak/','max_',num2str(max_id),'/'));
%     saveas(fig1,strcat(tot_plot_path_out,'/peak/','max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_peak_z',num2str(z),'_plot.png'));
 
%     mkdir(snan_tot_path_out,strcat('max_',num2str(max_id),'/'));
%      snan_zoom=[max(peaks_maxima) find(peaks_maxima==max(peaks_maxima),1)];
%     dlmwrite(strcat(snan_tot_path_out,strcat('max_',num2str(max_id),'/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_snan_box_z',num2str(z),'_data.txt'),peak,'delimiter','\t');    
% 
%     
%     %plot the zoom
%     
%     [xq,yq] = meshgrid(-angle_range/2:angle_range/n_angles_1d:angle_range/2, -angle_range/2:angle_range/n_angles_1d:angle_range/2);
%     peaks_fine_filtered_proj1d_angles=zeros(size(xq));  %matrix with the peaks to be computed
%     stn_fine_filtered_proj1d_angles=zeros(size(xq));  %matrix with the peaks to be computed
%     thetas_sq=zeros(size(xq));
%     phis_sq=zeros(size(xq));
%     [phis_q,thetas_q,~] = cart2sph(xq(:),yq(:),1);  %planar approximation
% %     thetas_q=thetas_q-pi/2                          %elevation to theta
% %     clearvars xq yq
%     [rx(1,:),rx(2,:),rx(3,:)] = sph2cart(phis_q,thetas_q,1);           %points in the sphere
%     
%     %rotation to centralize on peak
%     
%     Ry = [cos(thetas_max) 0 -sin(thetas_max); 0 1 0; sin(thetas_max) 0 cos(thetas_max)];
%     Rz = [cos(phis_max) -sin(phis_max) 0; sin(phis_max) cos(phis_max) 0; 0 0 1];
%     
%     rx=Rz*rx;
%     rx=Ry*rx;
%     
%     
%     [phis_q,thetas_q,~] = cart2sph(rx(1,:),rx(2,:),rx(3,:));
%     thetas_q=thetas_q-pi/2;                      %elevation to theta
% 
%     
%     proj1d_angles=zeros(length(bins)-1,length(phis_q));     %stores the plain 1d projections
%     
%     
%     for part_id = 1  :   part_in
%         
%         [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,~, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part_in,part_id);
%         
%         Pos=mod(Pos,nc);
%         
%         Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
%         Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
%         Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
%         
%         
%         lim_pre=(1/(lenght_factor))*nc;
%         
%         Pos(:,abs(Pos(1,:))>lim_pre|abs(Pos(2,:))>lim_pre|abs(Pos(3,:))>lim_pre)=[];
%         
%         histogr_1d_angles=zeros(length(bins)-1,length(phis_q));      %auxiliary variable that will stores all 1d projs
%         
%         
%         parfor angl_ind=1:length(phis_q)
%             
%             theta=thetas_q(angl_ind);
%             phi=phis_q(angl_ind);
%             
%             nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
%             dz=nz*Pos;
%             histogr_1d_angles(:,angl_ind)=histcounts(dz,bins);
%             
%             
%         end
%         
%         proj1d_angles=proj1d_angles+histogr_1d_angles;
%         
%     end
%     
%     %this part only computes a normalization factor, since diferent
%     %slices on the projection of the particles inside the cube have diferent
%     %areas. The normalization factor is this area.
%     
%     [v f] = createCube; v = (v-[0.5 0.5 0.5])*nc ;
%     
%     
%     parfor angl_ind=1:length(phis_q)
%         theta=thetas_q(angl_ind);
%         phi=phis_q(angl_ind);
%         nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
%         bin_sz=(length(bins)-1);
%         area=zeros(bin_sz,1);
%         for dist_1d=1:bin_sz
%             plane = createPlane(nz*((bins(1,dist_1d+1)+bins(1,dist_1d))/2), nz);
%             plane_overlap = intersectPlaneMesh(plane, v, f);
%             area(dist_1d,1)=polygonArea3d(plane_overlap);
%         end
%         areas(:,angl_ind)=area(:,1);
%     end
%     
%     %do the normalization
%     
%     proj1d_angles=proj1d_angles./areas;
%     
%     
%     %here the wavelet filter takes place (on both fileds, the particle number and the density contrast)
%     
%     filtered_proj1d_angles=zeros(length(bins)-1,length(phis_q));
%     level=floor(log2(length(bins)-1)); %level of the wavelet transformation
%     
%     parfor i=1:length(phis_q)
%         
%         proj1d=proj1d_angles(:,i);
%         [dwt_proj1d,levels] = wavedec(proj1d,level,dwbasis);
%         
%         D=zeros(length(bins)-1,1);
%         for lev_win = 1:length(level_window)
%             lvwin=level_window(lev_win);
%             D(:,1)=D(:,1)+wrcoef('d',dwt_proj1d,levels,dwbasis,lvwin);
%         end
%         filtered_proj1d_angles(:,i) = D(:,1);
%         
%     end
%     
%     %in this part the peak of each 1d projection is taked, Now we are
%     %performing these operations on the filtered 1d projections
%     
%     max_filtered_proj1d_angles=max(filtered_proj1d_angles);
%     average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
%     max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
%     
%     filtered_proj1d_index_max=zeros(1,length(phis_q));
%     
%     parfor angl=1:length(phis_q)
%         filtered_proj1d_index_max(1,angl)=find(filtered_proj1d_angles(:,angl)==max_filtered_proj1d_angles(1,angl),1);
%     end
%     filtered_proj1d_angles_snremoved=filtered_proj1d_angles;
%     
%     for angl=1:length(phis_q)
%         filtered_proj1d_angles_snremoved(filtered_proj1d_index_max(1,angl),angl)=average_filtered_proj1d_angles(angl);
%     end
%             
%     std_filtered_proj1d_angles=std(filtered_proj1d_angles_snremoved);
%     
%     stn_filtered_proj1d_angles=(max_amplitude_filtered_proj1d_angles(:))./std_filtered_proj1d_angles(:);
% 
%     
%     parfor i=1:length(phis_q)
%         
%         peaks_fine_filtered_proj1d_angles(i)=max_amplitude_filtered_proj1d_angles(i);
%         stn_fine_filtered_proj1d_angles(i)=stn_filtered_proj1d_angles(i);
% %         thetas_sq(i)=thetas_q(i);
% %         phis_sq(i)=phis_q(i);
%     end
%     %         peaks_maxima(max_id)=max(peaks_fine_filtered_proj1d_angles(:));
%     %     display(peaks_maxima)
%     
% %     fig2=figure;
%     fig2=figure('Visible', 'off');
%     %     mesh(thetas_sq,phis_sq,peaks_fine_filtered_proj1d_angles);
%     mesh(xq,yq,peaks_fine_filtered_proj1d_angles);
%     %     tot_plot_path_out_zoom=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
%     %     mkdir(tot_plot_path_out_zoom);
%     mkdir(tot_plot_path_out_zoom,strcat('max_',num2str(max_id),'/'));
%     saveas(fig2,strcat(tot_plot_path_out_zoom,'max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_peak_zoom_z',num2str(z),'_plot.png'));
% %     saveas(fig2,strcat(tot_plot_path_out_zoom,'level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_peak_zoom_z',num2str(z),'_plot'));
%     mkdir(snan_tot_path_out_zoom,strcat('max_',num2str(max_id)));
%     %     snan_zoom=[max(peaks_maxima) find(peaks_maxima==max(peaks_maxima),1)];
%     dlmwrite(strcat(snan_tot_path_out_zoom,strcat('max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_snan_box_zoom_z',num2str(z),'_peak.txt')),max(peaks_fine_filtered_proj1d_angles(:)),'delimiter','\t');
%     
%     
%     
%     
% %     fig3=figure;
%     fig3=figure('Visible', 'off');
%     %     mesh(thetas_sq,phis_sq,peaks_fine_filtered_proj1d_angles);
%     mesh(xq,yq,stn_fine_filtered_proj1d_angles);
%     %     tot_plot_path_out_zoom=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
%     %     mkdir(tot_plot_path_out_zoom);
%     %     mkdir(tot_plot_path_out_zoom,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/'));
%     saveas(fig3,strcat(tot_plot_path_out_zoom,'max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_stn_zoom_z',num2str(z),'_plot.png'));
% %     saveas(fig3,strcat(tot_plot_path_out_zoom,'level_window',mat2str(level_window(:)),'/peak_zoom/','max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_box_stn_zoom_z',num2str(z),'_plot'));
%     
% %      mkdir(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),mat2str(level_window(:)),'/peak_zoom/'));
%     %     snan_zoom=[max(peaks_maxima) find(peaks_maxima==max(peaks_maxima),1)];
%     dlmwrite(strcat(snan_tot_path_out_zoom,strcat('max_',num2str(max_id),'/_',num2str(find(str2num(char(redshift_list))==z)),'_snan_box_zoom_z',num2str(z),'_stn.txt')),max(stn_fine_filtered_proj1d_angles(:)),'delimiter','\t');
%     
%     
end

% snan_tot_path_out=strcat(root_snan_out,spec,aux_path,'snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv/','box/',dwbasis,'/anglRang_',num2str(angle_range),'/nangl_',num2str(n_angles_1d),'/');
% mkdir(snan_tot_path_out);
% mkdir(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'));
% snan_zoom=[max(peaks_maxima) find(peaks_maxima==max(peaks_maxima),1)];
% dlmwrite(strcat(snan_tot_path_out,strcat('level_window',mat2str(level_window(:)),'/peak_zoom/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_snan_box_zoom_z',num2str(z),'_data.txt'),snan_zoom,'delimiter','\t');    

cd('../wake_detection/box_statistics');

% delete(gcp('nocreate'))


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

function [f_index_max] = find_peaks(thetas, phis,f,NSIDE,number_of_maxima)

% number_of_maxima=4;

% peak=max(f);
% ave=mean(f);
% f_symmetric=f-ave;
% signal=max(f_symmetric);

f_size=length(f);
f(f_size/2:end)=[];

[f_sort,f_index_max]=sort(f,'descend');

% f_sort=flip(unique(sort(f,'descend')));


% f_symmetric_sort=fliplr(f_symmetric_sort);

% std_f=std(f);
% stn_f=signal/std_f;
% f_index_max=find(f_symmetric==signal);
% f_index_max=zeros(1,number_of_maxima);
% for pk=1:length(f_index_max)
%     f_index_max(1,pk)=find(f==f_sort(pk),1);
% end

phis_max=phis(f_index_max);
thetas_max=thetas(f_index_max);
% display(phis_max(1:10));
% display(thetas_max(1:10));

angle_cutoff=2*sqrt(4*pi/(12*(NSIDE/2)*(NSIDE/2)));

for f_verify =1:number_of_maxima
    
    f_index_max_aux=length(f_index_max);
    
%     display(f_index_max_aux)
        	v1=zeros(1,3);
            [v1(1,1),v1(1,2),v1(1,3)] = sph2cart(phis_max(f_verify),+pi/2-thetas_max(f_verify),1);
            v2=-v1;
    
    parfor f_compare=f_verify+1:f_index_max_aux

            w1=zeros(1,3);
            [w1(1,1),w1(1,2),w1(1,3)] = sph2cart(phis_max(f_compare),+pi/2-thetas_max(f_compare),1);
            w2=-w1;
            
            a=zeros(1,4);
            a(1,1)=acos(dot(v1, w1));
            a(1,2)=acos(dot(v1, w2));
            a(1,3)=acos(dot(v2, w1));
            a(1,4)=acos(dot(v2, w2));
            am=min(a);
            
            if (am<angle_cutoff)
                
                f_index_max(f_compare)=0;
                
            end
    end
    
    f_index_max(f_index_max==0) = [];
    
end

end



function [x, y] = ssht_mollweide(thetas, phis,thetas_shift, phis_shift)
% ssht_mollweide - Compute Mollweide projection
%
% Compute Mollweide projection of spherical coordinates.
%
% Usage is given by
%
%   [x,y] = ssht_mollweide(thetas, phis)
%
% where thetas and phis are spherical coordinates and x and y are the
% projected Mollweide coordinates.
%
% Author: Jason McEwen (www.jasonmcewen.org)

MAX_ITERATIONS = 1e5;
TOL = 1e-10;


% Convert theta to longitude.
thetas = pi/2 - thetas;
phis = phis - pi;


% recenter
[rx(1,:),rx(2,:),rx(3,:)] = sph2cart(phis,thetas,1);

% theta=pi/2+thetas_shift;
% phi=-phis_shift;

theta=pi/2;
phi=pi/4;

Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

rx=Rz*rx;
rx=Ry*rx;

[phis,thetas,~] = cart2sph(rx(1,:),rx(2,:),rx(3,:));


t = thetas;
for it = 1:MAX_ITERATIONS

   dt = (t + sin(t) - pi.*sin(thetas)) ./ (1 + cos(t));
   t = t - dt;
   
   if(max(abs(dt)) < TOL)
      break;
   end
   
end
t = t/2;
x = 2 .* sqrt(2) ./ pi .* phis .* cos(t);
y = sqrt(2) .* sin(t);
end