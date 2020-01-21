function [ anali,signal,equator_phi ] = slices_curv_2a3d_locangrealCanny_depth_all( root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth)

%(example)  [ anali] = slices_curv_2a3d_all('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','/home/asus/Dropbox/extras/storage/graham/small_res/anali_hpx/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,2,2 );
% 
% for i=3001:3010
%     spec='4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m';
%     aux_path=strcat('/sample',num2str(i),'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% for i=4001:4010
%     spec='16Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
%     aux_path=strcat('/sample',num2str(i),'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% for i=4001:4010
%     spec='16Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m';
%     aux_path=strcat('/sample',num2str(i),'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% 
% for i=4001:4010
%     spec='8Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
%     aux_path=strcat('/sample',num2str(i),'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% 
% for i=3001:3010
%     spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
%     aux_path=strcat('/sample',num2str(i),'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end
% 

% for i=4001:4010
%     spec='16Mpc_2048c_1024p_zi63_nowakem';    
%     aux_path=strcat('/sample',num2str(i),'/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% for i=4001:4010
%     spec='8Mpc_2048c_1024p_zi63_nowakem';    
%     aux_path=strcat('/sample',num2str(i),'/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end

% for i=3001:3010
%     spec='4Mpc_2048c_1024p_zi63_nowakem';    
%     aux_path=strcat('/sample',num2str(i),'/');
%     slices_curv_2a3d_locangrealCanny_depth_all(root,root_anali_2d_in,root_2d_anali_hpx,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,slice,NSIDE,analy_type,sum_depth);
% end


% root='/home/asus/Dropbox/extras/storage/graham/small_res/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/small_res/anali3/';
% root_2d_anali_hpx='/home/asus/Dropbox/extras/storage/graham/small_res/anali3_hpx/';
% spec='64Mpc_256c_128p_zi63_nowakem';
% aux_path='/sample2001/';
% aux_path_out='';
% filename='10.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% slice=32;
% NSIDE=8;
% analy_type=3;
% sum_depth=4;

% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclarhpxNSIDE8l2lr1na256_to_3dparcurv-2dclarhpxNSIDE8l1lr1_anali/';
% root_2d_anali_hpx='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclarhpxNSIDE8l2lr1na256_to_3dparcurv-2dclarhpxNSIDE8l1lr1_anali3_dp4_hpx/';
% spec='4Mpc_2048c_1024p_zi63_nowakem';
% aux_path='/sample3001/';
% % spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3001/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% slice=32;
% NSIDE=8;
% analy_type=3;
% sum_depth=4;
% 
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaraCannyl1lr1_anali/';
% root_2d_anali_hpx='/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaraCannyl1lr1_anali_hpx_Sa4t1_dp4/';
% spec='4Mpc_2048c_1024p_zi63_nowakem';
% % aux_path='/sample3001/';
% aux_path='/sample4001/';
% % spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3001/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=0.5;
% slice=32;
% NSIDE=4;
% analy_type=1;
% sum_depth=4;


angles_hpx(1,:) = dlmread(strcat('../../../python/angles',num2str(NSIDE),'_t.cvs'));
angles_hpx(2,:) = dlmread(strcat('../../../python/angles',num2str(NSIDE),'_p.cvs'));

cd('../../preprocessing');

[~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,1024,1);


z_glob=str2num(filename(1:end-7))
z=filename(1:end-7);

N_angles=12*NSIDE*NSIDE/2;
N_angles_t=12*NSIDE*NSIDE;



for angle_id=1:N_angles
    
    path1=strcat(root_anali_2d_in,spec,aux_path,'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(angle_id),'/')
    
    path2=dir(char(strcat(path1,"*pv*")));
    path2={path2.name}
    path2=path2(1);
    
    path3='/2dproj/dm/';
    
    path_in=string(strcat(strcat(root_anali_2d_in,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(angle_id),'/',path2,'/2dproj/dm/'))
    
    aux=dlmread(strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3_depth',num2str(sum_depth),'.txt'));
    
    anali(angle_id,:,:,:)=reshape(aux,slice/sum_depth,4,5);
    
end


% signal=reshape(permute(anali(:,:,2,1),[1,3,2,4,5]),[1,numel(anali(:,:,2,1))]);
% signal = reshape(signal,[N_angles,2])';
% signal = reshape(signal,[slice,N_angles]);
signal=anali(:,:,4,analy_type);
signal=sum(signal,2);
signal=signal';
% signal = max(signal);

signal_det=anali(:,:,4,analy_type);

for angl=N_angles+1:N_angles_t
    
    theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));    
    theta_inverse_indices=N_angles_t-theta_indices+1;    
    phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
    signal(:,angl)=flip(signal(:,phi_indice_inverse));
    
end


cd('../wake_detection/slices_curv_2d/');





%plot/data out

tot_anali_path_out=char(strcat(root_2d_anali_hpx,spec,aux_path,'anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/slices_curv_2a3d/hpx/'))
mkdir(tot_anali_path_out);
        
tot_anali_data_path_out=char(strcat(root_2d_anali_hpx,spec,aux_path,'data_anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/slices_curv_2a3d/hpx/'))
mkdir(tot_anali_data_path_out);


%molweide projection
    
    fig=figure('Visible', 'off');
% fig=figure
 [thetas, phis] = s2let_hpx_sampling_ring(NSIDE);
 
 [f_index_max] = find_peaks(thetas, phis,signal(:),NSIDE,1);
 peak=signal(f_index_max(1));
 
        display(peak)
        phis_max=phis(f_index_max(1));
        thetas_max=thetas(f_index_max(1));
        display(phis_max)
        display(thetas_max)
 
 
 [x, y] = ssht_mollweide(thetas, phis,0,0);
 
 
    gridDelaunay = delaunay(x,y);
    h = trisurf(gridDelaunay,x,y,signal(:)*0.0,signal(:));
    set(h, 'LineStyle', 'none')
    axis equal
    axis off
    campos([0 0 1])
    camup([0 1 0])
    
    mkdir(tot_anali_path_out,strcat('/molvp/'));
    saveas(fig,strcat(tot_anali_path_out,strcat('/molvp'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_plot.png'));
   
    mkdir(tot_anali_data_path_out,strcat('/molvd/'));
    dlmwrite(strcat(tot_anali_data_path_out,strcat('/molvd'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_data.txt'),[thetas,phis,signal'],'delimiter','\t')
    dlmwrite(strcat(tot_anali_data_path_out,strcat('/molvd'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_data_det.txt'),signal_det,'delimiter','\t')
close(fig)
    
%molweide projection details


    
    fig=figure('Visible', 'off');
% fig=figure
 
    set(gcf, 'Position', [0 0 1200 600]);
        ax1 = axes('Position',[0.01 0.13 0.7 0.7]);
    gridDelaunay = delaunay(x,y);
    h = trisurf(gridDelaunay,x,y,signal(:)*0.0,signal(:));
    set(h, 'LineStyle', 'none')
    hold on
    axis equal
    axis off
    campos([0 0 1])
    camup([0 1 0])
    [x_max, y_max] = ssht_mollweide(thetas_max, phis_max,0,0);
    colorbar('southoutside');
        scatter(x_max, y_max,255,'r');
        descr = {strcat('max = ',num2str(peak));
            strcat('mean = ',num2str(mean(signal(:))));
            strcat('std = ',num2str(std(signal(:))));
            strcat('satn = ',num2str((peak-mean(signal(:)))/std(signal(:))));
            strcat('kurt = ',num2str(kurtosis(signal(:))))};
            
        ax1 = axes('Position',[0 0 1 1],'Visible','off');        
        txt=text(0.75,0.5,descr);
        set(txt,'Parent',ax1,'interpreter', 'latex');
        hold off;
       
    
    mkdir(tot_anali_path_out,strcat('/molvp_det/'));
    saveas(fig,strcat(tot_anali_path_out,strcat('/molvp_det'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_plot.png'));
   
    dlmwrite(strcat(tot_anali_data_path_out,strcat('/molvd'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_data_summary.txt'),[thetas_max, phis_max,peak,std(signal(:)),peak/std(signal(:))],'delimiter','\t')
    
    close(fig)
    %molweide projection histogram


    
    fig=figure('Visible', 'off');
% fig=figure
 
    
    histogram(signal(:));
   
    hold on
    
        hold off;
       
    
    mkdir(tot_anali_path_out,strcat('/hist/'));
    saveas(fig,strcat(tot_anali_path_out,strcat('/hist'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_hist_hpx_slice_cuvr_2a3d_z',z,'_plot.png'));
   

    
    %plot equator
    
    equator_idx=find(abs(thetas-pi/2)==0);
    equator_signal=signal(equator_idx);
    equator_phi=phis(equator_idx);
    
    max_indx=find(equator_signal==max(equator_signal),1);
    close(fig)
    
        fig=figure('Visible', 'off');
% figure;
plot(equator_phi,equator_signal);
xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$max$', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

    mkdir(tot_anali_path_out,strcat('/equator/'));
    saveas(fig,strcat(tot_anali_path_out,strcat('/equator'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_equ_slice_cuvr_2d_z',z,'_plot.png'));
   
    title({strcat('max = ',num2str(max(equator_signal(:)))),strcat('std = ',num2str(std(equator_signal(:)))),strcat('max','/','std = ',num2str(max(equator_signal(:))/std(equator_signal(:))))});
    
    mkdir(tot_anali_path_out,strcat('/equator_det/'));
    saveas(fig,strcat(tot_anali_path_out,strcat('/equator_det'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_equ_d_slice_cuvr_2a3d_z',z,'_plot.png'));
   hold off
   
   mkdir(tot_anali_data_path_out,strcat('/equator_data/'));
    dlmwrite(strcat(tot_anali_data_path_out,strcat('/equator_data'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_equdata_slice_cuvr_2a3d_z',z,'_data.txt'),[equator_phi,equator_signal'],'delimiter','\t')
    dlmwrite(strcat(tot_anali_data_path_out,strcat('/equator_data'),'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_equdata_slice_cuvr_2a3d_z',z,'_data_summary.txt'),[equator_phi(max_indx),max(equator_signal(:)),std(equator_signal(:)),max(equator_signal(:))/std(equator_signal(:))],'delimiter','\t')
close(fig)
    
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

theta=pi/2+thetas_shift;
phi=-phis_shift;

% theta=pi/2;
% phi=pi/4;

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
    
%     parfor f_compare=f_verify+1:f_index_max_aux
    for f_compare=f_verify+1:f_index_max_aux

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
