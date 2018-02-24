function [ pks,locs,out_proj1d_angles ,out_filtered_proj1d_angles,f_max,thetas_max,phis_max] = box_statistics_dm_analysis( root,root_data_out,root_plot_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,lim,data_stream,info,analysis ,cutoff,projection_analysed,signal_analysed)

%(example)  [  ] = box_statistics_dm_analysis('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/data_test/','/home/asus/Dropbox/extras/storage/graham/small_res/test_plot_box/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','','0.000xv0.dat',2,1,[0,0,0],4,1,4,'minmax',[1,3],[0,1,2,3],1,10,[1,2,3,4],[1,2,3]);

%(example)  [  ] = box_statistics_dm_analysis('/home/asus/Dropbox/extras/storage/guillimin/', '/home/asus/Dropbox/extras/storage/guillimin/box_stat/','/home/asus/Dropbox/extras/storage/guillimin/box_stat/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','','10.000xv0.dat',2,1,[0,0,0],64,16,4,'minmax',[1],[0,1,2,3],1,0.8,[1,2,3,4],[1,2,3]);

% 
% path_total_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj/');
% path_analysis_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj_analysis/');

cd('../../preprocessing');

[~,redshift_list,~,size_box,nc,np,zi,~,~,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/wake_detection/box_statistics');

if ismember(0,data_stream)
    if (ismember(3,projection_analysed)||ismember(4,projection_analysed)||ismember(5,projection_analysed))
        [~,~,out_proj1d_angles,out_dc_proj1d_angles,out_filtered_proj1d_angles,out_filtered_dc_proj1d_angles,out_dc_filtered_proj1d_angles] =box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,0,cutoff);
    else
        [~,~,out_proj1d_angles,out_dc_proj1d_angles,~,~,~] =box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,0,cutoff);
    end
else
    path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    
    if ismember(2,data_stream)
        if ismember(3,data_stream)
            box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,2,cutoff);
        end
        out_proj1d_angles=dlmread(strcat(path_data,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
        out_dc_proj1d_angles=dlmread(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
    end
    if ismember(1,data_stream)
        if ismember(3,data_stream)
            box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,1,cutoff);
        end
%         display(strcat(path_data,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
        fileID = fopen(strcat(path_data,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
        out_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
        out_dc_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
        fclose(fileID);
    end
    
    if (ismember(3,projection_analysed)||ismember(4,projection_analysed)||ismember(5,projection_analysed))
        if ismember(2,data_stream)
            out_filtered_proj1d_angles=dlmread(strcat(path_data,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
            out_dc_filtered_proj1d_angles=dlmread(strcat(path_data,'dc/','cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
            out_filtered_dc_proj1d_angles=dlmread(strcat(path_data,'cutoff_',num2str(cutoff),'MpcCut/','dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
        end
        if ismember(1,data_stream) 
            fileID = fopen(strcat(path_data,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
            out_filtered_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
            fclose(fileID);
            
%             fileID = fopen(strcat(path_data,'cutoff_',num2str(cutoff),'MpcCut/','dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
%             out_dc_filtered_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
%             fclose(fileID);
            
            fileID = fopen(strcat(path_data,'dc/','cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
            out_filtered_dc_proj1d_angles=fread(fileID,[3,12*NSIDE^2],'float32','l');
            fclose(fileID);
        end
    end
end

addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab','/home/asus');
   
if ismember(1,projection_analysed)
    if ismember(1,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_proj1d_angles(1,:)),'Peak of the mass amplitude');
    end
    
    if ismember(2,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_proj1d_angles(2,:)),{'Standard deviation';'of the mass amplitude'})
    end
    
    if ismember(3,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_proj1d_angles(3,:)),{'Signal to noise';'of the mass amplitude'})
    end
    
end

if ismember(2,projection_analysed)
    if ismember(1,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_dc_proj1d_angles(1,:)),'Peak of the dc amplitude');
    end
    
    if ismember(2,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_dc_proj1d_angles(2,:)),{'Standard deviation';'of the dc amplitude'})
    end
    
    if ismember(3,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_dc_proj1d_angles(3,:)),{'Signal to noise';'of the dc amplitude'})
    end
    
end

if ismember(3,projection_analysed)
    if ismember(1,signal_analysed)
       [f_index_max,f_max,thetas_max,phis_max]= plot_molweide(transpose(out_filtered_proj1d_angles(1,:)),'Peak of the filtered mass amplitude');
    end
    
    if ismember(2,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_proj1d_angles(2,:)),{'Standard deviation';'of the filtered mass amplitude'})
    end
    
    if ismember(3,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_proj1d_angles(3,:)),{'Signal to noise';'of the filtered mass amplitude'})
    end
    
    %multiplication of the filtered amplitude and standard deviation
    
    if ismember(4,signal_analysed)
%         plot_molweide(transpose(out_filtered_proj1d_angles(1,:)).*transpose(out_filtered_proj1d_angles(2,:)),{'Peak of the filtered mass amplitude';'times Standard deviation';'of the filtered mass amplitude'})
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_proj1d_angles(1,:)).*transpose(out_filtered_proj1d_angles(2,:)),{'Original Map'})    
    end    
end


            
if ismember(4,projection_analysed)
    if ismember(1,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_dc_proj1d_angles(1,:)),'Peak of the filtered dc amplitude');
    end
    
    if ismember(2,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_dc_proj1d_angles(2,:)),{'Standard deviation';'of the filtered dc amplitude'})
    end
    
    if ismember(3,signal_analysed)
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_dc_proj1d_angles(3,:)),{'Signal to noise';'of the filtered dc amplitude'})
    end
    
    if ismember(4,signal_analysed)
%         plot_molweide(transpose(out_filtered_proj1d_angles(1,:)).*transpose(out_filtered_proj1d_angles(2,:)),{'Peak of the filtered mass amplitude';'times Standard deviation';'of the filtered mass amplitude'})
        [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(transpose(out_filtered_dc_proj1d_angles(1,:)).*transpose(out_filtered_dc_proj1d_angles(2,:)),{'Original Map'})    
    end
    
end

%multiplication of the filtered amplitude and standard deviation

% if ismember(5,projection_analysed)
%     if ismember(1,signal_analysed)
%         plot_molweide(transpose(out_dc_filtered_proj1d_angles(1,:)),'Peak of the dc filtered mass amplitude');
%     end
%     
%     if ismember(2,signal_analysed)
%         plot_molweide(transpose(out_dc_filtered_proj1d_angles(2,:)),{'Standard deviation';'of the dc filtered mass amplitude'})
%     end
%     
%     if ismember(3,signal_analysed)
%         plot_molweide(transpose(out_dc_filtered_proj1d_angles(3,:)),{'Signal to noise';'of the dc filtered mass amplitude'})
%     end
%     
% end

% if cutoff~=0
%     
%     plot_molweide(transpose(out_filtered_proj1d_angles(1,:)),'filtered data')
% 
% end

%save the wake possible position (angle and location)

n_angle_peaks=2;
n_1d_peaks=2
path_data_all=strcat(strcat(root_data_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');

% fileID = fopen(strcat(path_data_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
% proj1d_angles=fread(fileID,[np*resol_factor/lenght_factor ,12*NSIDE^2],'float32','l');
% fclose(fileID);
% 
% figure;
% plot(proj1d_angles(:,f_index_max(1:2)));
 
% fileID = fopen(strcat(path_data_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
% dc_proj1d_angles=fread(fileID,[np*resol_factor/lenght_factor ,12*NSIDE^2],'float32','l');
% fclose(fileID);
% 
% figure;
% plot(dc_proj1d_angles(:,f_index_max(1:2)));
% 
% 
% fileID = fopen(strcat(path_data_all,'/cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
% filtered_proj1d_angles=fread(fileID,[np*resol_factor/lenght_factor ,12*NSIDE^2],'float32','l');
% fclose(fileID);
% 
% figure;
% plot(filtered_proj1d_angles(:,f_index_max(2)));


fileID = fopen(strcat(path_data_all,'/cutoff_',num2str(cutoff),'MpcCut/dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
filtered_dc_proj1d_angles=fread(fileID,[np*resol_factor/lenght_factor ,12*NSIDE^2],'float32','l');
fclose(fileID);

figure;
plot(filtered_dc_proj1d_angles(:,f_index_max(1:n_angle_peaks)));

peaks_tot=[];
for pk_angl_indx=1:n_angle_peaks
    [pks,locs] =findpeaks(filtered_dc_proj1d_angles(:,f_index_max(pk_angl_indx)),'NPeaks',n_1d_peaks,'SortStr','descend');
%     peaks=[pks,((locs-np*resol_factor/(2*lenght_factor))*2/resol_factor),thetas_max(1:n_angle_peaks),phis_max(1:n_angle_peaks)];
    for pk_1d_indx=1:n_1d_peaks
        peaks(1,1)=pks(pk_1d_indx);
        peaks(1,2)=(locs(pk_1d_indx)-np*resol_factor/(2*lenght_factor))*2/resol_factor;
        peaks(1,3)=thetas_max(pk_angl_indx);
        peaks(1,4)=phis_max(pk_angl_indx);
        peaks_tot=[peaks_tot; peaks];
    end
end

path_data_out_curv=strcat(strcat(root_data_out(1,1:end-1),'_curv/',spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
mkdir(strcat(strcat(root_data_out(1,1:end-1),'_curv/',spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/'));
dlmwrite(strcat(path_data_out_curv,'_',num2str(find(str2num(char(redshift_list))==z)),'_peaks_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),peaks_tot,'delimiter','\t');

end


function [f_index_max,f_max,thetas_max,phis_max]=plot_molweide(f,string)

%     s2let_hpx_plot_mollweide(maxi);
    
    sz = size(f);
nsideguessed = sqrt(max(sz)/12);
    L = 2*nsideguessed;
    B=2;
    J_min=0;
    
    [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f,'B',B,'L',L,'J_min',J_min);
    
    % Plot
J = s2let_jmax(L, B);
zoomfactor = 1.2;
% zoomfactor = 20;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = 2  ;
nx = ceil((2+J-J_min)/2 );
figure('Position',[100 100 1300 1000])

subplot(nx, ny, 1);
s2let_hpx_plot_mollweide(f);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title(string)

% subplot(nx, ny, 2);
% s2let_hpx_plot_mollweide(f_scal);
% campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
% title('Scaling fct')

% display(J_min);
% display(J);
% display(ns);

for j = J_min:J
   subplot(nx, ny, j-J_min+2);
   s2let_hpx_plot_mollweide(f_wav{j-J_min+1});
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   title(['Wavelet scale : ',int2str(j)-J_min+1])
      
end 

% for j = J_min:J
%    
%          s2let_hpx_plot_mollweide_info(f_wav{j-J_min+1},1);
% 
% end


% s2let_hpx_plot_mollweide_info(f_wav{10},1);
s2let_hpx_plot_mollweide_info(f,1);
% s2let_hpx_plot_mollweide(f);
% 
% figure;
% histogram(f);

% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1},1);

% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1}+f_wav{J-J_min-2}+f_wav{J-J_min-3},1); %best for gmu2t10m7
% 
% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1}+f_wav{J-J_min-2},1); %best for gmu4t10m7
% 
% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1},1); %best for Gmu8t10m7

% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min},1);

s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1},1);


% s2let_hpx_plot_mollweide_info(f,1); %best for za gmu10m7,nside32, (3,1)


% s2let_hpx_plot_mollweide_info(f_wav{6},1); %no signal for gmu2t10m7


% s2let_hpx_plot_mollweide_info(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1}+f_wav{J-J_min-2}+f_wav{J-J_min-3}+f_wav{J-J_min-4}+f_wav{J-J_min-5},1); %best for gmu2t10m7
% [f_index_max,f_max,thetas_max,phis_max] = s2let_hpx_plot_mollweide_findpeaks(f_wav{J-J_min+1}+f_wav{J-J_min}+f_wav{J-J_min-1}+f_wav{J-J_min-2}+f_wav{J-J_min-3}+f_wav{J-J_min-4},1); 

% [f_index_max,f_max,thetas_max,phis_max] = s2let_hpx_plot_mollweide_findpeaks(f,1); 

% s2let_hpx_plot_mollweide_info(f_wav{5},1); %for the cubic scheme nside64, 3.4
% s2let_hpx_plot_mollweide_info(f_wav{5}+f_wav{4},1); %same as above
% s2let_hpx_plot_mollweide_info(f_wav{5}+f_wav{4},1); %same as above, but
% for 4,1
% 
% s2let_hpx_plot_mollweide_info(f_wav{6},1);
% 
% figure;
% 
% s2let_hpx_plot_mollweide_info(f_wav{10}+f_wav{9}+f_wav{8}+f_wav{7}+f_wav{6}+f_wav{5});

% display(thetas_max);
% display(phis_max);

end