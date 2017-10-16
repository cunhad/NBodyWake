function [ periods dc_cwt  ] = wavelets_1d_dm_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,cutoff)
 
%(example) [ proj1d ] = wavelets_1d_dm_data_out('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','','0.000xv0.dat',1,1,[0,0,0],[0,0],10);

%p = parpool(num_cores);
%tic;
path_in=strcat(root,spec,aux_path);

cd('../../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,filename);

cd('../1dproj');

%mkdir(root_out);
%mkdir(root_out,strcat(spec,aux_path));
%mkdir(strcat(root_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/'));
%path_out=strcat(root_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');


path_orig_data=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');

path_out=strcat(path_orig_data,'wavelet/dc/');
mkdir(path_orig_data,strcat('wavelet/dc/'));
mkdir(path_orig_data,strcat('wavelet/dc/scale/'));

proj1d_dm_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);
proj1d=dlmread(char(strcat(path_orig_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.txt')));

average=mean2(proj1d);
proj1d=(proj1d-average)/average;

%     %wavelet analysis
%     
% 
%     %hold on;
%     fig=figure('Visible', 'off');
%     

[dc_cwt,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);


dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_data.txt'),dc_cwt,'delimiter','\t');
dlmwrite(strcat(path_out,'scale/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_scale.txt'),seconds(periods),'delimiter','\t');


%proj1d with the cutoff

i_dc_cwt = icwt(dc_cwt,periods,[periods(1) seconds(cutoff)],'waveletparameters',[3 3.01]);
     
mkdir(path_orig_data,strcat('wavelet/dc/filter_1dproj/'));
dlmwrite(strcat(path_orig_data,strcat('wavelet/dc/filter_1dproj/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.txt'),i_dc_cwt,'delimiter','\t');

%filtered wavelet coeficients

[filtered_dc_cwt,periods] = cwt(i_dc_cwt,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);

mkdir(path_orig_data,strcat('wavelet/dc/filtered_cwt/'));
dlmwrite(strcat(path_orig_data,strcat('wavelet/dc/filtered_cwt/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.txt'),filtered_dc_cwt,'delimiter','\t');


     
% hp = pcolor( 0:size_box/length(proj1d):size_box-size_box/length(proj1d),seconds(periods),abs(mcwt)); hp.EdgeColor = 'none';
% colorbar;
% set(gca,'YScale','log');




% %     hp = pcolor( 0:size_box/length(count_sum):size_box-size_box/length(count_sum),seconds(periods),abs(mcwt)); 
% %     hp.EdgeColor = 'none';
% 
 %imagesc(0:size_box/length(proj1d):size_box-size_box/length(proj1d),seconds(periods),abs(mcwt));
 
 

 
% 
% set(gca,'YScale','log');
% xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20); ylabel('Scale parameter (Mpc)', 'interpreter', 'latex', 'fontsize', 20);
% title(strcat('Continuous wavelet transformation (Morse)'),'interpreter', 'latex', 'fontsize', 20);
% mkdir(path_analysis_out,strcat('wavelets/mtwt/'));
% path_file_out=strcat(path_analysis_out,'wavelets/mtwt/','_',num2str(rds),'_mcwt_z',num2str(z),'.png');
% saveas(fig,path_file_out);
% hold off;
% 
% fig=figure('Visible', 'off');
% clims = [0 1];
% imagesc(0:size_box/length(count_sum):size_box-size_box/length(count_sum),seconds(periods),abs(mcwt),clims);
% colorbar;
% set(gca,'YScale','log');
% xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20); ylabel('Scale parameter (Mpc)', 'interpreter', 'latex', 'fontsize', 20);
% title(strcat('Continuous wavelet transformation (Morse)'),'interpreter', 'latex', 'fontsize', 20);
% mkdir(path_analysis_out,strcat('wavelets/mtwt_lims/'));
% path_file_out=strcat(path_analysis_out,'wavelets/mtwt_lims/','_',num2str(rds),'_mcwt_lims_z',num2str(z),'.png');
% saveas(fig,path_file_out);
% hold off;
%     
%     
% end
% 

cd('../wake_detection/lets');

% 
% %toc;
% 
% %delete(gcp('nocreate'))

end

