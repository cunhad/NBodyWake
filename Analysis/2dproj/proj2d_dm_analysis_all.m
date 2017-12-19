function [  ] = proj2d_dm_analysis_all( root,root_data_out,root_plot_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,lenght_factor,resol_factor,pivot,rot_angle,lim,data_stream,info,z_id_range)
% reads data of the 2d projections aconding to the input specifications and
% plot the result for all files in the folder

% (example) proj2d_dm_analysis_all('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','',1,1,[0,0,0],[0,0],'minmax',[1,3],[0,1,2,3],'all');

cd('../preprocessing');

[~,redshift_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

if ischar(z_id_range)
    z_id_range=[1 : length(redshift_list)];
end

cd('../2dproj');

for rds = z_id_range
    
    filename = strcat(char(redshift_list(rds)),'xv0','.dat');
    proj2d_dm_analysis( root,root_data_out,root_plot_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim,data_stream,info);
    
end

end