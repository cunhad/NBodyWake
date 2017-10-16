function [  ] = wavelets_1d_dm_plot_all( root,root_data_out,root_out,spec,aux_path,aux_path_out,lenght_factor,resol_factor,pivot,rot_angle,lim1d,lim_cwt,cutoff  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (example) wavelets_1d_dm_plot_all('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',1,1,[0,0,0],[0,0],'minmax','minmax',10);

cd('../../preprocessing');

path_in=strcat(root,spec,aux_path);
[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path);

cd('../wake_detection/lets');

for rds = 1 : length(redshift_list)
    
    filename = dir(strcat(path_in,char(redshift_list(rds)),'xv0','.dat'));
    filename=filename.name;
    wavelets_1d_dm_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim1d,lim_cwt,cutoff );
    

end

