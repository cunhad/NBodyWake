function [  ] = proj2d_halos_plot_all( root,root_data_out,root_out,spec,aux_path,aux_path_out,lenght_factor,resol_factor,pivot,rot_angle,lim )
% reads data of the 2d projections aconding to the input specifications and
% plot the result for all files in the folder

% (example) proj2d_halos_plot_all('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',1,1,[0,0,0],[0,0],'minmax');


cd('../preprocessing');

path_in=strcat(root,spec,aux_path);
[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path);

cd('../2dproj');

for rds = 1 : length(redshift_list)
    
    filename = dir(strcat(path_in,char(redshift_list(rds)),'halo0','.dat'));
    filename=filename.name;
    proj2d_halos_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim);
    
end

end