function [  ] = proj1d_halos_plot_all( root,root_data_out,root_out,spec,aux_path,aux_path_out,lenght_factor,resol_factor,pivot,rot_angle,lim )
% reads data of the 1d projections aconding to the input specifications and
% plot the result for all files in the folder

% (example) proj1d_halos_plot_all('/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','',1,1,[0,0,0],[0,0],'minmax');
% (example) proj1d_halos_plot_all('/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',1,1,[0,0,0],[0,0],[-1 2]);

cd('../preprocessing');

path_in=strcat(root,spec,aux_path);
[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path);

cd('../1dproj');

for rds = 1 : length(redshift_list)
    
    filename = dir(strcat(path_in,char(redshift_list(rds)),'halo0','.dat'));
    filename=filename.name;
    proj1d_halos_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim);
    
end

end