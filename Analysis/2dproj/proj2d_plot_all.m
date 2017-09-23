function [  ] = proj2d_plot_all( path,spec,aux_path,lenght_factor,resol_factor,pivot,rot_angle,lim )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cd('../preprocessing');

path_in=strcat(path,spec,aux_path);
[ nodes_list redshift_list ] = preprocessing_many_nodes(path,spec,aux_path);

cd('../2dproj');

for rds = 1 : length(redshift_list)
    
    filename = dir(strcat(path_in,char(redshift_list(rds)),'xv0','.dat'));
    filename=filename.name;
    proj2d_plot( path,spec,aux_path,filename,lenght_factor,resol_factor,pivot,rot_angle,lim);
    
end

end