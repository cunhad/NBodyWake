function [  ] = wavelets_1d_dm_analysis_all( root,root_data_out,root_plot_out,root_snan_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,aux_path_snan_out,lenght_factor,resol_factor,pivot,rot_angle,lim1d,lim_cwt,cutoff,data_stream,info,analysis,z_id_range)
% reads (and/or generate) data of the filtered 1d projections aconding to the input specifications and
% plot the result for all files in the folder

% (example) wavelets_1d_dm_analysis_all('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','/home/asus/Dropbox/extras/storage/graham/small_res/snan/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','','',1,1,[0,0,0],[0,0],'minmax','minmax',10,[1,3],[0,2,3],1,'all');
% (example) wavelets_1d_dm_analysis_all('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/data/','/home/asus/Dropbox/extras/storage/guillimin/plot/','/home/asus/Dropbox/extras/storage/guillimin/snan/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','','',1,1,[0,0,0],[0,0],'minmax','minmax',0.8,[1,3],[0,2,3],1,'all');

% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% plot will be stored in  root_plot_out+spec+aux_path+aux_path_plot_out

% if specified, data will be stored in  root_data_out+spec+aux_path+aux_path_data_out

% if specified, signal to noise analysis will be stored in  root_snan_out+spec+aux_path+aux_path_snan_out

% filename is the output file from the nbody simulation

% lenght_factor = the analysis cube will have a lateral size given by the
% lateral size of the simulation cube divided by this number

% resol_factor= the bin will hte the particle bin size divided by this
%number

% pivot = a 3d array containing the translation wrt to the center of the
% cube (in grid cell units)

%rot_angle = 2d aray containing the theta and phy spherical angles pointing
%to the direction where the new z axis will be rotated

%cutoff is the lenght scale wich the fluctuations will be removed if above
%that. In Mpc.

%lim1d= limits on the y axis of the plot, in array format. If set to 'minmax'
%will display between the min and max values

%lim_cwt= limits the display of the absolute values of the cont wave transf. If set to 'minmax'
%will display between the min and max values


% data_stream=[1,2,3]
% if data_stream = 0, no data output generated and readed
% if data_stream = 1, reads data binaries 
% if data_stream = 2, reads data text 
% if data_stream = 3, generates the data output in binary or text if 2 or 3 options are given, respectively


% info=[0,1,2,3]
% if info=0, histogram of each plot is generated
% if info=1, minimal plots are generated
% if info=2 complete plots are generated

% analysis=1 -> create a textfile with signal to noise data (peak, std, peak/std)

% z_id_range = array with the redshift id of the requested plots and
% analysis, which starts
% with the highest one equals to 1 and decreasing by unit as the redshift
% is decreased for the id convention. If set to "all" will do for every
% redshift


cd('../../preprocessing');

path_in=strcat(root,spec,aux_path);
[~,redshift_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

cd('../wake_detection/lets');

if ischar(z_id_range)
    z_id_range=[1 : length(redshift_list)];
end

for rds = z_id_range
    
    filename = strcat(char(redshift_list(rds)),'xv0','.dat');
    wavelets_1d_dm_analysis( root,root_data_out,root_plot_out,root_snan_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,aux_path_snan_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim1d,lim_cwt,cutoff,data_stream,info,analysis);
    

end