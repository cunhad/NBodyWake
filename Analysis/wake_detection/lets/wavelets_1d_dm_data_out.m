function [ dc_cwt,periods,i_dc_cwt,filtered_dc_cwt] = wavelets_1d_dm_data_out( root,root_data,spec,aux_path,aux_path_data,filename,lenght_factor,resol_factor,pivot,rot_angle,cutoff,data_stream)
 
%Computes the filetered 1d projections aconding to the input specifications and stores (and/or returns) the resulting data


%(example) [ dc_cwt periods i_dc_cwt filtered_dc_cwt] = wavelets_1d_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',1,1,[0,0,0],[0,0],10,[1]);
%(example) wavelets_1d_dm_data_out('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','15.000xv0.dat',1,1,[0,0,0],[0,0],0.8);

% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% if specified, data will be stored in  root_out+spec+aux_path+aux_path_out

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

% data_stream=[0,1,2]
% if data_stream = 0, no output
% if data_stream = 1, binaries generated
% if data_stream = 2, text files generated


%tic;
path_in=strcat(root,spec,aux_path);

cd('../../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,filename);

cd('../1dproj');

if ismember(0,data_stream)    
   [proj1d] = proj1d_dm_data_out( root,root_data,spec,aux_path,aux_path_data,filename,lenght_factor,resol_factor,pivot,rot_angle,0);
else    
    path_orig_data=strcat(strcat(root_data,spec,aux_path),'data/',aux_path_data,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');
    path_out=strcat(path_orig_data,'wavelet/dc/');
    mkdir(path_orig_data,strcat('wavelet/dc/'));
    mkdir(path_orig_data,strcat('wavelet/dc/scale/'));
    if ismember(2,data_stream)
        if ismember(3,data_stream)
            proj1d_dm_data_out( root,root_data,spec,aux_path,aux_path_data,filename,lenght_factor,resol_factor,pivot,rot_angle,2);
        end
        proj1d=dlmread(char(strcat(path_orig_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.txt')));    
    end
    if ismember(1,data_stream)
        if ismember(3,data_stream)
            proj1d_dm_data_out( root,root_data,spec,aux_path,aux_path_data,filename,lenght_factor,resol_factor,pivot,rot_angle,1);
        end
        fileID = fopen(strcat(path_orig_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.bin'));
        proj1d=fread(fileID,'float32','l');
        fclose(fileID);
    end
end

average=mean2(proj1d);
proj1d=(proj1d-average)/average;   

[dc_cwt,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);



if ~ismember(0,data_stream)
    if ismember(1,data_stream)        
        fileID1 = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_data.bin'),'w');
        fwrite(fileID1,dc_cwt, 'float32','l');
        fclose(fileID1);
        fileID2 = fopen(strcat(path_out,'scale/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_scale.bin'),'w');
        fwrite(fileID2,seconds(periods), 'float32','l');
        fclose(fileID2);
    end

    if ismember(2,data_stream)
        dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_data.txt'),dc_cwt,'delimiter','\t');
        dlmwrite(strcat(path_out,'scale/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_scale.txt'),seconds(periods),'delimiter','\t');
    end
    
end

%proj1d with the cutoff

if resol_factor>=1
    low_pass=resol_factor;
else
    low_pass=1;
end

i_dc_cwt = icwt(dc_cwt,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
     
if ~ismember(0,data_stream)
    mkdir(path_orig_data,strcat('wavelet/dc/filter_1dproj_',num2str(cutoff),'MpcCut/'));
    if ismember(1,data_stream)        
        fileID1 = fopen(strcat(path_orig_data,strcat('wavelet/dc/filter_1dproj_',num2str(cutoff),'MpcCut/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.bin'),'w');
        fwrite(fileID1,i_dc_cwt, 'float32','l');
        fclose(fileID1);
    end

    if ismember(2,data_stream)
        dlmwrite(strcat(path_orig_data,strcat('wavelet/dc/filter_1dproj_',num2str(cutoff),'MpcCut/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.txt'),i_dc_cwt,'delimiter','\t');
    end
    
end

%filtered wavelet coeficients

[filtered_dc_cwt,periods] = cwt(i_dc_cwt,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);

if ~ismember(0,data_stream)
    mkdir(path_orig_data,strcat('wavelet/dc/filtered_cwt_',num2str(cutoff),'MpcCut/'));
    if ismember(1,data_stream)        
        fileID1 = fopen(strcat(path_orig_data,strcat('wavelet/dc/filtered_cwt_',num2str(cutoff),'MpcCut/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.bin'),'w');
        fwrite(fileID1,filtered_dc_cwt, 'float32','l');
        fclose(fileID1);
    end
    if ismember(2,data_stream)
        dlmwrite(strcat(path_orig_data,strcat('wavelet/dc/filtered_cwt_',num2str(cutoff),'MpcCut/'),'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_wavelet_filter_data.txt'),filtered_dc_cwt,'delimiter','\t');
    end
    
end

     

cd('../wake_detection/lets');

% 
% %toc;



end

