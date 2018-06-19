function [  root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,particl_part,angle_part,angle_p,num_cores,data_stream,level_window,dwbasis ] = dm_box_cubic_fast_files_to_file(  root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,particl_part,angle_part,angle_p,num_cores,data_stream,level_window,dwbasis )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%(example)  [] = dm_box_cubic_fast_files_to_file('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','10.000xv0.dat',2,1,[0,0,0],64,1,1,4,4,[1,2],[1],'sym6');
% 

angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));

cd ../preprocessing

[~,redshift_list,~,size_box,nc,np,~,~,~,Gmu,~] = preprocessing_info(root,spec,aux_path );
[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);

bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];
[~,number_of_angle_nuple_hpx] = size(angles_hpx);


if ~ismember(0,data_stream)
    
    path_out_all=strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/');
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/');
    
    if ismember(1,data_stream)
        
        %         fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        %         fwrite(fileID,out_proj1d_angles, 'float32','l');
        %         fclose(fileID);
        
        files_out_proj1d_angles_list=dir(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid*_NSIDE',num2str(NSIDE),'.bin'));
        files_out_proj1d_angles_list={files_out_proj1d_angles_list.name};
        
        fid_o = fopen(strcat(path_out,files_out_proj1d_angles_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'),'w');
        
        for file_idx=1:angle_p
            filename=strcat(path_out,files_out_proj1d_angles_list{file_idx});
            fid = fopen(filename);
            %             fread(fid, 'float32','l');
            %             fseek(fid_o, 0, 'eof');
            fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
            fclose(fid);
        end
        fclose(fid_o);
        
        fid = fopen(strcat(path_out,files_out_proj1d_angles_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'));
        fid_o = fopen(strcat(path_out,files_out_proj1d_angles_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        
        fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
        
        for angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx
            
            theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));
            theta_inverse_indices=number_of_angle_nuple_hpx-theta_indices+1;
            phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
            fseek(fid,(3*(phi_indice_inverse-1)*4),'bof');
            fwrite(fid_o,fread(fid,3, 'float32','l'),'float32','l');
        end
        
        fclose(fid);
        fclose(fid_o);
        
        
        
        %
        %         fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        %         fwrite(fileID,out_dc_proj1d_angles, 'float32','l');
        %         fclose(fileID);
        %
        if level_window~=0
            
            %             fileID = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            %             fwrite(fileID,out_filtered_proj1d_angles, 'float32','l');
            %             fclose(fileID);
            
            files_out_filtered_1dproj_angle_z_list=dir(strcat(path_out,'level_window',mat2str(level_window(:)),'/_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid*_NSIDE',num2str(NSIDE),'.bin'));
            files_out_filtered_1dproj_angle_z_list={files_out_filtered_1dproj_angle_z_list.name};
            
            fid_o = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/',files_out_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'),'w');
            
            for file_idx=1:angle_p
                filename=strcat(path_out,'level_window',mat2str(level_window(:)),'/',files_out_filtered_1dproj_angle_z_list{file_idx});
                fid = fopen(filename);
                %             fread(fid, 'float32','l');
                %             fseek(fid_o, 0, 'eof');
                fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
                fclose(fid);
            end
            fclose(fid_o);
            
            fid = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/',files_out_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'));
            fid_o = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/',files_out_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
            
            for angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx
                
                theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));
                theta_inverse_indices=number_of_angle_nuple_hpx-theta_indices+1;
                phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
                fseek(fid,(3*(phi_indice_inverse-1)*4),'bof');
                fwrite(fid_o,fread(fid,3, 'float32','l'),'float32','l');
            end
            
            fclose(fid);
            fclose(fid_o);
           
%             
%             fileID = fopen(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%             fwrite(fileID,out_filtered_dc_proj1d_angles, 'float32','l');
%             fclose(fileID);
%             
%             fileID = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%             fwrite(fileID,filtered_proj1d_angles, 'float32','l');
%             fclose(fileID);
%             

            files_filtered_1dproj_angle_z_list=dir(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid*_NSIDE',num2str(NSIDE),'.bin'));
            files_filtered_1dproj_angle_z_list={files_filtered_1dproj_angle_z_list.name};
            
            fid_o = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/',files_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'),'w');
            
            for file_idx=1:angle_p
                filename=strcat(path_out_all,'level_window',mat2str(level_window(:)),'/',files_filtered_1dproj_angle_z_list{file_idx});
                fid = fopen(filename);
                %             fread(fid, 'float32','l');
                %             fseek(fid_o, 0, 'eof');
                fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
                fclose(fid);
            end
            fclose(fid_o);
            
            fid = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/',files_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'_half.bin'));
            fid_o = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/',files_filtered_1dproj_angle_z_list{1}(1:end-20),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fid_o,fread(fid, 'float32','l'),'float32','l');
            
            for angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx
                
                theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));
                theta_inverse_indices=number_of_angle_nuple_hpx-theta_indices+1;
                phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
                fseek(fid,((length(bins)-1)*(phi_indice_inverse-1)*4),'bof');
                fwrite(fid_o,fread(fid,(length(bins)-1), 'float32','l'),'float32','l');
            end
            
            fclose(fid);
            fclose(fid_o);

%             fileID = fopen(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%             fwrite(fileID,filtered_dc_proj1d_angles, 'float32','l');
%             fclose(fileID);
%             
        end
%         
%         fileID = fopen(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%         fwrite(fileID,proj1d_angles, 'float32','l');
%         fclose(fileID);
%         
%         fileID = fopen(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z','_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%         fwrite(fileID,dc_proj1d_angles, 'float32','l');
%         fclose(fileID);
%         
    end
%     
%     if ismember(2,data_stream)
%             dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_dc_proj1d_angles,'delimiter','\t');
%         
%         if level_window~=0 
%             dlmwrite(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),filtered_proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),filtered_dc_proj1d_angles,'delimiter','\t');
% 
%         end
%         
%             dlmwrite(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),dc_proj1d_angles,'delimiter','\t');
% 
%         
%     end
%     
end

cd ../processing


end