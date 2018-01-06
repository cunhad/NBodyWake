function [ stn_proj1d_angles,stn_filtered_proj1d_angles] = box_statistics_dm_data_out_fo( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part_mem_per_outputfile,part_fileout_number,part_fileout_id,num_cores,data_stream,filter,cutoff)

%(example)  box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',4,0,1,1);
%(example)  for i=1:8; box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',4,0,8,i); end;

%(example)  [stn_proj1d_angles,stn_filtered_proj1d_angles] = box_statistics_dm_data_out_fo('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/data_test/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',2,1,[0,0,0],4,1,2,1,4,[1,2],[0,1],10);

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filter=[0,1,2,3]

% if filter = 0 , no wavelet filter
% if filter = 1
% if filter = 2



myCluster = parcluster('local');
myCluster.NumWorkers=num_cores;
saveProfile(myCluster);


p = parpool(num_cores);
tic;

part=part_mem_per_outputfile*part_fileout_number;

cd('../../preprocessing');

 [ ~,redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

% z_string=filename(1:strfind(filename,'xv')-1);
% 
% 
%     path_in=strcat(root,spec,aux_path);
%     file_name_list = dir(strcat(path_in,z_string,'xv*','.dat'));
%     file_name_list={file_name_list.name};
%
% cd('../processing');
%
% file_name_list=sort_nat(file_name_list);
%
% cd('../preprocessing');
%
angles = dlmread(strcat('../../python/angles',num2str(NSIDE),'.txt'));
[~,number_of_angle_nuple] = size(angles);
%
% mkdir(root_out);
% mkdir(root_out,strcat(spec,aux_path));
[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,part,1);
bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];
proj1d_angles=zeros(number_of_angle_nuple,length(bins)-1);

if ismember(0,filter)
    %for k = 3:-1:1
    
    for part_men_per_out_id= 1 : part_mem_per_outputfile
        
        part_id = part_mem_per_outputfile*(part_fileout_id-1)+part_men_per_out_id;
        %for k = 1  : 1
        
        cd('../preprocessing');
        
        [ ~, nc, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part,part_id);
        
        Pos=mod(Pos,nc);
        
        
        
        parfor i=1:number_of_angle_nuple
            %     for i=1:number_of_angle_nuple
            %   for i=1:1
            
            
            
            theta=angles(1,i);
            phi=angles(2,i);
            
            hist1d_cor=zeros(1,length(bins)-1);
            
            % for j=1:1000
            
            rx=[];
            
            rx(1,:)=Pos(1,:)-(nc/2)-pivot(1);
            rx(2,:)=Pos(2,:)-(nc/2)-pivot(2);
            rx(3,:)=Pos(3,:)-(nc/2)-pivot(3);
            
            Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
            Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
            %
            rx=Rz*rx;
            rx=Ry*rx;
            
            liminf=-(1/(2*lenght_factor))*nc;
            limsup= (1/(2*lenght_factor))*nc;
            conditionsx=rx(1,:)<=liminf|rx(1,:)>=limsup;
            conditionsy=rx(2,:)<=liminf|rx(2,:)>=limsup;
            conditionsz=rx(3,:)<=liminf|rx(3,:)>=limsup;
            conditions=conditionsx|conditionsy|conditionsz;
            rx(:,conditions)=[];
            
            rx=transpose(rx);
            
            %display(rx);
            
            if(~isempty(rx))
                
                [count, ~ ,~ ,~] = histcn(rx,1,1,bins);
                % display(count);
                % display(length(bins));
                count=count(1:1,1:1,1:length(bins)-1);
                %     average=mean2(count);
                %     count=(count-average)/average;
                count=squeeze(count);
                count=squeeze(count);
                
                count=transpose(count);
                
                % display(count);
                
                
                
                hist1d_cor=count;
                
                % display(hist1d_cor);
                
                
                
                
                
            end
            
            % proj=transpose(proj);
            
            %display(proj);
            
            % [proj edges mid loc] = histcn(proj,bins);
            
            
            %proj=transpose(proj);
            
            % display(size(proj1d_angles));
            % display(size(hist1d_cor));
            
            % average=mean2(hist1d_cor);
            % hist1d_cor=(hist1d_cor-average)/average;
            
            proj1d_angles(i,:)=proj1d_angles(i,:)+hist1d_cor(1,:);
            
            fprintf('done for z= %f and  i= %d\n',z, i);
            %display(proj);
            
        end
    end
    
end

proj1d_angles=transpose(proj1d_angles);

max_proj1d_angles=max(proj1d_angles);
average=mean(proj1d_angles,1);

std_proj1d_angles=std(proj1d_angles);
stn_proj1d_angles=(max_proj1d_angles(:)-average(:))./std_proj1d_angles(:);

if ismember(1,filter)
    
    filtered_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple);
    
    if resol_factor>=1
        low_pass=resol_factor;
    else
        low_pass=1;
    end
    
    parfor i=1:number_of_angle_nuple
        %average=mean2(proj1d_angles(:,i));
        %proj1d=(proj1d_angles(:,i)-average)/average;
        proj1d=proj1d_angles(:,i);
        [dc_cwt,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
        filtered_proj1d_angles(:,i) = icwt(dc_cwt,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
    end
end

if ismember(1,filter)
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    std_filtered_proj1d_angles=std(filtered_proj1d_angles);
    stn_filtered_proj1d_angles=max_filtered_proj1d_angles(:)./std_filtered_proj1d_angles(:);
else
    stn_filtered_proj1d_angles=zeros(1,number_of_angle_nuple);
end

if ~ismember(0,data_stream)
    
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/'));
    
    mkdir(path_out,'box_1dproj/');
    mkdir(path_out,strcat('box_1dproj/cutoff_',num2str(cutoff),'MpcCut/'));

    
    if ismember(1,data_stream)
%         fileID = fopen(strcat(path_out,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%         fwrite(fileID,proj1d_angles, 'float32','l');
%         fclose(fileID);
        if ismember(0,filter)

        fileID = fopen(strcat(path_out,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dproj_angle_z',num2str(z),'_fon',num2str(part_fileout_number),'_foid',num2str(part_fileout_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,stn_proj1d_angles, 'float32','l');
        fclose(fileID);
        end
        if ismember(1,filter)
        fileID = fopen(strcat(path_out,'box_1dproj/cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dwproj_angle_z',num2str(z),'_fon',num2str(part_fileout_number),'_foid',num2str(part_fileout_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,stn_filtered_proj1d_angles, 'float32','l');
        fclose(fileID);
        end
    end
    
    if ismember(2,data_stream)
        if ismember(0,filter)
        %dlmwrite(strcat(path_out,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
        dlmwrite(strcat(path_out,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dproj_angle_z',num2str(z),'_fon',num2str(part_fileout_number),'_foid',num2str(part_fileout_id),'_NSIDE',num2str(NSIDE),'.txt'),stn_proj1d_angles,'delimiter','\t');
        end
        if ismember(1,filter)
        %dlmwrite(strcat(path_out,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
        dlmwrite(strcat(path_out,'box_1dproj/cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dwproj_angle_z',num2str(z),'_fon',num2str(part_fileout_number),'_foid',num2str(part_fileout_id),'_NSIDE',num2str(NSIDE),'.txt'),stn_proj1d_angles,'delimiter','\t');
        end
    end
    
end

%     fileID = fopen(strcat(path_out,'npart_per_node_1dproj/','1dproj_angle_z',num2str(z),'_node',num2str(node),'_partID',num2str(part_id),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%     fwrite(fileID,proj1d_angles, 'float32','l');
%     fclose(fileID);

% display(strcat(path_out,'npart_per_node_1dproj/','1dproj_angle_z',num2str(z),'_node',num2str(node),'_partID',num2str(part_id),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));


cd('../wake_detection/box_statistics');

toc;
delete(gcp('nocreate'))


end

