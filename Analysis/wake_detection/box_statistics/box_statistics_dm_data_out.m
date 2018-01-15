function [ proj1d_angles,filtered_proj1d_angles,out_proj1d_angles,out_dc_proj1d_angles,out_filtered_proj1d_angles,out_filtered_dc_proj1d_angles,out_dc_filtered_proj1d_angles] = box_statistics_dm_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,data_stream,cutoff)

%(example)  box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',4,0,1,1);
%(example)  for i=1:8; box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',4,0,8,i); end;

%(example)  [proj1d_angles,filtered_proj1d_angles,out_proj1d_angles,out_dc_proj1d_angles,out_filtered_proj1d_angles,out_filtered_dc_proj1d_angles,out_dc_filtered_proj1d_angles] = box_statistics_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/data_test/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',2,1,[0,0,0],4,1,4,[1,2],10);

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

cd('../../preprocessing');

 [ ~,redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

angles = dlmread(strcat('../../python/angles',num2str(NSIDE),'.txt'));
[~,number_of_angle_nuple] = size(angles);
%
% mkdir(root_out);
% mkdir(root_out,strcat(spec,aux_path));
[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,part,1);
bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];
proj1d_angles=zeros(number_of_angle_nuple,length(bins)-1);

%for k = 3:-1:1
for part_id = 1  :   part
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

proj1d_angles=transpose(proj1d_angles);
max_proj1d_angles=max(proj1d_angles);
average_proj1d_angles=mean(proj1d_angles,1);
max_amplitude_proj1d_angles=max_proj1d_angles(:)-average_proj1d_angles(:);
std_proj1d_angles=std(proj1d_angles);
stn_proj1d_angles=(max_amplitude_proj1d_angles(:))./std_proj1d_angles(:);

out_proj1d_angles=[transpose(max_amplitude_proj1d_angles);std_proj1d_angles;transpose(stn_proj1d_angles)];

parfor angl=1:number_of_angle_nuple
dc_proj1d_angles(:,angl)=(proj1d_angles(:,angl)-average_proj1d_angles(angl))./average_proj1d_angles(angl);
end

max_dc_proj1d_angles=max(dc_proj1d_angles);
max_amplitude_dc_proj1d_angles=max_dc_proj1d_angles(:);
std_dc_proj1d_angles=std(dc_proj1d_angles);
stn_dc_proj1d_angles=(max_amplitude_dc_proj1d_angles(:))./std_dc_proj1d_angles(:);

out_dc_proj1d_angles=[transpose(max_amplitude_dc_proj1d_angles);std_dc_proj1d_angles;transpose(stn_dc_proj1d_angles)];

if cutoff~=0        
    filtered_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple);
    filtered_dc_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple);
    
    if resol_factor>=1
        low_pass=resol_factor;
    else
        low_pass=1;
    end    
    parfor i=1:number_of_angle_nuple
%         average_proj1d_angles=mean2(proj1d_angles(:,i));
%         proj1d=(proj1d_angles(:,i)-average_proj1d_angles)/average_proj1d_angles;
        proj1d=proj1d_angles(:,i);
        [cwt_proj1d,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
        filtered_proj1d_angles(:,i) = icwt(cwt_proj1d,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
        
        dc_proj1d=dc_proj1d_angles(:,i);
        [cwt_dc_proj1d,periods] = cwt(dc_proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
        filtered_dc_proj1d_angles(:,i) = icwt(cwt_dc_proj1d,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
        
    end
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
    max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
    std_filtered_proj1d_angles=std(filtered_proj1d_angles);
    stn_filtered_proj1d_angles=(max_amplitude_filtered_proj1d_angles(:))./std_filtered_proj1d_angles(:);
    
    out_filtered_proj1d_angles=[transpose(max_amplitude_filtered_proj1d_angles);std_filtered_proj1d_angles;transpose(stn_filtered_proj1d_angles)];
    
    parfor angl=1:number_of_angle_nuple
        dc_filtered_proj1d_angles(:,angl)=(filtered_proj1d_angles(:,angl)-average_filtered_proj1d_angles(angl))./average_filtered_proj1d_angles(angl);
    end
    
    max_amplitude_dc_filtered_proj1d_angles=max(dc_filtered_proj1d_angles);
    std_dc_filtered_proj1d_angles=std(dc_filtered_proj1d_angles);
    stn_dc_filtered_proj1d_angles=(max_amplitude_dc_filtered_proj1d_angles(:))./std_dc_filtered_proj1d_angles(:);
    
    out_dc_filtered_proj1d_angles=[(max_amplitude_dc_filtered_proj1d_angles);std_dc_filtered_proj1d_angles;transpose(stn_dc_filtered_proj1d_angles)];
    
    max_filtered_dc_proj1d_angles=max(filtered_dc_proj1d_angles);
    average_filtered_dc_proj1d_angles=mean(filtered_dc_proj1d_angles,1);
    max_amplitude_filtered_dc_proj1d_angles=max_filtered_dc_proj1d_angles(:)-average_filtered_dc_proj1d_angles(:);
    std_filtered_dc_proj1d_angles=std(filtered_dc_proj1d_angles);
    stn_filtered_dc_proj1d_angles=(max_amplitude_filtered_dc_proj1d_angles(:))./std_filtered_dc_proj1d_angles(:);
    
    out_filtered_dc_proj1d_angles=[transpose(max_amplitude_filtered_dc_proj1d_angles);std_filtered_dc_proj1d_angles;transpose(stn_filtered_dc_proj1d_angles)];
    
end



if ~ismember(0,data_stream)
    
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc/'));
    
    if cutoff~=0 
        mkdir(path_out,strcat('/cutoff_',num2str(cutoff),'MpcCut/'));
        mkdir(path_out,strcat('/cutoff_',num2str(cutoff),'MpcCut/dc/'));
        mkdir(path_out,strcat('/dc/cutoff_',num2str(cutoff),'MpcCut/'));
    end
    
    if ismember(1,data_stream)

        fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        if cutoff~=0 
            
            fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'dc/','cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_dc_filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
        end
    end
    
    if ismember(2,data_stream)
            dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_dc_proj1d_angles,'delimiter','\t');
        
        if cutoff~=0 
            dlmwrite(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),filtered_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_dc_filtered_proj1d_angles,'delimiter','\t');
        end
        
    end
    
end


cd('../wake_detection/box_statistics');

toc;
delete(gcp('nocreate'))


end

