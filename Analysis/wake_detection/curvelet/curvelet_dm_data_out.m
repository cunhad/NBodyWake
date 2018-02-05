function [ proj1d_angles] = curvelet_dm_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor_in,resol_factor_in,lenght_factor_out,resol_factor_out,pivot_in,NSIDE_in,part_in,num_cores_in,NSIDE_out,part_out,num_cores_out,data_stream,cutoff)

%(example)  curvelet_dm_per_node_part('/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',4,0,1,1);
%(example)  for i=1:8; curvelet_dm_per_node_part('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',4,0,8,i); end;

%(example)  [proj1d_angles,filtered_proj1d_angles,out_proj1d_angles,out_dc_proj1d_angles,out_filtered_proj1d_angles,out_filtered_dc_proj1d_angles,out_dc_filtered_proj1d_angles] = curvelet_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/data_test/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',2,1,[0,0,0],4,1,4,[1,2],10);
%(example) curvelet_dm_data_out('/home/asus/Dropbox/extras/storage/guillimin/', '/home/asus/Dropbox/extras/storage/guillimin/box_stat/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','10.000xv0.dat',2,2,4,2,[0,0,0],64,16,4,1,16,4,[1,2],0.4);


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filter=[0,1,2,3]

% if filter = 0 , no wavelet filter
% if filter = 1
% if filter = 2

myCluster = parcluster('local');
myCluster.NumWorkers=num_cores_out;
saveProfile(myCluster);

p = parpool(num_cores_out);
tic;

cd('../../preprocessing');

%  [ ~,redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
 [~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

angles(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE_out),'_t.cvs'));
angles(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE_out),'_p.cvs'));
[~,number_of_angle_nuple] = size(angles);
%
% mkdir(root_out);
% mkdir(root_out,strcat(spec,aux_path));


[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,part_out,1);

root_out=strcat(root_out(1,1:end-1),'_curv/');

path_data_in_curv=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor_in),'lf_',num2str(resol_factor_in),'rf_',strcat(num2str(pivot_in(1)),'-',num2str(pivot_in(2)),'-',num2str(pivot_in(3))),'pv','/','stat/box_statistics/dm/');
peaks=dlmread(strcat(path_data_in_curv,'_',num2str(find(str2num(char(redshift_list))==z)),'_peaks_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_in),'_NSIDE',num2str(NSIDE_in),'.txt'));


distance_to_center=peaks(:,2);
theta=peaks(:,3);
phi=peaks(:,4);

pivot(1)=distance_to_center(4)*sin(theta(4))*cos(phi(4));
pivot(2)=distance_to_center(4)*sin(theta(4))*cos(phi(4));
pivot(3)=distance_to_center(4)*cos(theta(4));

bins=[-(nc/(2*lenght_factor_out)):nc/(np*resol_factor_out):(nc/(2*lenght_factor_out))];
proj1d_angles=zeros(number_of_angle_nuple,length(bins)-1);

%for k = 3:-1:1
for part_id = 1  :   part_out
    %for k = 1  : 1
    
    cd('../preprocessing');
    
    [ ~, nc, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part_out,part_id);
    
    Pos=mod(Pos,nc);
    
    
    
    parfor i=1:number_of_angle_nuple
        %     for i=1:number_of_angle_nuple
        %   for i=1:1
        
        
        
        theta=angles(1,i);
        phi=angles(2,i);
        
%         theta=mod(theta+pi/2,pi)
%         phi=mod(phi+pi/2,2*pi)
        
%         theta=theta-pi/2;
%         phi=phi-pi/4;
        
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
        
        liminf=-(1/(2*lenght_factor_out))*nc;
        limsup= (1/(2*lenght_factor_out))*nc;
        
        liminf_slice=-4;
        limsup_slice= 4;
        
        conditionsx=rx(1,:)<=liminf|rx(1,:)>=limsup;
        conditionsy=rx(2,:)<=liminf|rx(2,:)>=limsup;
        conditionsz=rx(3,:)<=liminf_slice|rx(3,:)>=limsup_slice;
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
        
        fprintf('done for z= %f,part=%d and  i= %d\n',z,part_id, i);
        %display(proj);
        
    end
end

proj1d_angles=transpose(proj1d_angles);
max_proj1d_angles=max(proj1d_angles);
average_proj1d_angles=mean(proj1d_angles,1);
max_amplitude_proj1d_angles=max_proj1d_angles(:)-average_proj1d_angles(:);

proj1d_index_max=zeros(1,number_of_angle_nuple);

parfor angl=1:number_of_angle_nuple
proj1d_index_max(1,angl)=find(proj1d_angles(:,angl)==max_proj1d_angles(1,angl),1);
end
proj1d_angles_snremoved=proj1d_angles;

for angl=1:number_of_angle_nuple
proj1d_angles_snremoved(proj1d_index_max(1,angl),angl)=average_proj1d_angles(angl);
end



std_proj1d_angles=std(proj1d_angles_snremoved);

clearvars proj1d_angles_snremoved;

stn_proj1d_angles=(max_amplitude_proj1d_angles(:))./std_proj1d_angles(:);


out_proj1d_angles=[transpose(max_amplitude_proj1d_angles);std_proj1d_angles;transpose(stn_proj1d_angles)];

parfor angl=1:number_of_angle_nuple
dc_proj1d_angles(:,angl)=(proj1d_angles(:,angl)-average_proj1d_angles(angl))./average_proj1d_angles(angl);
end

max_dc_proj1d_angles=max(dc_proj1d_angles);
max_amplitude_dc_proj1d_angles=max_dc_proj1d_angles(:);

dc_proj1d_index_max=zeros(1,number_of_angle_nuple);

parfor angl=1:number_of_angle_nuple
dc_proj1d_index_max(1,angl)=find(dc_proj1d_angles(:,angl)==max_dc_proj1d_angles(1,angl),1);
end
dc_proj1d_angles_snremoved=dc_proj1d_angles;

for angl=1:number_of_angle_nuple
dc_proj1d_angles_snremoved(dc_proj1d_index_max(1,angl),angl)=0;
end

std_dc_proj1d_angles=std(dc_proj1d_angles_snremoved);
clearvars dc_proj1d_angles_snremoved;
stn_dc_proj1d_angles=(max_amplitude_dc_proj1d_angles(:))./std_dc_proj1d_angles(:);

out_dc_proj1d_angles=[transpose(max_amplitude_dc_proj1d_angles);std_dc_proj1d_angles;transpose(stn_dc_proj1d_angles)];

if cutoff~=0        
    filtered_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple);
    filtered_dc_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple);
    
    if resol_factor_out>=1
        low_pass=resol_factor_out;
    else
        low_pass=1;
    end    
    parfor i=1:number_of_angle_nuple
%         average_proj1d_angles=mean2(proj1d_angles(:,i));
%         proj1d=(proj1d_angles(:,i)-average_proj1d_angles)/average_proj1d_angles;
        proj1d=proj1d_angles(:,i);
        [cwt_proj1d,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor_out)),'waveletparameters',[3 3.01]);
        filtered_proj1d_angles(:,i) = icwt(cwt_proj1d,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
        
        dc_proj1d=dc_proj1d_angles(:,i);
        [cwt_dc_proj1d,periods] = cwt(dc_proj1d,seconds(size_box/(np*resol_factor_out)),'waveletparameters',[3 3.01]);
        filtered_dc_proj1d_angles(:,i) = icwt(cwt_dc_proj1d,periods,[periods(low_pass) seconds(cutoff)],'waveletparameters',[3 3.01]);
        
    end
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
    max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
    
    filtered_proj1d_index_max=zeros(1,number_of_angle_nuple);
    
    parfor angl=1:number_of_angle_nuple
        filtered_proj1d_index_max(1,angl)=find(filtered_proj1d_angles(:,angl)==max_filtered_proj1d_angles(1,angl),1);
    end
    filtered_proj1d_angles_snremoved=filtered_proj1d_angles;
    
    for angl=1:number_of_angle_nuple
        filtered_proj1d_angles_snremoved(filtered_proj1d_index_max(1,angl),angl)=average_filtered_proj1d_angles(angl);
    end
            
    std_filtered_proj1d_angles=std(filtered_proj1d_angles_snremoved);
    clearvars filtered_proj1d_angles_snremoved;
    stn_filtered_proj1d_angles=(max_amplitude_filtered_proj1d_angles(:))./std_filtered_proj1d_angles(:);
    
    out_filtered_proj1d_angles=[transpose(max_amplitude_filtered_proj1d_angles);std_filtered_proj1d_angles;transpose(stn_filtered_proj1d_angles)];
    
%     parfor angl=1:number_of_angle_nuple
%         dc_filtered_proj1d_angles(:,angl)=(filtered_proj1d_angles(:,angl)-average_filtered_proj1d_angles(angl))./average_filtered_proj1d_angles(angl);
%     end
    
%     max_dc_filtered_proj1d_angles=max(dc_filtered_proj1d_angles);
%     max_amplitude_dc_filtered_proj1d_angles=max_dc_filtered_proj1d_angles(:);
    
%     max_amplitude_dc_filtered_proj1d_angles=max(dc_filtered_proj1d_angles);

    
%     dc_filtered_proj1d_index_max=zeros(1,number_of_angle_nuple);
    
%     parfor angl=1:number_of_angle_nuple
%         dc_filtered_proj1d_index_max(1,angl)=find(dc_filtered_proj1d_angles(:,angl)==max_amplitude_dc_filtered_proj1d_angles(1,angl));
%     end
%     dc_filtered_proj1d_angles_snremoved=dc_filtered_proj1d_angles;
%     
%     for angl=1:number_of_angle_nuple
%         dc_filtered_proj1d_angles_snremoved(dc_filtered_proj1d_index_max(1,angl),angl)=0;
%     end
%        
    
%     std_dc_filtered_proj1d_angles=std(dc_filtered_proj1d_angles_snremoved);
%     clearvars dc_filtered_proj1d_angles_snremoved;
%     stn_dc_filtered_proj1d_angles=(max_amplitude_dc_filtered_proj1d_angles(:))./std_dc_filtered_proj1d_angles(:);
%     
%     out_dc_filtered_proj1d_angles=[(max_amplitude_dc_filtered_proj1d_angles);std_dc_filtered_proj1d_angles;transpose(stn_dc_filtered_proj1d_angles)];
%     
    max_filtered_dc_proj1d_angles=max(filtered_dc_proj1d_angles);
    average_filtered_dc_proj1d_angles=mean(filtered_dc_proj1d_angles,1);
    
    filtered_dc_proj1d_index_max=zeros(1,number_of_angle_nuple);
    
    parfor angl=1:number_of_angle_nuple
        filtered_dc_proj1d_index_max(1,angl)=find(filtered_dc_proj1d_angles(:,angl)==max_filtered_dc_proj1d_angles(1,angl),1);
    end
    filtered_dc_proj1d_angles_snremoved=filtered_dc_proj1d_angles;
    
    for angl=1:number_of_angle_nuple
        filtered_dc_proj1d_angles_snremoved(filtered_dc_proj1d_index_max(1,angl),angl)=0;
    end
        
    max_amplitude_filtered_dc_proj1d_angles=max_filtered_dc_proj1d_angles(:)-average_filtered_dc_proj1d_angles(:);
    std_filtered_dc_proj1d_angles=std(filtered_dc_proj1d_angles_snremoved);
    clearvars filtered_dc_proj1d_angles_snremoved;
    
    stn_filtered_dc_proj1d_angles=(max_amplitude_filtered_dc_proj1d_angles(:))./std_filtered_dc_proj1d_angles(:);
    
    out_filtered_dc_proj1d_angles=[transpose(max_amplitude_filtered_dc_proj1d_angles);std_filtered_dc_proj1d_angles;transpose(stn_filtered_dc_proj1d_angles)];
    
end



if ~ismember(0,data_stream)
    
    path_out_all=strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor_out),'lf_',num2str(resol_factor_out),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor_out),'lf_',num2str(resol_factor_out),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    mkdir(strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor_out),'lf_',num2str(resol_factor_out),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor_out),'lf_',num2str(resol_factor_out),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor_out),'lf_',num2str(resol_factor_out),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc/'));
    
    if cutoff~=0 
        mkdir(path_out,strcat('/cutoff_',num2str(cutoff),'MpcCut/'));
        mkdir(path_out,strcat('/cutoff_',num2str(cutoff),'MpcCut/dc/'));
        mkdir(path_out_all,strcat('/cutoff_',num2str(cutoff),'MpcCut/'));
        mkdir(path_out_all,strcat('/cutoff_',num2str(cutoff),'MpcCut/dc/'));
%         mkdir(path_out,strcat('/dc/cutoff_',num2str(cutoff),'MpcCut/'));
    end
    
    if ismember(1,data_stream)

        fileID = fopen(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
        fwrite(fileID,proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
        fwrite(fileID,dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
        fwrite(fileID,out_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
        fwrite(fileID,out_dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        if cutoff~=0 
            
            fileID = fopen(strcat(path_out_all,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
            fwrite(fileID,filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out_all,'cutoff_',num2str(cutoff),'MpcCut/dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
            fwrite(fileID,filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
            fwrite(fileID,out_filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.bin'),'w');
            fwrite(fileID,out_filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
%             fileID = fopen(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
%             fwrite(fileID,out_dc_filtered_proj1d_angles, 'float32','l');
%             fclose(fileID);
            
        end
    end
    
    if ismember(2,data_stream)
            dlmwrite(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),out_dc_proj1d_angles,'delimiter','\t');
        
        if cutoff~=0 
            dlmwrite(strcat(path_out_all,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),filtered_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'cutoff_',num2str(cutoff),'MpcCut/dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),filtered_dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),out_proj1d_angles,'delimiter','\t');
%             dlmwrite(strcat(path_out,'dc/','cutoff_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'cutoff_',num2str(cutoff),'MpcCut/dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part_out),'_NSIDE',num2str(NSIDE_out),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
        end
        
    end
    
end


cd('../wake_detection/curvelet');

toc;
delete(gcp('nocreate'))


end

