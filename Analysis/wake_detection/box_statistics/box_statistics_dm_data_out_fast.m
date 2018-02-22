function [ proj1d_angles ] = box_statistics_dm_data_out_fast( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,data_stream,level_window,dwbasis)

%(example)  box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',4,0,1,1);
%(example)  for i=1:8; box_statistics_dm_per_node_part('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',4,0,8,i); end;

%(example)  [proj1d_angles] = box_statistics_dm_data_out_fast('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_fast/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','10.000xv0.dat',2,1,[0,0,0],2,1,4,[1,2],[1],'sym6');

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filter=[0,1,2,3]

% if filter = 0 , no wavelet filter
% if filter = 1
% if filter = 2

% %passivetransf melhorou1/3
% 
myCluster = parcluster('local');
myCluster.NumWorkers=num_cores;
saveProfile(myCluster);

p = parpool(num_cores);



tic;
% 
cd('../../preprocessing');

[~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));

% angles_hpx_t = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
% angles_hpx_p= dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));


[~,number_of_angle_nuple_hpx] = size(angles_hpx);

number_of_angle_nuple_hpx=number_of_angle_nuple_hpx/2;

n_angle_per_node=ceil(number_of_angle_nuple_hpx/num_cores);

for cr=1:num_cores+1
    angl_indx(cr)= n_angle_per_node*(cr-1)+1;
end
for cr=1:num_cores
    n_angl_indx(cr)= -angl_indx(cr)+angl_indx(cr+1);
    
end

for cr=1:num_cores
     angl_ind_start=angl_indx(cr);
     angl_ind_end=angl_indx(cr+1)-1;
     angl_chunk_t{cr}=angles_hpx(1,angl_ind_start:angl_ind_end);
     angl_chunk_p{cr}=angles_hpx(2,angl_ind_start:angl_ind_end);   
end

clearvars angles_hpx;

% N_particles=0;
% 
% for part_pre_id = 1 : part_pre
%     [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,~, ~, ~, n_size_id ] = preprocessing_part_aux(root,spec,aux_path,filename,lenght_factor,pivot,part_pre,part_pre_id);
%     N_particles=N_particles+n_size_id;
% end
% 
% N_particles_chunk_exact=N_particles/part;




[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,part,1);
bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];
proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);

% histogr_1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);

% % part_id=1;

% while part_id <   part +1
for part_id = 1  :   part

    
    cd('../preprocessing');
    
    [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part,part_id);
    
        Pos=mod(Pos,nc);
 
        Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
        Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
        Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
        
%         lenght_factor_pre=lenght_factor/2;

        
        lim_pre=(1/(lenght_factor))*nc;
        lim_pre_sq= lim_pre^2;

        
%         conditionsx_pre=Pos(1,:)<=liminf_pre|Pos(1,:)>=limsup_pre;
%         conditionsy_pre=Pos(2,:)<=liminf_pre|Pos(2,:)>=limsup_pre;
%         conditionsz_pre=Pos(3,:)<=liminf_pre|Pos(3,:)>=limsup_pre;
%         conditions_pre=conditionsx_pre|conditionsy_pre|conditionsz_pre;
        
        Pos(:,Pos(1,:).*Pos(1,:)+Pos(2,:).*Pos(2,:)+Pos(3,:).*Pos(3,:)>=lim_pre_sq)=[];
%     
%     [~,N_particle_in_part]=size(Pos);
%     
%     while (N_particle_in_part<N_particles_chunk_exact&&part_id<part)
%     
%     [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos_t  ] = preprocessing_part(root,spec,aux_path,filename,part,part_id);
%     
%     [ N_particle_in_part
%     
    
histogr_1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);
    
        
    parfor cor=1:num_cores
%         
        angl_ind_start=angl_indx(cor);
        angl_ind_end=angl_indx(cor+1)-1;
                
%         angles=angles_hpx_t(angl_ind_start:angl_ind_end);
        
        angles_t=angl_chunk_t{cor};
        angles_p=angl_chunk_p{cor};

%         rx=Pos;
        
%         rx(1,:)=Pos(1,:)-(nc/2)-pivot(1);
%         rx(2,:)=Pos(2,:)-(nc/2)-pivot(2);
%         rx(3,:)=Pos(3,:)-(nc/2)-pivot(3);
        
        histogr_1d_angles1=zeros(length(bins)-1,number_of_angle_nuple_hpx);

        
        for i=angl_ind_start:angl_ind_end
                    
%         theta=angles_hpx(1,i);
%         phi=angles_hpx(2,i);

        theta=angles_t(i-angl_ind_start+1);
        phi=angles_p(i-angl_ind_start+1);

        nx=[transpose(cos(theta).*cos(phi)) ,transpose(cos(theta).*sin(phi)), transpose(-sin(theta))];
        ny=[transpose(-sin(phi)) ,transpose(cos(phi)), 0];
        nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
        
        dz=nz*Pos;
        
        liminf=-(1/(2*lenght_factor))*nc;
        limsup= (1/(2*lenght_factor))*nc;        
        limsup_sq= 2*limsup^2;
        
        dz((nx*Pos).*(nx*Pos)+(ny*Pos).*(ny*Pos)>=limsup_sq|dz<=liminf|dz>=limsup)=[];
                
        histogr_1d_angles1(:,i)=histcounts(dz,bins);
%          histogr(i,:,:) = histcounts(dz,bins);
%         proj1d_angles(:,i)=proj1d_angles(:,i)+transpose(histogr(1,:));
%         
        end
        
        histogr_1d_angles=histogr_1d_angles+histogr_1d_angles1;
        
    end
    
    proj1d_angles=proj1d_angles+histogr_1d_angles;
    
%     part_id=part_id+1;
    
end
clearvars angl_chunk_t angl_chunk_p;


toc;




tic;
% 
%  angles(1,:)=angles_hpx(1,:);
%  angles(2,:)=angles_hpx(2,:);
%  number_of_angle_nuple_hpx = size(angles_hpx(1,:),2);

max_proj1d_angles=max(proj1d_angles);
average_proj1d_angles=mean(proj1d_angles,1);
max_amplitude_proj1d_angles=max_proj1d_angles(:)-average_proj1d_angles(:);

proj1d_index_max=zeros(1,number_of_angle_nuple_hpx);

parfor angl=1:number_of_angle_nuple_hpx
proj1d_index_max(1,angl)=find(proj1d_angles(:,angl)==max_proj1d_angles(1,angl),1);
end
proj1d_angles_snremoved=proj1d_angles;

for angl=1:number_of_angle_nuple_hpx
proj1d_angles_snremoved(proj1d_index_max(1,angl),angl)=average_proj1d_angles(angl);
end



std_proj1d_angles=std(proj1d_angles_snremoved);

clearvars proj1d_angles_snremoved;

stn_proj1d_angles=(max_amplitude_proj1d_angles(:))./std_proj1d_angles(:);


out_proj1d_angles=[transpose(max_amplitude_proj1d_angles);std_proj1d_angles;transpose(stn_proj1d_angles)];

% parfor angl=1:number_of_angle_nuple_hpx    %memory problem
for angl=1:number_of_angle_nuple_hpx
dc_proj1d_angles(:,angl)=(proj1d_angles(:,angl)-average_proj1d_angles(angl))./average_proj1d_angles(angl);
end

max_dc_proj1d_angles=max(dc_proj1d_angles);
max_amplitude_dc_proj1d_angles=max_dc_proj1d_angles(:);

dc_proj1d_index_max=zeros(1,number_of_angle_nuple_hpx);

% parfor angl=1:number_of_angle_nuple_hpx %memory problem
for angl=1:number_of_angle_nuple_hpx
dc_proj1d_index_max(1,angl)=find(dc_proj1d_angles(:,angl)==max_dc_proj1d_angles(1,angl),1);
end
dc_proj1d_angles_snremoved=dc_proj1d_angles;

for angl=1:number_of_angle_nuple_hpx
dc_proj1d_angles_snremoved(dc_proj1d_index_max(1,angl),angl)=0;
end

std_dc_proj1d_angles=std(dc_proj1d_angles_snremoved);
clearvars dc_proj1d_angles_snremoved;
stn_dc_proj1d_angles=(max_amplitude_dc_proj1d_angles(:))./std_dc_proj1d_angles(:);

out_dc_proj1d_angles=[transpose(max_amplitude_dc_proj1d_angles);std_dc_proj1d_angles;transpose(stn_dc_proj1d_angles)];

if level_window~=0        
    filtered_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);
    filtered_dc_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);
    
     level=floor(log2(length(bins)-1));

    
    parfor i=1:number_of_angle_nuple_hpx
%         average_proj1d_angles=mean2(proj1d_angles(:,i));
%         proj1d=(proj1d_angles(:,i)-average_proj1d_angles)/average_proj1d_angles;
        proj1d=proj1d_angles(:,i);
        [dwt_proj1d,levels] = wavedec(proj1d,level,dwbasis);

        dc_proj1d=dc_proj1d_angles(:,i);
        [dwt_dc_proj1d,levels] = wavedec(dc_proj1d,level,dwbasis);


        
%         [cwt_proj1d,periods] = cwt(proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
%         filtered_proj1d_angles(:,i) = icwt(cwt_proj1d,periods,[periods(low_pass) seconds(level_window)],'waveletparameters',[3 3.01]);
%         
%         dc_proj1d=dc_proj1d_angles(:,i);
%         [cwt_dc_proj1d,periods] = cwt(dc_proj1d,seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
        
        D=zeros(length(bins)-1,1);
        D_dc=zeros(length(bins)-1,1);
        for lev_win = 1:length(level_window)
            lvwin=level_window(lev_win);
            D(:,1)=D(:,1)+wrcoef('d',dwt_proj1d,levels,dwbasis,lvwin);
            D_dc(:,1)=D(:,1)+wrcoef('d',dwt_dc_proj1d,levels,dwbasis,lvwin);
        end
        filtered_proj1d_angles(:,i) = D(:,1);
        filtered_dc_proj1d_angles(:,i) = D_dc(:,1);
        
    end
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
    max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
    
    filtered_proj1d_index_max=zeros(1,number_of_angle_nuple_hpx);
    
    parfor angl=1:number_of_angle_nuple_hpx
        filtered_proj1d_index_max(1,angl)=find(filtered_proj1d_angles(:,angl)==max_filtered_proj1d_angles(1,angl),1);
    end
    filtered_proj1d_angles_snremoved=filtered_proj1d_angles;
    
    for angl=1:number_of_angle_nuple_hpx
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
    
    filtered_dc_proj1d_index_max=zeros(1,number_of_angle_nuple_hpx);
    
    parfor angl=1:number_of_angle_nuple_hpx
        filtered_dc_proj1d_index_max(1,angl)=find(filtered_dc_proj1d_angles(:,angl)==max_filtered_dc_proj1d_angles(1,angl),1);
    end
    filtered_dc_proj1d_angles_snremoved=filtered_dc_proj1d_angles;
    
    for angl=1:number_of_angle_nuple_hpx
        filtered_dc_proj1d_angles_snremoved(filtered_dc_proj1d_index_max(1,angl),angl)=0;
    end
        
    max_amplitude_filtered_dc_proj1d_angles=max_filtered_dc_proj1d_angles(:)-average_filtered_dc_proj1d_angles(:);
    std_filtered_dc_proj1d_angles=std(filtered_dc_proj1d_angles_snremoved);
    clearvars filtered_dc_proj1d_angles_snremoved;
    
    stn_filtered_dc_proj1d_angles=(max_amplitude_filtered_dc_proj1d_angles(:))./std_filtered_dc_proj1d_angles(:);
    
    out_filtered_dc_proj1d_angles=[transpose(max_amplitude_filtered_dc_proj1d_angles);std_filtered_dc_proj1d_angles;transpose(stn_filtered_dc_proj1d_angles)];
    
end

angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));

[~,number_of_angle_nuple_hpx] = size(angles_hpx);

proj1d_angles2=proj1d_angles;
dc_proj1d_angles2=dc_proj1d_angles;
filtered_dc_proj1d_angles2=filtered_dc_proj1d_angles;
filtered_proj1d_angles2=filtered_proj1d_angles;
out_proj1d_angles2=out_proj1d_angles;
out_dc_proj1d_angles2=out_dc_proj1d_angles;
out_filtered_proj1d_angles2=out_filtered_proj1d_angles;
out_filtered_dc_proj1d_angles2=out_filtered_dc_proj1d_angles;


% parfor angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx
% /memory problem
for angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx

    
    theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));

 theta_inverse_indices=number_of_angle_nuple_hpx-theta_indices+1;
 
%  phi_indice_inverse=theta_inverse_indices(find((mod(angles_hpx(2,angl)+pi,2*pi))==angles_hpx(2,theta_inverse_indices)));
 phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
 
    
    proj1d_angles(:,angl)=flip(proj1d_angles2(:,phi_indice_inverse));
    dc_proj1d_angles(:,angl)=flip(dc_proj1d_angles2(:,phi_indice_inverse));
    filtered_dc_proj1d_angles(:,angl)=flip(filtered_dc_proj1d_angles2(:,phi_indice_inverse));
    filtered_proj1d_angles(:,angl)=flip(filtered_proj1d_angles2(:,phi_indice_inverse));
    out_proj1d_angles(:,angl)=out_proj1d_angles2(:,phi_indice_inverse);
    out_dc_proj1d_angles(:,angl)=out_dc_proj1d_angles2(:,phi_indice_inverse);
    out_filtered_proj1d_angles(:,angl)=out_filtered_proj1d_angles2(:,phi_indice_inverse);
    out_filtered_dc_proj1d_angles(:,angl)=out_filtered_dc_proj1d_angles2(:,phi_indice_inverse);
    
end


if ~ismember(0,data_stream)
    
    path_out_all=strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/');
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/');
    mkdir(strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/dc/'));
    
    if level_window~=0 
        mkdir(path_out,strcat('/','level_window',mat2str(level_window(:)),'/'));
        mkdir(path_out,strcat('/dc','/','level_window',mat2str(level_window(:)),'/'));
%         mkdir(path_out,strcat('/cutoff_',num2str(cutoff),'MpcCut/dc/'));       
        mkdir(path_out_all,strcat('/','level_window',mat2str(level_window(:)),'/'));
        mkdir(path_out_all,strcat('/dc/','level_window',mat2str(level_window(:)),'/'));
        
    end
    
    if ismember(1,data_stream)
        
        fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        if level_window~=0             
            
            fileID = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
        end
        
        fileID = fopen(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
    end
    
    if ismember(2,data_stream)
            dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_dc_proj1d_angles,'delimiter','\t');
        
        if level_window~=0 
            dlmwrite(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),filtered_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),filtered_dc_proj1d_angles,'delimiter','\t');

        end
        
            dlmwrite(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'),dc_proj1d_angles,'delimiter','\t');

        
    end
    
end


cd('../wake_detection/box_statistics');
% 
toc;
delete(gcp('nocreate'))


end

