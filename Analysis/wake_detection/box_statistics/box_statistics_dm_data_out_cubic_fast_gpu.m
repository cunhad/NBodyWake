function [ proj1d_angles ] = box_statistics_dm_data_out_cubic_fast_gpu( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,data_stream,level_window,dwbasis)

%Computes for each pair of spherical angles, the projection of the
%particles positions on the corresponding axis, then construct a histogram
%of the resulting projection and filters it by extrating the first wavelet
%coefficient. The process is repeated for each different angle. 
%A paralelization setup is used in which all worker loads the same particle
%catalogue but each different work project to a different angle. If the
%angle plus particles does not fit in the worker, the particle catalogue is
%divided in "part" parts and at the end the histogram is added for each
%part.


%(example)  [proj1d_angles] = box_statistics_dm_data_out_cubic_fast_gpu('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','10.000xv0.dat',2,1,[0,0,0],2,1,4,[1,2],[1],'sym6');

% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% if specified, data will be stored in  root_out+spec+aux_path+aux_path_out

% "filename" is the output file from the nbody simulation

% lenght_factor = the analysis cube will have a lateral size given by the
% lateral size of the simulation cube divided by this number

% resol_factor= the bin of the histogram will the the particle bin size divided by this
%number

% pivot = a 3d array containing the translation wrt to the center of the
% cube (in grid cell units)

%NSIDE determines the number of angles probed, according to the healpix
%scheme: number_of_angles=12*NSIDE^2

% "part" specifies the partition of the particle catalogue done in order to
% it to fit in the memory of each worker.

%num_cores specifies the number of cores used

% data_stream set to 1, so it will generate binary outputs instead of usual
% txt outputs

%level_window is the wavelet filter coefficint considered. Usualy set to 1

%dwbasis is the discrete wavelet basis, usualy set to 'sym6', whici is the
%Symlets family. They are a modified version of Daubechies wavelets with increased symmetry.

%setup the parallel process specifications

myCluster = parcluster('local');
myCluster.NumWorkers=num_cores;
saveProfile(myCluster);

p = parpool(num_cores);

addpath(genpath('../../processing/'));


tic;
% 
cd('../../preprocessing');

[~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));


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
     angl_chunk_t{cr}=single(angles_hpx(1,angl_ind_start:angl_ind_end));
     angl_chunk_p{cr}=single(angles_hpx(2,angl_ind_start:angl_ind_end));   
end


[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,part,1);

%the variables below are stored as single, to better use of the gpus

if lenght_factor==1
    bins=single([1-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))-1]);
else
    bins=single([-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))]);
end

proj1d_angles=single(zeros(length(bins)-1,number_of_angle_nuple_hpx));

nc=single(nc);
pivot=single(pivot);
lenght_factor=single(lenght_factor);

for part_id = 1  :   part

    
    cd('../preprocessing');
    
    [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,part,part_id);
    
        Pos=single(mod(Pos,nc));
        
        
 
        Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
        Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
        Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
        
        
        lim_pre=single((1/(lenght_factor))*nc);
    
        Pos(:,abs(Pos(1,:))>lim_pre|abs(Pos(2,:))>lim_pre|abs(Pos(3,:))>lim_pre)=[];

    
histogr_1d_angles=single(zeros(length(bins)-1,number_of_angle_nuple_hpx));
    
        
    parfor cor=1:num_cores
        angl_ind_start=angl_indx(cor);
        angl_ind_end=angl_indx(cor+1)-1;
                
        %transfer data to the gpus
        
        angles_t=gpuArray(angl_chunk_t{cor});
        angles_p=gpuArray(angl_chunk_p{cor});
        
        pos=gpuArray(Pos);
        Bins=gpuArray(bins);

        
        histogr_1d_angles1=gpuArray(single(zeros(length(bins)-1,number_of_angle_nuple_hpx)));

        
        for i=angl_ind_start:angl_ind_end
                    

        theta=angles_t(i-angl_ind_start+1);
        phi=angles_p(i-angl_ind_start+1);

        nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
        
        dz=nz*pos;        
                
        histogr_1d_angles1(:,i)=histcounts(dz,Bins);
        
        end
        
        %gather data from the gpus
        
        histogr_1d_angles=histogr_1d_angles+gather(histogr_1d_angles1);
        
    end
    
    proj1d_angles=proj1d_angles+histogr_1d_angles;
    
    
end
clearvars angl_chunk_t angl_chunk_p;

toc;

tic;

[v f] = createCube; v = (v-[0.5 0.5 0.5])*nc ;


parfor i=1:number_of_angle_nuple_hpx        
    theta=angles_hpx(1,i);
    phi=angles_hpx(2,i);
    nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
    
    bin_sz=(length(bins)-1);
    area=zeros(bin_sz,1);
    for dist_1d=1:bin_sz
        plane = createPlane(nz*((bins(1,dist_1d+1)+bins(1,dist_1d))/2), nz);
        plane_overlap = intersectPlaneMesh(plane, v, f);
        area(dist_1d,1)=polygonArea3d(plane_overlap);
    end
    areas(:,i)=area(:,1);    
end

proj1d_angles=proj1d_angles./areas;

toc;


tic;

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

        proj1d=proj1d_angles(:,i);
        [dwt_proj1d,levels] = wavedec(proj1d,level,dwbasis);

        dc_proj1d=dc_proj1d_angles(:,i);
        [dwt_dc_proj1d,levels] = wavedec(dc_proj1d,level,dwbasis);

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

