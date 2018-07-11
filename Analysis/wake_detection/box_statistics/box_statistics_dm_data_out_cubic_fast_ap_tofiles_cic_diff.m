function [ proj1d_angles,filtered_proj1d_angles ] = box_statistics_dm_data_out_cubic_fast_ap_tofiles_cic( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,NSIDE,particl_part,angle_part,angle_p,angle_id,num_cores,data_stream,level_window,dwbasis)
%Computes for each pair of spherical angles, the projection of the
%particles positions on the corresponding axis, then construct a histogram
%of the resulting projection and filters it by extrating the first wavelet
%coefficient. The process is repeated for each different angle. 
%A paralelization setup is used in which all worker loads the same particle
%catalogue but each different work project to a different angle. If the
%angle plus particles does not fit in the worker, the particle catalogue is
%divided in "part" parts and at the end the histogram is added for each
%part.

%


%(example)  [proj1d_angles,filtered_proj1d_angles] = box_statistics_dm_data_out_cubic_fast_ap_tofiles_cic('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/box_stat_cubic_fast_cic/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','10.000xv0.dat',2,1,[0,0,0],2,1,4,4,1,1,[1,2],[1],'sym6');

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

% "part_p" specifies the partition of the particle catalogue done in order to
% it to fit in the memory of each worker.

% "angle_p" specifies the partition of the angles catalogue done in order to
% it to fit in the memory of each worker. This should be not samml than the number of workers. It is recommended that this
% quantity is an integer times the numbers of workers.


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

% load the redshift of the simulation

[~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

%load the angles to be used in the analysis and creates cells which
%contains the partition of those angles to be distributed in each node

n_angle_per_node=ceil(6*NSIDE*NSIDE/angle_p);

start_angl_indx=n_angle_per_node*(angle_id-1)+1;
end_angl_indx=n_angle_per_node*(angle_id);
if end_angl_indx>6*NSIDE*NSIDE
    end_angl_indx=6*NSIDE*NSIDE;
end

angles_hpx(1,:)  = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'),' ',[start_angl_indx 0 end_angl_indx 0]);
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'),' ',[start_angl_indx 0 end_angl_indx 0]);
angl_chuck_size=length(angles_hpx(1,:));
number_of_angle_nuple_hpx=angl_chuck_size;

% 
% angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
% angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));


% [~,number_of_angle_nuple_hpx] = size(angles_hpx);



% number_of_angle_nuple_hpx=number_of_angle_nuple_hpx/2; %we use just half of the angles since there is a reflection symmetry (because the projection in one axis is the same as the projection on the oposite axis)

n_angle_per_node=ceil(number_of_angle_nuple_hpx/angle_part);

for cr=1:angle_part+1
    angl_indx(cr)= n_angle_per_node*(cr-1)+1;
    if(angl_indx(cr)>number_of_angle_nuple_hpx)
        angl_indx(cr)=number_of_angle_nuple_hpx+1;
    end
end
for cr=1:angle_part
    n_angl_indx(cr)= -angl_indx(cr)+angl_indx(cr+1);
    
end

for cr=1:angle_part
     angl_ind_start=angl_indx(cr);
     angl_ind_end=angl_indx(cr+1)-1;
     angl_chunk_t{cr}=angles_hpx(1,angl_ind_start:angl_ind_end);
     angl_chunk_p{cr}=angles_hpx(2,angl_ind_start:angl_ind_end);   
end

clearvars angles_hpx;

%stores other information about the simulation (4 numbers only)

[ size_box, nc, np, ~, ~ ,~ ,~ ,~ ,z, ~, ~  ] = preprocessing_part(root,spec,aux_path,filename,particl_part,1);

%creates the bins to be used in the histograms (it has 2 fewer points if the volume analysed is maximum)

% if lenght_factor==1
%     bins=[1-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))-1];
% else
%     bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];
% end

bins=[-np*resol_factor/(2*lenght_factor):(np*resol_factor/(2*lenght_factor))];

% proj1d_angles(a,b) will store the 1d projection at point a for each angle
% a

proj1d_angles=zeros((np*resol_factor/lenght_factor),number_of_angle_nuple_hpx);

% this first loop only exist to avoid overload the memory of the workers,
% so the particle catalogue is partitioned, the projections are performed and the resulted is summed afterwards

% for angle_id = 1: angl_p

for part_id = 1  :   particl_part

    
    cd('../preprocessing');
    
    %loads part of the particles of the simulation
    
    [ ~, ~, ~, ~, ~ ,~ ,~ ,~ ,z, ~, Pos  ] = preprocessing_part(root,spec,aux_path,filename,particl_part,part_id);
    
    %remove particles that are outside the analysis box
    
        Pos=Pos*(np*resol_factor)/(nc);
    
    Pos=mod(Pos,np*resol_factor);
    
    
    Pos(1,:)=Pos(1,:)-(np*resol_factor/2)-pivot(1);
    Pos(2,:)=Pos(2,:)-(np*resol_factor/2)-pivot(2);
    Pos(3,:)=Pos(3,:)-(np*resol_factor/2)-pivot(3);
    
    
    lim_pre=(1/(lenght_factor))*np*resol_factor;
    
    Pos(:,abs(Pos(1,:))>lim_pre|abs(Pos(2,:))>lim_pre|abs(Pos(3,:))>lim_pre)=[];

    %this is an auxiliary variable that will be passed to the main variable
    %proj1d_angles
        
histogr_1d_angles=zeros((np*resol_factor/lenght_factor),number_of_angle_nuple_hpx);
    

ticBytes(gcp);

    %here the projections are performed by the workers
        %matlab says "Using sliced variables can reduce communication
        %between the client and workers. Only those slices needed by a
        %worker are sent to it when it starts working on a particular range
        %of indices.", but this seems not to be the case
    for cor=1:angle_part
        angl_ind_start=angl_indx(cor);
        angl_ind_end=angl_indx(cor+1)-1;
                
        %load the angles in the workers
        
        angles_t=angl_chunk_t{cor};
        angles_p=angl_chunk_p{cor};

           
        %created the auxiliary variable to histogr_1d_angles
        
        histogr_1d_angles1=zeros((np*resol_factor/lenght_factor),number_of_angle_nuple_hpx);

        %do the loop for each angle to perform the projections
        
        parfor i=angl_ind_start:angl_ind_end
                    
        %here are the angles

        theta=angles_t(i-angl_ind_start+1);
        phi=angles_p(i-angl_ind_start+1);
        
        %here is the unit vector associated with the angles above

        nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
        
        %here the projection is performed
        
        dz=nz*Pos;
        
        %and the histogram is computed
        
        %         histogr_1d_angles1(:,i)=histcounts(dz,bins);
        
        dz=dz+(np*resol_factor/(2*lenght_factor))+pivot(3);
        
        [histogr_1d_anglesaux]=cic(dz,(np*resol_factor/lenght_factor));
        
        histogr_1d_angles1(:,i)=histogr_1d_anglesaux;

        end
        
        %histogram is stored in the auxiliary variable
        
        histogr_1d_angles=histogr_1d_angles+histogr_1d_angles1;
        
    end
    
    %histogram is stored in the auxiliary variable
    
    proj1d_angles=proj1d_angles+histogr_1d_angles;
    
    tocBytes(gcp)
end

%we don't need to store the partition of angles anymore

clearvars angl_chunk_t angl_chunk_p;

angles_hpx(1,:)  = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'),' ',[start_angl_indx 0 end_angl_indx 0]);
angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'),' ',[start_angl_indx 0 end_angl_indx 0]);


toc;

tic;

%this part only computes a normalization factor, since diferent
%slices on the projection of the particles inside the cube have diferent
%areas. The normalization factor is this area.

[v f] = createCube; v = (v-[0.5 0.5 0.5])*2*np*resol_factor/lenght_factor ;


parfor i=1:number_of_angle_nuple_hpx        
%     display(i)
    theta=angles_hpx(1,i);
    phi=angles_hpx(2,i);
    nz=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
    
    bin_sz=(length(bins)-1);
    area=zeros(bin_sz,1);
    for dist_1d=1:bin_sz
        plane = createPlane(nz*((bins(1,dist_1d+1)+bins(1,dist_1d))/2), nz);
        plane_overlap = intersectPlaneMesh(plane, v+1e-10, f);
        area(dist_1d,1)=polygonArea3d(plane_overlap);
    end
    areas(:,i)=area(:,1);    
end

% %remove unanted boundaries
% 
% proj1d_angles(1:8,:)=0;
% proj1d_angles(end-7:end,:)=0;

%do the normalization

proj1d_angles=proj1d_angles./areas;

toc;


tic;

%in this part the peak of each 1d projection is taked, as well the standard
%deviation (with the peak removed)

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


%the same as above, but instead of the number of particles in each slice,
%the density contrast is used

% parfor angl=1:number_of_angle_nuple_hpx    %memory problem
for angl=1:number_of_angle_nuple_hpx
dc_proj1d_angles(:,angl)=(proj1d_angles(:,angl)-average_proj1d_angles(angl))./average_proj1d_angles(angl);
end

max_dc_proj1d_angles=max(dc_proj1d_angles);
max_amplitude_dc_proj1d_angles=max_dc_proj1d_angles(:);

dc_proj1d_index_max=zeros(1,number_of_angle_nuple_hpx);

%could not do parfor here. Matlab complains memory problem
% parfor angl=1:number_of_angle_nuple_hpx 


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

%here the wavelet filter takes place (on both fileds, the particle number and the density contrast)

if level_window~=0        
    filtered_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);
    filtered_dc_proj1d_angles=zeros(length(bins)-1,number_of_angle_nuple_hpx);
    
     level=floor(log2(length(bins)-1)); %level of the wavelet transformation

    
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
    
    %here we compute the absolute difference
    
    filtered_proj1d_angles=diff(filtered_proj1d_angles);
    filtered_proj1d_angles(end+1,:)=0;
    filtered_proj1d_angles=(abs(filtered_proj1d_angles));
    
    filtered_dc_proj1d_angles=diff(filtered_dc_proj1d_angles);
    filtered_dc_proj1d_angles(end+1,:)=0;
    filtered_dc_proj1d_angles=(abs(filtered_dc_proj1d_angles));
    
    
    %try this
    
    windowSize = 2;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    filtered_proj1d_angles = filter(b,a,filtered_proj1d_angles);
    
    windowSize = 2;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    filtered_dc_proj1d_angles = filter(b,a,filtered_dc_proj1d_angles);
    
    
    %avoid boundary problem
    
    filtered_proj1d_angles(1:8,:)=0;
    filtered_proj1d_angles(end-7:end,:)=0;
    filtered_dc_proj1d_angles(1:8,:)=0;
    filtered_dc_proj1d_angles(end-7:end,:)=0;
    
    %in this part the peak of each 1d projection is taked, as well the standard
    %deviation (with the peak removed). The diference is that now we are
    %performing these operations on the filtered 1d projections
    
    max_filtered_proj1d_angles=max(filtered_proj1d_angles);
    average_filtered_proj1d_angles=mean(filtered_proj1d_angles,1);
    %didn't removed the mean here
    %     max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:)-average_filtered_proj1d_angles(:);
    max_amplitude_filtered_proj1d_angles=max_filtered_proj1d_angles(:);
    
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
        
 %didn't removed the mean here  
%     max_amplitude_filtered_dc_proj1d_angles=max_filtered_dc_proj1d_angles(:)-average_filtered_dc_proj1d_angles(:);
    max_amplitude_filtered_dc_proj1d_angles=max_filtered_dc_proj1d_angles(:);
    std_filtered_dc_proj1d_angles=std(filtered_dc_proj1d_angles_snremoved);
    clearvars filtered_dc_proj1d_angles_snremoved;
    
    stn_filtered_dc_proj1d_angles=(max_amplitude_filtered_dc_proj1d_angles(:))./std_filtered_dc_proj1d_angles(:);
    
    out_filtered_dc_proj1d_angles=[transpose(max_amplitude_filtered_dc_proj1d_angles);std_filtered_dc_proj1d_angles;transpose(stn_filtered_dc_proj1d_angles)];
    
end


% 
% angles_hpx(1,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_t.cvs'));
% angles_hpx(2,:) = dlmread(strcat('../../python/angles',num2str(NSIDE),'_p.cvs'));
% 
% [~,number_of_angle_nuple_hpx] = size(angles_hpx);
% 
% proj1d_angles2=proj1d_angles;
% dc_proj1d_angles2=dc_proj1d_angles;
% filtered_dc_proj1d_angles2=filtered_dc_proj1d_angles;
% filtered_proj1d_angles2=filtered_proj1d_angles;
% out_proj1d_angles2=out_proj1d_angles;
% out_dc_proj1d_angles2=out_dc_proj1d_angles;
% out_filtered_proj1d_angles2=out_filtered_proj1d_angles;
% out_filtered_dc_proj1d_angles2=out_filtered_dc_proj1d_angles;
% 
% 
% %this part recover the result that the data for a given projection axis
% %(associated with a pair of angles) is the same as the data for the oposite
% %axis (oposite angles).
% 
% %could not do parfor here. Matlab complains memory problem
% for angl=number_of_angle_nuple_hpx/2+1:number_of_angle_nuple_hpx
% 
%     
%     theta_indices=find(angles_hpx(1,angl)==angles_hpx(1,:));
% 
%  theta_inverse_indices=number_of_angle_nuple_hpx-theta_indices+1;
%  
%  phi_indice_inverse=theta_inverse_indices(abs((mod(angles_hpx(2,angl)+pi,2*pi))-angles_hpx(2,theta_inverse_indices))<10E-10);
%  
%     
%     proj1d_angles(:,angl)=flip(proj1d_angles2(:,phi_indice_inverse));
%     dc_proj1d_angles(:,angl)=flip(dc_proj1d_angles2(:,phi_indice_inverse));
%     filtered_dc_proj1d_angles(:,angl)=flip(filtered_dc_proj1d_angles2(:,phi_indice_inverse));
%     filtered_proj1d_angles(:,angl)=flip(filtered_proj1d_angles2(:,phi_indice_inverse));
%     out_proj1d_angles(:,angl)=out_proj1d_angles2(:,phi_indice_inverse);
%     out_dc_proj1d_angles(:,angl)=out_dc_proj1d_angles2(:,phi_indice_inverse);
%     out_filtered_proj1d_angles(:,angl)=out_filtered_proj1d_angles2(:,phi_indice_inverse);
%     out_filtered_dc_proj1d_angles(:,angl)=out_filtered_dc_proj1d_angles2(:,phi_indice_inverse);
%     
% end
% 
%now the data is saved

if ~ismember(0,data_stream)
    
    path_out_all=strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/');
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/');
    mkdir(strcat(strcat(root_out(1,1:end-1),'_all/',spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/parts/'));
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/',dwbasis,'/dc/parts/'));
    
    if level_window~=0 
        mkdir(path_out,strcat('/','level_window',mat2str(level_window(:)),'/'));
        mkdir(path_out,strcat('/dc','/','level_window',mat2str(level_window(:)),'/'));
        mkdir(path_out_all,strcat('/','level_window',mat2str(level_window(:)),'/'));
        mkdir(path_out_all,strcat('/dc/','level_window',mat2str(level_window(:)),'/'));
        
    end
    
    if ismember(1,data_stream)
        
        fileID = fopen(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,out_dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
        if level_window~=0             
            
            fileID = fopen(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,out_filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,filtered_proj1d_angles, 'float32','l');
            fclose(fileID);
            
            fileID = fopen(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
            fwrite(fileID,filtered_dc_proj1d_angles, 'float32','l');
            fclose(fileID);
            
        end
        
        fileID = fopen(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,proj1d_angles, 'float32','l');
        fclose(fileID);
        
        fileID = fopen(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z','_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.bin'),'w');
        fwrite(fileID,dc_proj1d_angles, 'float32','l');
        fclose(fileID);
        
    end
    
    if ismember(2,data_stream)
            dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z)),'_out_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_dc_proj1d_angles,'delimiter','\t');
        
        if level_window~=0 
            dlmwrite(strcat(path_out,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_out_filtered_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),out_filtered_dc_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),filtered_proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'dc/','level_window',mat2str(level_window(:)),'/','_',num2str(find(str2num(char(redshift_list))==z)),'_filtered_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),filtered_dc_proj1d_angles,'delimiter','\t');

        end
        
            dlmwrite(strcat(path_out_all,'_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
            dlmwrite(strcat(path_out_all,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_dc_1dproj_angle_z',num2str(z),'_anglparts',num2str(angle_p),'_anglid',num2str(angle_id),'_NSIDE',num2str(NSIDE),'.txt'),dc_proj1d_angles,'delimiter','\t');

        
    end
    
end

cd('../wake_detection/box_statistics');
% 
toc;
delete(gcp('nocreate'))


end


function [histogr_1d_angles1]=cic(dz,size)

histogr_1d_angles1=zeros(1,size);
    
    for particle=1:length(dz)
        
       x=dz(1,particle)-0.5;
       i1=floor(x)+1;
       i2=i1+1;
       dx1=i1-x;
       dx2=1-dx1;
       
        if (i1>0 && i2<= size) 
         histogr_1d_angles1(i1)=histogr_1d_angles1(i1)+dx1;
         histogr_1d_angles1(i2)=histogr_1d_angles1(i2)+dx2;
        end
        
    end


    
end

