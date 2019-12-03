function [  ] = claraCanRid_3dpart_data_anali_out_from_2d_clocangreal_atg_out( root,root_data_2d_in,root_data_2d_out,root_anali_2d_out,root_visual_2d,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,slices,lev_3d,lev_3drig,sigma,step_of_degree,wavel_removal_factor,snapshot,visual_type,visual_in_or_out,partition2d,partition3rd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 

% 
% 
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclar_l2lr1na256_data/';
% % % root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps256_1024_3dcurv_s5lv3_data/';
% root_data_2d_out='';
% root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclar-l2lr1na256_to_3dparcurv-l1lr1_anali/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclar-l2lr1na256_to_3dparcurv-l1lr1_visual/';
% % spec='4Mpc_2048c_1024p_zi63_nowakem';
% spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3003/';
% aux_path='/sample3007/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% pivot=[0,0,0];
% rot_angle=[1.5708,0,0];
% slices=32;
% % lev_2d=2;
% % lev_2drig=1;
% lev_3d=1;
% lev_3drig=1;
% sigma=5;        %not necessary
% step_of_degree=1*(180/256);
% wavel_removal_factor=1/2;
% % snapshot=[];
% snapshot=[23];
% % snapshot=[13,28]*(256/32);
% % snapshot=[9,29]*(128/32);
% visual_type=[1:2]; %if 1, shows the 2d proj; if 2 shows the ridgelet transformation
% % visual_in_or_out=[1,2,3]; %if 1 do visualization of the input, if 2 of
% % the output, if 3 does also canny output
%  visual_in_or_out=[2];
%  partition2d=1;
% partition3rd=2;
% %  sum_depth=1;
% %  sum_depth=4;
% 
% root='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/'
% root_data_2d_in='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps128_512_2dclara_l3lr2_data/'
% root_data_2d_out=''
% root_anali_2d_out='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps128_512_2dclara-l3lr2na1024_to_3dparclar_p4d2_rid-l1lr3_anali/'
% root_visual_2d='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps128_512_2dclara-l3lr2na1024_to_3dparclar_p4d2_rid-l1lr3_visual/'
% spec='4Mpc_2048c_1024p_zi63_nowakem'
% aux_path='/sample3001/'
% aux_path_out=''
% filename='3.000xv0.dat'
% lenght_factor=1
% resol_factor=0.5
% pivot=[0,0,0]
% rot_angle=[pi/2,0,0]
% slices=128
% lev_3d=1
% lev_3drig=3
% sigma=5
% step_of_degree=1*(180/1024)
% wavel_removal_factor=1/2
% snapshot=[]
% visual_type=[1:2]
% visual_in_or_out=[1:3]
% partition2d=4
% partition3rd=1

sum_depth=slices;

% filename_read_path=_2dproj_z3_data_sl;
% % 
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));
% 
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

run('../../../ridgelet/ppft3_nomex/ppft3/initpath.m')

% record the path to data in

% path_data_in=string(strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));


%making path to data out

if ~isempty(root_data_2d_out)
    
    path_data_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
    % string(strcat(root_data_2d_out,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    mkdir(char(strcat(root_data_2d_out,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
    
end


%making path to analisis out

path_anali_out=string(strcat(strcat(root_anali_2d_out,spec,aux_path),'anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
% string(strcat(root_data_2d_anali,spec,aux_path))
% string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
mkdir(char(strcat(root_anali_2d_out,spec,aux_path)),char(strcat('anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));


%making path to analisis depth out

if ~ismember(1,sum_depth)
    
    path_anali_depth_out=string(strcat(strcat(root_anali_2d_out,spec,aux_path),'anali_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
    % string(strcat(root_data_2d_anali,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    mkdir(char(strcat(root_anali_2d_out,spec,aux_path)),char(strcat('anali_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
    
end

%making path to visualization out

if ~isempty(snapshot)
    
    
    if ismember(1,visual_in_or_out)
        
        path_visual_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
        
        
    end
    
    if ismember(2,visual_in_or_out)
        
        path_visual=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
        
        
    end
    
    if ismember(3,visual_in_or_out)
        
        path_visual_canny=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_canny/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_canny/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
        
        
    end
    
    
    
    if (ismember(2,visual_type))
        
        %         if ismember(1,visual_in_or_out)
        %
        %             path_visual_rid_in=strcat(path_visual_in,'rid/');
        %             mkdir(char(path_visual_in),'rid/')
        %
        %         end
        
        if ismember(2,visual_in_or_out)
            
            path_visual_rid=strcat(path_visual,'rid/');
            mkdir(char(path_visual),'rid/')
            
        end
        
        if ismember(3,visual_in_or_out)
            
            path_visual_canny_rid=strcat(path_visual_canny,'rid/');
            mkdir(char(path_visual_canny),'rid/')
            
        end
        
    end
    
    
    
    if ~ismember(1,sum_depth)
        
        if ismember(1,visual_in_or_out)
            
            path_visual_depth_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
            
            
        end
        
        if ismember(2,visual_in_or_out)
            
            path_visual_depth=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
            
            
        end
        
        if ismember(3,visual_in_or_out)
            
            path_visual_canny_depth=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_canny_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_canny_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')));
            
            
        end
        
        
        
        if (ismember(2,visual_type))
            
            %         if ismember(1,visual_in_or_out)
            %
            %             path_visual_rid_in=strcat(path_visual_in,'rid/');
            %             mkdir(char(path_visual_in),'rid/')
            %
            %         end
            
            if ismember(2,visual_in_or_out)
                
                path_visual_rid_depth=strcat(path_visual_depth,'rid_depth/');
                mkdir(char(path_visual_depth),'rid_depth/')
                
            end
            
            if ismember(3,visual_in_or_out)
                
                path_visual_canny_rid_depth=strcat(path_visual_canny_depth,'rid_depth/');
                mkdir(char(path_visual_canny_depth),'rid_depth/')
                
            end
            
        end
        
        
    end
    
    
    
end



cd('../../../../preprocessing');


[~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

% [ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
% [~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );


% [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,nc,1);



% display(z)
z_string=char(filename);
z_string=z_string(1:end-7);
z=str2num(z_string);
z_glob=z;

filename_read=strcat('_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_data_sl');

% psi=rot_angle(3);
% theta=rot_angle(2);
% phi=rot_angle(1);

nb=np*resol_factor/lenght_factor;
nc_anal3d=nb/partition2d;


F=zeros(nb);
C_zero = fdct_wrapping(F,0);
 
% F2 = ones(nb,nb,slices/partition3rd);
% % F2 = ones(new_nc,new_nc,slices);
X2 = zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
X2(1+nc_anal3d/2,1+nc_anal3d/2,1+slices/(2*partition3rd))=(nc_anal3d*nc_anal3d*slices/partition3rd)^(1/3);
% X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
C2 = fdct3d_forward(X2);
E2 = cell(size(C2));
for s=1:length(C2)
    E2{s} = cell(size(C2{s}));
    for w=1:length(C2{s})
        A2 = C2{s}{w};
        E2{s}{w} = sqrt(sum(sum(sum(A2.*conj(A2)))) / prod(size(A2)));
    end
end
%
%
%
F2=zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
C_zero2 = fdct3d_forward(F2);

%radon

% nb_m=max(nb,slices);
nb_m=min(nb,slices);
% nb_m=sqrt(nb*slices);

Rad_norm = radon3(ones(nb_m,nb_m,nb_m));

map_3d_slices=zeros(nb,nb,slices);
map_3d_slices_filt3d=zeros(nb,nb,slices);
map_3d_slices_filtCanny3d=zeros(nb,nb,slices);


for slice_id=1:slices
    
    
    
    filename_read_path=strcat(strcat(root_data_2d_in,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/',filename_read,num2str(slice_id),'.bin')
    fid = fopen(filename_read_path);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    %             map = imresize(map,new_nc/nc,'triangle');
    
    
%     map(map<=1)=1;%to remove problem with holes
    
    
    map_3d_slices(:,:,slice_id)=map;
    
    
end


for partition=0:partition3rd-1
    for partition_x=0:partition2d-1
        for partition_y=0:partition2d-1
            
            C_aux=map_3d_slices((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
%       Caux=map_3d_slices_filt2d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);

            
            C = fdct3d_forward(C_aux);
            
            Ct = C;
            Ct2=C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero2{s}{w};
                    Ct2{s}{w} = C_zero2{s}{w};
                end
            end
            
            aux_count=1;
            for s =length(Ct)-lev_3d:length(Ct)-1
                %                 thresh=0;
                thresh = sigma + sigma*(s == length(C));
                
                
                As=0;
                Bs=0;
                Cs=0;
                Ds=0;
                
                size_nt=0;
                
                for part=0:5
                    
                    %do for each one of the 6 parts
                    
                    for w = 1+part*length(C{s})/6:(part+1)*length(C{s})/6
                        sz{w}=size(C{s}{w});
                    end
                    for  i=1:sz{1+part*length(C{s})/6}(1)
                        for j=1:sz{1+part*length(C{s})/6}(2)
                            for k=1:sz{1+part*length(C{s})/6}(3)
                                for w = 1+part*length(C{s})/6:(1+part)*length(C{s})/6
                                    i_c=ceil(i*sz{w}(1)/sz{1+part*length(C{s})/6}(1));
                                    j_c=ceil(j*sz{w}(2)/sz{1+part*length(C{s})/6}(2));
                                    k_c=ceil(k*sz{w}(3)/sz{1+part*length(C{s})/6}(3));
                                    a(w-part*length(C{s})/6)=abs(C{s}{w}(i_c,j_c,k_c))/E2{s}{w};
                                end
                                %                         std_a=std(a);
                                avr_a=mean(a);
                                %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                                %             treash_logic=double(a==max(a));
                                %             normalization=(a-mean(a)/std_a);
                                normalization=(a-avr_a)/avr_a;
                                for w = 1+part*length(C{s})/6:(part+1)*length(C{s})/6
                                    i_c=ceil(i*sz{w}(1)/sz{1+part*length(C{s})/6}(1));
                                    j_c=ceil(j*sz{w}(2)/sz{1+part*length(C{s})/6}(2));
                                    k_c=ceil(k*sz{w}(3)/sz{1+part*length(C{s})/6}(3));
                                    %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                                    Ct{s}{w}(i_c,j_c,k_c)=max(normalization(w-part*length(C{s})/6),Ct{s}{w}(i_c,j_c,k_c));
                                end
                            end
                        end
                    end
                end
                
                for w = 1:length(C{s})
                    
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                    
%                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E2{s}{w})
                    %                                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w};
                    %                   curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                    
                    
                    ave=mean(abs(Ct2{s}{w}(:))/E2{s}{w});
                    sig=std(abs(Ct2{s}{w}(:))/E2{s}{w})^2;
                    del=skewness(abs(Ct2{s}{w}(:))/E2{s}{w})*(sig^(3/2));
                    rho=kurtosis(abs(Ct2{s}{w}(:))/E2{s}{w})*(sig^(2));
                    
                    A_=ave;
                    B_=sig+A_^2;
                    C_=del+3*B_*A_-2*A_^3;
                    D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
                    
                    A_(isnan(A_))=0;
                    B_(isnan(B_))=0;
                    C_(isnan(C_))=0;
                    D_(isnan(D_))=0;
                    
                    size_n=prod(size(Ct2{s}{w}(:)));
                    size_nt=size_nt+size_n;
                    
                    As=A_*size_n+As;
                    Bs=B_*size_n+Bs;
                    Cs=C_*size_n+Cs;
                    Ds=D_*size_n+Ds;
                    
                    
                end
                %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                
                AT=As/size_nt;
                BT=Bs/size_nt;
                CT=Cs/size_nt;
                DT=Ds/size_nt;
                
                
                sigma_t=BT-AT^2;
                delta_t=CT-3*BT*AT+2*AT^3;
                rho_t=DT-4*CT*AT+6*BT*AT^2-3*AT^4;
                
                var_t=sigma_t;
                skew_t=delta_t/var_t^(3/2);
                kurt_t=rho_t/var_t^2;
                
                curv_1(aux_count)=AT;
                curv_2(aux_count)=var_t;
                curv_3(aux_count)=skew_t;
                curv_4(aux_count)=kurt_t;
                curv_5(aux_count)=rho_t;
                
                aux_count=aux_count+1;
            end
            
            %             map_3d_slices_filt2a3d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2) = real(fdct3d_inverse(Ct));
            map_3d_slices_filt3d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct2));
%           map_3d_slices_filt2a3d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct2));

            anali_curv(:,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,5)=curv_5(:);
            anali_curv(:,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,1)=curv_1(:);
            anali_curv(:,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,2)=curv_2(:);
            anali_curv(:,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,3)=curv_3(:);
            anali_curv(:,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,4)=curv_4(:);
            
        end
    end
end

anali3_curv1(:)=sum(anali_curv(:,:,1),2);
anali3_curv2(:)=sum(anali_curv(:,:,2),2);
anali3_curv3(:)=sum(anali_curv(:,:,3),2);
anali3_curv4(:)=sum(anali_curv(:,:,4),2);
anali3_curv5(:)=sum(anali_curv(:,:,5),2);

            anali3_curv(1,:)=anali3_curv1;
            anali3_curv(2,:)=anali3_curv2;
            anali3_curv(3,:)=anali3_curv3;
            anali3_curv(4,:)=anali3_curv4;
            anali3_curv(5,:)=anali3_curv5;

%             dc=(map-mean(map(:)))/mean(map(:));
%         dc=map;
%         dc(dc>cut)=cut;

%         nc_red=nc/red;
%         conv_=ones(red);
%         dc_red = conv2(dc,conv_,'valid');
%         dc = dc_red(1:red:end,1:red:end)/(red*red);
%

%
%
%             dc_cut=dc;
%             dc_cut(dc_cut>cut)=cut;
%             %         dc_cut(dc_cut>cut)=-1;
%             %         dc_cut = edge(dc_cut,'canny');
%             %         dc=double(dc_cut);
%
%             thresh = multithresh(dc_cut,trsh);
%             seg_I = imquantize(dc_cut,thresh);
%             %         test=dc_cut;
%             %         test=seg_I;
%             test=log(seg_I);
%     map(map<=1)=1;%to remove problem with holes
%     map_3d_slices(:,:,slice_id)=log(map);
%
%             map(map<=1)=1;%to remove problem with holes
%             map=map+1;
%             test=log(log(map));

%             if false
%     if ismember(slice_id,display_slice{find(sample_id_range==sample)})
%
%         figure; imagesc((size_mpc/new_nc)*[1:new_nc],(size_mpc/new_nc)*[1:new_nc],map_3d_slices(:,:,slice_id)); colorbar; axis('image');
%         xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%         ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%         set(gca,'FontName','FixedWidth');
%         set(gca,'FontSize',16);
%         set(gca,'linewidth',2);
%         title(strcat('log for sample ',num2str(sample),' slice ',num2str(slice_id)));
%     end


%3d curv filter

for slice_id=1:slices
    
    
    if ismember(slice_id,snapshot)&(ismember(1,visual_type))
        
        if ismember(1,visual_in_or_out)
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices(:,:,slice_id)); colorbar; axis('image');
            xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('filt2 info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_in','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                                close(fig);

        end
        
        if ismember(2,visual_in_or_out)
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt3d(:,:,slice_id)); colorbar; axis('image');
            xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('filt2 info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                                close(fig);

        end
    end
    
%     this=map_3d_slices_filt3d(:,:,slice_id);
    
%     anali(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    
    % theta = 0:180/nc:180;
    
    
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filt3d(:,:,slice_id),theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
    
    %     boudary_removal_factor=2048/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R_nor(1,:))
        %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
        %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
        %                 D(1:floor((448+200)/boudary_removal_factor))=0;
        D(length(xp)-nb*wavel_removal_factor:end)=0;
        D(1:nb*wavel_removal_factor)=0;
        
        R_nor_filt(:,i)=D;
    end
    %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
    %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
    R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
    R_nor_filt(1:nb*wavel_removal_factor,:)=[];
    
    
    if ismember(slice_id,snapshot)&(ismember(2,visual_type))&ismember(2,visual_in_or_out)
        
        %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
        fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
        xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
        ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
        set(gca,'FontName','FixedWidth');
        set(gca,'FontSize',10);
        set(gca,'linewidth',2);
        title(strcat('ridg filt2 for ',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
        
        saveas(fig,char(strcat(path_visual_rid','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                            close(fig);

    end
    
        this=map_3d_slices_filt3d(:,:,slice_id);
    
    anali(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    anali(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
%     anali3_curv(slice_id,1,:)=curv_1;
%     anali3_curv(slice_id,2,:)=curv_2;
%     anali3_curv(slice_id,3,:)=curv_3;
%     anali3_curv(slice_id,4,:)=curv_4;
%     anali3_curv(slice_id,5,:)=curv_5;
    
    %     path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    % string(strcat(root_data_2d_anali,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
    
    % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
    
    if ~isempty(root_data_2d_out)
        
        fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
        fwrite(fileID,map_3d_slices_filt3d(:,:,slice_id), 'float32','l');
        fclose(fileID);
        
    end
    
    
end


%now canny filter

for slice_id=1:slices
    
	[this,thresOut] = edge(map_3d_slices_filt3d(:,:,slice_id),'Canny',0.5);
    
    map_3d_slices_filtCanny3d(:,:,slice_id)=this;
    
    if ismember(slice_id,snapshot)&(ismember(1,visual_type))
                
        
        if ismember(3,visual_in_or_out)
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filtCanny3d(:,:,slice_id)); colorbar; axis('image');
            xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('filt3 plus canny, info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_canny','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                                close(fig);

        end
    end
    
%     this=map_3d_slices_filt3d(:,:,slice_id);
    
%     anali(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    
    % theta = 0:180/nc:180;
    
    
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filtCanny3d(:,:,slice_id),theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
    
    %     boudary_removal_factor=2048/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R_nor(1,:))
        %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
        %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
        %                 D(1:floor((448+200)/boudary_removal_factor))=0;
        D(length(xp)-nb*wavel_removal_factor:end)=0;
        D(1:nb*wavel_removal_factor)=0;
        
        R_nor_filt(:,i)=D;
    end
    %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
    %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
    R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
    R_nor_filt(1:nb*wavel_removal_factor,:)=[];
    
    
    if ismember(slice_id,snapshot)&(ismember(2,visual_type))&ismember(3,visual_in_or_out)
        
        %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
        fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
        xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
        ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
        set(gca,'FontName','FixedWidth');
        set(gca,'FontSize',10);
        set(gca,'linewidth',2);
        title(strcat('ridg filt3 and canny for ',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
        
        saveas(fig,char(strcat(path_visual_canny_rid,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                            close(fig);

    end
    
        this=map_3d_slices_filtCanny3d(:,:,slice_id);
    
    analiCanny(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    analiCanny(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    analiCanny(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    analiCanny(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
%     anali3_curv(slice_id,1,:)=curv_1;
%     anali3_curv(slice_id,2,:)=curv_2;
%     anali3_curv(slice_id,3,:)=curv_3;
%     anali3_curv(slice_id,4,:)=curv_4;
%     anali3_curv(slice_id,5,:)=curv_5;
    
    %     path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    % string(strcat(root_data_2d_anali,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
    
    % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
    
    if ~isempty(root_data_2d_out)
        
        fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
        fwrite(fileID,map_3d_slices_filtCanny3d(:,:,slice_id), 'float32','l');
        fclose(fileID);
        
    end
    
    
end

dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_analiCanny.txt'),analiCanny,'delimiter','\t');
dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali3_curv.txt'),anali3_curv,'delimiter','\t');




% now we do the depth analysis

if ~ismember(1,sum_depth)
    
    slices_depth=slices/sum_depth;
    map_3d_slices_depth=zeros(nb,nb,slices_depth);
    map_3d_slices_filt3d_depth=zeros(nb,nb,slices_depth);
    map_3d_slices_filt3d_canny_depth=zeros(nb,nb,slices_depth);


%         map_3d_slices_depth=map_3d_slices_depth+map_3d_slices;
%         map_3d_slices_filt3d_depth=map_2d_slices_filt2d_depth+map_3d_slices_filt2d;
%         slice_depth_id=floor(slice_id/sum_depth);

    for slice_depth_id=1:slices_depth
        
        map_3d_slices_depth(:,:,slice_depth_id)=sum(map_3d_slices(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        map_3d_slices_filt3d_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt3d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        map_3d_slices_filt3d_canny_depth(:,:,slice_depth_id)=sum(map_3d_slices_filtCanny3d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        
        if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(1,visual_type))
            
            if ismember(1,visual_in_or_out)
                
                fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_depth(:,:,slice_depth_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',10);
                set(gca,'linewidth',2);
                title(strcat('filt2 info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_depth_in','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                                    close(fig);

            end
            
            if ismember(2,visual_in_or_out)
                
                fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt3d_depth(:,:,slice_depth_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',10);
                set(gca,'linewidth',2);
                title(strcat('filt2 info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                                    close(fig);

            end
            
            if ismember(3,visual_in_or_out)
                
                fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt3d_canny_depth(:,:,slice_depth_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',10);
                set(gca,'linewidth',2);
                title(strcat('filt2 canny info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_canny_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                                    close(fig);

            end
            
        end
        %
        %         this=map_3d_slices_filt3d_depth(:,:,slice_depth_id);
        %
        %         anali_depth(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
        %
        % theta = 0:180/nc:180;
        
        
        
        %         theta = 0:step_of_degree:180;
        %         [R,xp] = radon(map_3d_slices_filt3d_depth(:,:,slice_depth_id),theta);
        
        R = radon3(imresize3(map_3d_slices_filt3d,[nb_m nb_m nb_m]));
        
        %         unit=ones(nb);
        %         [R_u,xp] = radon(unit,theta);
        
        %         frac_cut=0.5;
        %         R_nor=R;
        %         R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
        %         R_nor(R_u<=nb*frac_cut)=0;
        
        frac_cut=0.5;
        R_nor=R;
        R_nor(Rad_norm>nb_m*nb_m*frac_cut)=R_nor(Rad_norm>nb_m*nb_m*frac_cut)./Rad_norm(Rad_norm>nb_m*nb_m*frac_cut);
        R_nor(Rad_norm<=nb_m*nb_m*frac_cut)=0;
        
        %         %     boudary_removal_factor=2048/nb;
        %
        %         n_levels=floor(log2(length(R_nor(:,1))));
        %         R_nor_filt=zeros(size(R_nor));
        %         for i=1:length(R_nor(1,:))
        %             %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        %             [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        %             D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
        %             %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
        %             %                 D(1:floor((448+200)/boudary_removal_factor))=0;
        %             D(length(xp)-nb*wavel_removal_factor:end)=0;
        %             D(1:nb*wavel_removal_factor)=0;
        %
        %             R_nor_filt(:,i)=D;
        %         end
        %         %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
        %         %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
        %         R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
        %         R_nor_filt(1:nb*wavel_removal_factor,:)=[];
        
        WDEC = wavedec3(squeeze(R_nor(1,:,:,:)),1,'db1') ;
        R_nor_filt_x = waverec3(WDEC,'d',lev_3drig);
        WDEC = wavedec3(squeeze(R_nor(2,:,:,:)),1,'db1') ;
        R_nor_filt_y = waverec3(WDEC,'d',lev_3drig);
        WDEC = wavedec3(squeeze(R_nor(3,:,:,:)),1,'db1') ;
        R_nor_filt_z = waverec3(WDEC,'d',lev_3drig);
        
        R_nor_filt=[R_nor_filt_x;R_nor_filt_y;R_nor_filt_z];
        
        
        if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(2,visual_type))&ismember(2,visual_in_or_out)
            
            %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
            fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],sum(R_nor_filt_z,3));colorbar;
            xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('ridg filt2 for ',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_rid_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                                close(fig);

        end
        
         
        this=map_3d_slices_filt3d_depth(:,:,slice_depth_id);
        
        anali_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(this(:))];      
        anali_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(R(:))];
        anali_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(R_nor(:))];
        anali_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(R_nor_filt(:))];
        
        
        %now for Canny
        
        
%         
%         theta = 0:step_of_degree:180;
%         [R,xp] = radon(map_3d_slices_filt3d_canny_depth(:,:,slice_depth_id),theta);
        
        R = radon3(imresize3(map_3d_slices_filtCanny3d,[nb_m nb_m nb_m]));
        
%         unit=ones(nb);
%         [R_u,xp] = radon(unit,theta);
        
        
%         frac_cut=0.5;
%         R_nor=R;
%         R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
%         R_nor(R_u<=nb*frac_cut)=0;
        
        frac_cut=0.5;
        R_nor=R;
        R_nor(Rad_norm>nb_m*nb_m*frac_cut)=R_nor(Rad_norm>nb_m*nb_m*frac_cut)./Rad_norm(Rad_norm>nb_m*nb_m*frac_cut);
        R_nor(Rad_norm<=nb_m*nb_m*frac_cut)=0;
        
        %     boudary_removal_factor=2048/nb;
        
%         n_levels=floor(log2(length(R_nor(:,1))));
%         R_nor_filt=zeros(size(R_nor));
%         for i=1:length(R_nor(1,:))
%             %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%             [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%             D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
%             %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
%             %                 D(1:floor((448+200)/boudary_removal_factor))=0;
%             D(length(xp)-nb*wavel_removal_factor:end)=0;
%             D(1:nb*wavel_removal_factor)=0;
%             
%             R_nor_filt(:,i)=D;
%         end
%         %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
%         %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
%         R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
%         R_nor_filt(1:nb*wavel_removal_factor,:)=[];
        
        WDEC = wavedec3(squeeze(R_nor(1,:,:,:)),1,'db1') ;
        R_nor_filt_x = waverec3(WDEC,'d',lev_3drig);
        WDEC = wavedec3(squeeze(R_nor(2,:,:,:)),1,'db1') ;
        R_nor_filt_y = waverec3(WDEC,'d',lev_3drig);
        WDEC = wavedec3(squeeze(R_nor(3,:,:,:)),1,'db1') ;
        R_nor_filt_z = waverec3(WDEC,'d',lev_3drig);
        
        R_nor_filt=[R_nor_filt_x;R_nor_filt_y;R_nor_filt_z];


         
        
        if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(2,visual_type))&ismember(3,visual_in_or_out)
            
            %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
            fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],sum(R_nor_filt_z,3));colorbar;
            xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('ridg filt2 canny for ',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_canny_rid_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                                close(fig);

        end
        
         
        this=map_3d_slices_filt3d_canny_depth(:,:,slice_depth_id);
        
        analiCanny_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(this(:))];      
        analiCanny_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(R(:))];
        analiCanny_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(R_nor(:))];
        analiCanny_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(R_nor_filt(:))];
        
        
        
        %     path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
        % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
        
        % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
        
        
        
    end

    dlmwrite(strcat(path_anali_depth_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali3_depth',num2str(sum_depth),'.txt'),anali_depth,'delimiter','\t');
    dlmwrite(strcat(path_anali_depth_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_anali_depth',num2str(sum_depth),'.txt'),analiCanny_depth,'delimiter','\t');

end



cd('../wake_detection/curvelet/curvelab/3d/');

end

