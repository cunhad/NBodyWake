function [ map,anali ] = slices_curv_2a3d_locangrealaRid_inf( root,root_data_2d_in,root_data_2d_out,root_data_2d_anali_out,root_visual_2d,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,numb_rand,slices,lev,lev_2drid,lev_3d,lev_3drid,step_of_degree,wavel_removal_factor,NSIDE,partition2d,partition3rd,sum_depth,snapshot,visual_type,visual_in_or_out,stage,partition2d_dept,partition3rd_dept)


% Ridgelet done in the deep analysis after 3d filter


% (example) [ map ] = slices_curv_2d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,8,4,8,2,5 );

%(example) [ map ] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/ht/','/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/','/home/asus/Dropbox/extras/storage/graham/ht/analy_cps32_1024/','4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m','/sample3001/','','3.000xv0.dat',1,1,1,32,8,2,5 );
%(example) for i=1:24; [ count_sum] = proj2d_cic_dm_data_out_part_rand_hpx_slices('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,1,1,i,2,2); end;
%(example) for i=1:24; [ map ,anali] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,i,2,2,2,5 ); end;
%  
% % 
% root='/home/asus/Dropbox/extras/storage/graham/small_res/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/small_res/data_test3/';
% root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/small_res/data3/';
% root_data_2d_anali_out='/home/asus/Dropbox/extras/storage/graham/small_res/anali3/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/visual3';
% spec='64Mpc_256c_128p_zi63_nowakem';
% aux_path='/sample2001/';
% % spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3007/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='10.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% numb_rand=1;
% slices=32;
% lev=2;
% lev_2drid=2;
% lev_3d=1;
% lev_3drid=2;
% step_of_degree=(64/180);
% wavel_removal_factor=1/2;
% NSIDE=8;
% partition2d=1;
% partition3rd=2;
% sum_depth=4;
% % snapshot=[];
% snapshot=[23];
% % snapshot=[13,28]*(256/32);
% % snapshot=[9,29]*(128/32);
% visual_type=[1,2,3,4,5];        %if 1, shows the 2d proj; %if 2 shows the ridgelet transformation    %if 3 does everything above also for the depth                          
% visual_in_or_out=[1:3]; %if 1 do visualization of the input,%  if 2 of the 2d_filtered    %  if 3 of the 3d filtered                            
% %  visual_in_or_out=[2];
%  stage=[3]              % for doing the 3d
% % %  info=[0,1,2]             %0,for just simplest visualisation
% %                             %1 with minimal information
% %                             %2 wthi some information
% %                             %3for maximal information
%   


% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_hpx_2d/';
% root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_tst1/';
% root_data_2d_anali_out='/home/asus/Dropbox/extras/storage/graham/ht/data_tst_anali1/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_tst_anali_visual1/';
% spec='4Mpc_2048c_1024p_zi63_nowakem';
% aux_path='/sample3004/';
% % spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3004/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% numb_rand=95;
% slices=32;
% lev=2;
% lev_2drid=1;
% lev_3d=1;
% lev_3drid=1;
% step_of_degree=(180/256);
% wavel_removal_factor=1/2;
% NSIDE=8;
% partition2d=1;
% partition3rd=2;
% sum_depth=4;
% % snapshot=[];
% snapshot=[23];
% % snapshot=[13,28]*(256/32);
% % snapshot=[9,29]*(128/32);
% visual_type=[1:5];        %if 1, shows the 2d proj; %if 2 shows the ridgelet transformation    %if 3 does everything above also for the depth % 
                            %if 4 produce figures to be used in nn analysis
                            %if 5 produce figures to be used in nn analysis
% visual_in_or_out=[1:3]; %if 1 do visualization of the input,%  if 2 of the 2d_filtered    %  if 3 of the 3d filtered                            
% %  visual_in_or_out=[2];
% stage=[3];              % for doing the 3d
% %  info=[0,1,2]             %0,for just simplest visualisation
%                             %1 with minimal information
%                             %2 with some information
%                             %3 for maximal information
                           

%%%%Problems

%%%%    anali3_curv definition

% 
% 
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

% 
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

 addpath(genpath('../ridgelet/'));

cd('../../preprocessing');

[~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,1024,1);


z_glob=str2num(filename(1:end-7))


nb=np*resol_factor/lenght_factor;
nc_anal3d=nb/partition2d;


% filename='_2dproj_z3_data_sl';    

% stage=[2,3]
%     if ismember(1,visual_in_or_out)

F = ones(nb);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%X = F * sqrt(prod(size(F)));
%C = fdct_wrapping(X,0);
C = fdct_wrapping(X,0,2);
%C = fdct_wrapping(F,0);
E = cell(size(C));
for s=1:length(C)
    E{s} = cell(size(C{s}));
    for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
    end
end

F=zeros(nb);
C_zero = fdct_wrapping(F,0);

% F2 = ones(nb,nb,slices/partition3rd);
% % F2 = ones(new_nc,new_nc,slices);
% X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));

X2 = zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
X2(1+nc_anal3d/2,1+nc_anal3d/2,1+slices/(2*partition3rd))=(nc_anal3d*nc_anal3d*slices/partition3rd)^(1/3);

if ismember(3,stage)
    
    C2 = fdct3d_forward(X2);
    E2 = cell(size(C2));
    for s=1:length(C2)
        E2{s} = cell(size(C2{s}));
        for w=1:length(C2{s})
            A2 = C2{s}{w};
            E2{s}{w} = sqrt(sum(sum(sum(A2.*conj(A2)))) / prod(size(A2)));
        end
    end
    
    F2=zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
    C_zero2 = fdct3d_forward(F2);
    
end

anali=zeros(slices,4,5);
anali_curv=zeros(slices,5,lev);
% anali3_curv=zeros(5,(lev_3d)*partition2d*partition2d*partition3rd);
anali2a3=zeros(slices,4,5);
anali2a3_curv=zeros(slices,5,lev);

if ~ismember(1,sum_depth)
    anali_depth=zeros(slices/sum_depth,4,5);
    anali2a3_depth=zeros(partition2d_dept*partition2d_dept*partition3rd_dept,4);
    anali2a3_curv_depth=zeros(slices/sum_depth,5,lev);
end

mkdir(root_data_2d_anali_out);
mkdir(root_data_2d_anali_out,strcat(spec,aux_path));


path1=strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/');

path2=dir(char(strcat(path1,"*pv*")));
path2={path2.name};
path2=path2(1);

path3='/2dproj/dm/';

%making path to analisis out

path_data_2d_anali_out=string(strcat(strcat(root_data_2d_anali_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
mkdir(char(strcat(root_data_2d_anali_out,spec,aux_path)),char(strcat('data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));

%making path to data out

if ~isempty(root_data_2d_out)
    
    mkdir(root_data_2d_out);
    mkdir(root_data_2d_out,strcat(spec,aux_path));
        
    path_data_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    mkdir(char(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    if ismember(3,stage)
        path_data3_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt3_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        mkdir(char(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt3_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    end
end


%making path to visualization out

if ~isempty(snapshot)
    
    
    if ismember(1,visual_in_or_out)
        
        path_visual_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
        
        
    end
    
    if ismember(2,visual_in_or_out)
        
        path_visual_2dfilt=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));                
        
        if (ismember(4,visual_type)|ismember(5,visual_type))
            
            path_visual_2dfilt_fig=string(strcat(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path),'visual_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path)),char(strcat('visual_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
            
        end
        
    end
    
    if ismember(3,visual_in_or_out)
        
        path_visual_3dfilt=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));                
        
        if (ismember(4,visual_type)|ismember(5,visual_type))
            
            path_visual_3dfilt_fig=string(strcat(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path),'visual_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path)),char(strcat('visual_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
            
        end
        
    end
    
    
    if (ismember(2,visual_type))
        
        %         if ismember(1,visual_in_or_out)
        %
        %             path_visual_rid_in=strcat(path_visual_in,'rid/');
        %             mkdir(char(path_visual_in),'rid/')
        %
        %         end
        
        if ismember(2,visual_in_or_out)
            
            path_visual_rid_2dfilt=strcat(path_visual_2dfilt,'rid/');
            mkdir(char(path_visual_2dfilt),'rid/')
            
        end
        
        if ismember(3,visual_in_or_out)
            
            path_visual_rid_3dfilt=strcat(path_visual_3dfilt,'rid/');
            mkdir(char(path_visual_3dfilt),'rid/')
            
        end
        
    end
    
    
    
    if ~ismember(1,sum_depth)
        
        if ismember(1,visual_in_or_out)
            
            path_visual_depth_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
            
            
        end
        
        if ismember(2,visual_in_or_out)
            
            path_visual_depth_2dfilt=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
            
            if (ismember(4,visual_type)|ismember(5,visual_type))
                
                path_visual_depth_2dfilt_fig=string(strcat(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path),'visual_depth_',num2str(sum_depth),'_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
                % string(strcat(root_data_2d_anali,spec,aux_path))
                % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
                mkdir(char(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_2dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
                
            end
            
        end
        
        if ismember(3,visual_in_or_out)
            
            path_visual_depth_3dfilt=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
            
            if (ismember(4,visual_type)|ismember(5,visual_type))
                
                path_visual_depth_3dfilt_fig=string(strcat(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path),'visual_depth_',num2str(sum_depth),'_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
                % string(strcat(root_data_2d_anali,spec,aux_path))
                % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
                mkdir(char(strcat(root_visual_2d(1:end-1),'_fig/',spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_3dfilt/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
                
            end
            
        end
        
        
        
        
        if (ismember(2,visual_type))
            
            %         if ismember(1,visual_in_or_out)
            %
            %             path_visual_rid_in=strcat(path_visual_in,'rid/');
            %             mkdir(char(path_visual_in),'rid/')
            %
            %         end
            
            if ismember(2,visual_in_or_out)
                
                path_visual_rid_depth_2dfilt=strcat(path_visual_depth_2dfilt,'rid_depth/');
                mkdir(char(path_visual_depth_2dfilt),'rid_depth/')
                
            end
            
            if ismember(3,visual_in_or_out)
                
                path_visual_rid_depth_3dfilt=strcat(path_visual_depth_3dfilt,'rid_depth/');
                mkdir(char(path_visual_depth_3dfilt),'rid_depth/')
                
            end
        end
        
        
    end
    
    
    
end


map_3d_slices_pre=zeros(nb,nb,slices);
map_3d_slices=zeros(nb,nb,slices);
map_3d_slices_filt2d=zeros(nb,nb,slices);
map_3d_slices_filt3d=zeros(nb,nb,slices);

% for slice_id=1:1
for slice_id=1:slices

    
    filename_2d=string(strcat(path1,path2,path3,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'))
    fid = fopen(filename_2d);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    
%     map2=map;
%     map(map<=1)=1;%to remove problem with holes
%     map_3d_slices(:,:,slice_id)=log(map);
map_3d_slices_pre(:,:,slice_id)=map;

end

for slice_id=1:slices
    
    map=map_3d_slices_pre(:,:,slice_id);
    dc=(map-mean(map_3d_slices_pre(:)))/mean(map_3d_slices_pre(:));
    %             dc=(map-mean(map(:)))/mean(map(:));
    %             map_3d_slices(:,:,slice_id)=(atan(dc+1));
    map_3d_slices(:,:,slice_id)=(atan((dc+1)*16));

    if ismember(slice_id,snapshot)&(ismember(1,visual_type))&ismember(1,visual_in_or_out)
        
        
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
    
    C = fdct_wrapping(map_3d_slices(:,:,slice_id),0);
    Ct = C;
    Ct2=C;    
    for s = 1:length(C)
        for w = 1:length(C{s})
            Ct{s}{w} = C_zero{s}{w};
            Ct2{s}{w} = C_zero{s}{w};
        end
    end
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        
        %right part
                sz=cell(length(C{s}),1);
                for w = 1:length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1}(1)
                    for j=1:sz{1}(2)
                        for w = 1:length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1}(2));
                            a(w)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a);
%                         treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1:length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w)*(C{s}{w}(i_c,j_c))*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1:length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                
                %right2 part
                for w = 1+length(C{s})/4:2*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+length(C{s})/4}(1)
                    for j=1:sz{1+length(C{s})/4}(2)
                        for w = 1+length(C{s})/4:2*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+length(C{s})/4}(2));
                            a(w-length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+length(C{s})/4:2*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-length(C{s})/4),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1+length(C{s})/4:2*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                %right3 part
                for w = 1+2*length(C{s})/4:3*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+2*length(C{s})/4}(1)
                    for j=1:sz{1+2*length(C{s})/4}(2)
                        for w = 1+2*length(C{s})/4:3*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+2*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+2*length(C{s})/4}(2));
                            a(w-2*length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+2*length(C{s})/4:3*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+2*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+2*length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-2*length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-2*length(C{s})/4),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1+2*length(C{s})/4:3*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                %right4 part
                for w = 1+3*length(C{s})/4:4*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+3*length(C{s})/4}(1)
                    for j=1:sz{1+3*length(C{s})/4}(2)
                        for w = 1+3*length(C{s})/4:4*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+3*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+3*length(C{s})/4}(2));
                            a(w-3*length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+3*length(C{s})/4:4*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+3*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+3*length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-3*length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-3*length(C{s})/4),Ct{s}{w}(i_c,j_c));
                            
                        end
                    end
                end
                
                for w = 1+3*length(C{s})/4:4*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                a=[];
                treash_logic=[];
                normalization=[];
                
        
        
        for w = 1:length(C{s})
            %                     Ct{s}{w} = C{s}{w};
            %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
            %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
%             Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
            %                 Ct{s}{w} = C{s}{w};
            %              curv(slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
            
%             ave=mean(abs(Ct{s}{w}(:)));
%             sig=std(abs(Ct{s}{w}(:)))^2;
%             del=skewness(abs(Ct{s}{w}(:)))*(sig^(3/2));
%             rho=kurtosis(abs(Ct{s}{w}(:)))*(sig^(2));

            ave=mean(abs(Ct2{s}{w}(:))/E{s}{w});
            sig=std(abs(Ct2{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(Ct2{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(Ct2{s}{w}(:))/E{s}{w})*(sig^(2));

            A_=ave;
            B_=sig+A_^2;
            C_=del+3*B_*A_-2*A_^3;
            D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
            
            size_n=prod(size(Ct2{s}{w}(:)));
            size_nt=size_nt+size_n;
            
            As=A_*size_n+As;
            Bs=B_*size_n+Bs;
            Cs=C_*size_n+Cs;
            Ds=D_*size_n+Ds;
            
            
            
        end
        %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
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
    
    
    
    
    
    
    
    
    map_3d_slices_filt2d(:,:,slice_id) = real(ifdct_wrapping(Ct2,0)); 
    
    
    if ismember(slice_id,snapshot)&(ismember(1,visual_type))&ismember(2,visual_in_or_out)
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt2d(:,:,slice_id)); colorbar; axis('image');
            xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('filt2d info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_2dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                        close(fig);
                        
    end
    
    %boundary removal
    
    map_3d_slices_filt2d(1:16,:,slice_id)=0;
    map_3d_slices_filt2d(end-16:end,:,slice_id)=0;
    map_3d_slices_filt2d(:,1:16,slice_id)=0;
    map_3d_slices_filt2d(:,end-16:end,slice_id)=0;
    
    
    if ~isempty(root_data_2d_out)
    
    fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
    fwrite(fileID,map_3d_slices_filt2d(:,:,slice_id), 'float32','l');
    fclose(fileID);
    end
    
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filt2d(:,:,slice_id),theta);
    
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
        D = wrcoef('d',dc_dwt,levels,'db1',lev_2drid);
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
        
        saveas(fig,char(strcat(path_visual_rid_2dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                    close(fig);

    end
    
    
    this=map_3d_slices_filt2d(:,:,slice_id);
    anali(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    anali(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
    anali_curv(slice_id,1,:)=curv_1;
    anali_curv(slice_id,2,:)=curv_2;
    anali_curv(slice_id,3,:)=curv_3;
    anali_curv(slice_id,4,:)=curv_4;
    anali_curv(slice_id,5,:)=curv_5;
    
end

dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_kurt.txt'),anali_curv,'delimiter','\t');



% now we do the depth analysis of the in and 2dfilter types

if ~ismember(1,sum_depth)
    
    slices_depth=slices/sum_depth;
    map_3d_slices_depth=zeros(nb,nb,slices_depth);
    map_3d_slices_filt2d_depth=zeros(nb,nb,slices_depth);

%         map_3d_slices_depth=map_3d_slices_depth+map_3d_slices;
%         map_3d_slices_filt3d_depth=map_2d_slices_filt2d_depth+map_3d_slices_filt2d;
%         slice_depth_id=floor(slice_id/sum_depth);

    for slice_depth_id=1:slices_depth
        
        map_3d_slices_depth(:,:,slice_depth_id)=sum(map_3d_slices(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        map_3d_slices_filt2d_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt2d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3); 
        
        if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(3,visual_type))
            
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
                
                fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt2d_depth(:,:,slice_depth_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',10);
                set(gca,'linewidth',2);
                title(strcat('filt2 info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_depth_2dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                            close(fig);

            end
        end
        
        C = fdct_wrapping(map_3d_slices_filt2d_depth(:,:,slice_depth_id),0);
%     Ct = C;
%     Ct2=C;    
   
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        for w = 1:length(C{s})

            ave=mean(abs(C{s}{w}(:))/E{s}{w});
            sig=std(abs(C{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(C{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(C{s}{w}(:))/E{s}{w})*(sig^(2));

            A_=ave;
            B_=sig+A_^2;
            C_=del+3*B_*A_-2*A_^3;
            D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
            
            size_n=prod(size(C{s}{w}(:)));
            size_nt=size_nt+size_n;
            
            As=A_*size_n+As;
            Bs=B_*size_n+Bs;
            Cs=C_*size_n+Cs;
            Ds=D_*size_n+Ds;
            
            
            
        end
        %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
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
        
        curv2_1_depth(aux_count)=AT;
        curv2_2_depth(aux_count)=var_t;
        curv2_3_depth(aux_count)=skew_t;
        curv2_4_depth(aux_count)=kurt_t;
        curv2_5_depth(aux_count)=rho_t;
        
        aux_count=aux_count+1;
    end
        
        theta = 0:step_of_degree:180;
        [R,xp] = radon(map_3d_slices_filt2d_depth(:,:,slice_depth_id),theta);
        
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
            D = wrcoef('d',dc_dwt,levels,'db1',lev_2drid);
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
                      
        if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(2,visual_type))&ismember(2,visual_in_or_out)
            
            %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
            fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
            xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('ridg filt2d for ',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_rid_depth_2dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                        close(fig);

        end
         
        this=map_3d_slices_filt2d_depth(:,:,slice_depth_id);
        
        anali_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];      
        anali_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
        anali_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
        anali_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];                            
        
        anali_curv_depth(slice_id,1,:)=curv2_1_depth;
        anali_curv_depth(slice_id,2,:)=curv2_2_depth;
        anali_curv_depth(slice_id,3,:)=curv2_3_depth;
        anali_curv_depth(slice_id,4,:)=curv2_4_depth;
        anali_curv_depth(slice_id,5,:)=curv2_5_depth;
    end

    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_depth',num2str(sum_depth),'.txt'),anali_depth,'delimiter','\t');
    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_kurt_depth',num2str(sum_depth),'.txt'),anali_curv_depth,'delimiter','\t');
    
end

if ismember(3,stage)
    
    
    for partition=0:partition3rd-1
        for partition_x=0:partition2d-1
            for partition_y=0:partition2d-1
                
                C_aux=map_3d_slices_filt2d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
                
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
                    %                 thresh = sigma + sigma*(s == length(C));
                    
                    
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
%                 map_3d_slices_filt3d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct));
                map_3d_slices_filt3d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct2));

                anali3_curv_aux(1,partition+1,partition_x+1,partition_y+1,:)=curv_1(:);
                anali3_curv_aux(2,partition+1,partition_x+1,partition_y+1,:)=curv_2(:);
                anali3_curv_aux(3,partition+1,partition_x+1,partition_y+1,:)=curv_3(:);
                anali3_curv_aux(4,partition+1,partition_x+1,partition_y+1,:)=curv_4(:);
                anali3_curv_aux(5,partition+1,partition_x+1,partition_y+1,:)=curv_5(:);
                
            end
        end
    end
    
    anali3_curv1_=anali3_curv_aux(1,:,:,:,:);
    anali3_curv2_=anali3_curv_aux(2,:,:,:,:);
    anali3_curv3_=anali3_curv_aux(3,:,:,:,:);
    anali3_curv4_=anali3_curv_aux(4,:,:,:,:);
    anali3_curv5_=anali3_curv_aux(5,:,:,:,:);
    
    
    anali3_curv1=anali3_curv1_(:);
    anali3_curv2=anali3_curv2_(:);
    anali3_curv3=anali3_curv3_(:);
    anali3_curv4=anali3_curv4_(:);
    anali3_curv5=anali3_curv5_(:);
    
    
    
    anali3_curv(1,:)=anali3_curv1;
    anali3_curv(2,:)=anali3_curv2;
    anali3_curv(3,:)=anali3_curv3;
    anali3_curv(4,:)=anali3_curv4;
    anali3_curv(5,:)=anali3_curv5;
    
    
    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali3_kurt.txt'),anali3_curv,'delimiter','\t');
    
    
    %now the analysis of the 3d filtered slices
    
    for slice_id=1:slices
        
        if ismember(slice_id,snapshot)&(ismember(1,visual_type))
            
            
            if ismember(3,visual_in_or_out)
                
                fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt3d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',10);
                set(gca,'linewidth',2);
                title(strcat('filt2 info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_3dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                            close(fig);

            end
        end
        
        C = fdct_wrapping(map_3d_slices_filt3d(:,:,slice_id),0);
        Ct = C;
        Ct2=C;
        
        
        aux_count=1;
        for s = length(C)-lev:length(C)-1
            
            As=0;
            Bs=0;
            Cs=0;
            Ds=0;
            
            size_nt=0;
            
            
            for w = 1:length(C{s})
                
                ave=mean(abs(C{s}{w}(:))/E{s}{w});
                sig=std(abs(C{s}{w}(:))/E{s}{w})^2;
                del=skewness(abs(C{s}{w}(:))/E{s}{w})*(sig^(3/2));
                rho=kurtosis(abs(C{s}{w}(:))/E{s}{w})*(sig^(2));
                
                A_=ave;
                B_=sig+A_^2;
                C_=del+3*B_*A_-2*A_^3;
                D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
                
                size_n=prod(size(C{s}{w}(:)));
                size_nt=size_nt+size_n;
                
                As=A_*size_n+As;
                Bs=B_*size_n+Bs;
                Cs=C_*size_n+Cs;
                Ds=D_*size_n+Ds;
                
                
                
            end
            %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
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
            
            curv3_1(aux_count)=AT;
            curv3_2(aux_count)=var_t;
            curv3_3(aux_count)=skew_t;
            curv3_4(aux_count)=kurt_t;
            curv3_5(aux_count)=rho_t;
            
            aux_count=aux_count+1;
        end
        
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
            D = wrcoef('d',dc_dwt,levels,'db1',lev_3drid);
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
            title(strcat('ridg filt2 for ',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual_rid_3dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                        close(fig);

        end
        
        this=map_3d_slices_filt3d(:,:,slice_id);
        
        anali2a3(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
        anali2a3(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
        anali2a3(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
        anali2a3(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
        
        anali2a3_curv(slice_id,1,:)=curv3_1;
        anali2a3_curv(slice_id,2,:)=curv3_2;
        anali2a3_curv(slice_id,3,:)=curv3_3;
        anali2a3_curv(slice_id,4,:)=curv3_4;
        anali2a3_curv(slice_id,5,:)=curv3_5;
        
        %     if ~isempty(root_data_2d_out)
        %
        %         fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
        %         fwrite(fileID,map_3d_slices_filt3d(:,:,slice_id), 'float32','l');
        %         fclose(fileID);
        %
        %     end
        
        
        if ~isempty(root_data_2d_out)
            
            fileID = fopen(strcat(path_data3_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3a2dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
            fwrite(fileID,map_3d_slices_filt3d(:,:,slice_id), 'float32','l');
            fclose(fileID);
        end
        
        
    end
    
    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3.txt'),anali2a3,'delimiter','\t');
    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3_kurt.txt'),anali2a3_curv,'delimiter','\t');
    
    
    % now we do the depth analysis
    
    if ~ismember(1,sum_depth)
        
        slices_depth=slices/sum_depth;
        %     map_3d_slices_depth=zeros(nb,nb,slices_depth);
        map_3d_slices_filt3d_depth=zeros(nb,nb,slices_depth);
        
        %         map_3d_slices_depth=map_3d_slices_depth+map_3d_slices;
        %         map_3d_slices_filt3d_depth=map_2d_slices_filt2d_depth+map_3d_slices_filt2d;
        %         slice_depth_id=floor(slice_id/sum_depth);
        
        for partition=0:partition3rd_dept-1
                for partition_x=0:partition2d_dept-1
                    for partition_y=0:partition2d_dept-1                                                
                        %             C_aux=map_3d_slices_filt2d(:,:,:);
                        %             C_aux=map_3d_slices((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
                        map_3d_slices_filt3d_depth_=map_3d_slices_filt3d((partition_x)*nb/partition2d_dept+1:(partition_x)*nb/partition2d_dept+nb/partition2d_dept,(partition_y)*nb/partition2d_dept+1:(partition_y)*nb/partition2d_dept+nb/partition2d_dept,(partition)*slices/partition3rd_dept+1:(partition)*slices/partition3rd_dept+slices/partition3rd_dept);
                        
                        % %                         3d ridgelet implementation
%                         
%                         
%                         
%                         
%                         map_3d_slices_depth = repelem(map_3d_slices_depth,2,2,2);                        
                        
%                         [ Radon_Z_ones,interp_Z_info,sample_points_Z_id,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(ones(size(map_3d_slices_depth)));
%                         [ Radon_Z_,interp_Z_info,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(map_3d_slices_depth);
%                         Radon_Z_norm_ = real(Radon_Z_)./real(Radon_Z_ones);
%                         Ridgelet_Z = Ridgelet3d_fromRadon_dev3(Radon_Z_norm_,sample_points_Z_id,lev_3drig);
% %                         dataZ = Radon_Z_norm_(1,:);
%                         Ridgelet_Z(isnan(Ridgelet_Z(:))) = [];
%                         dataZ = Ridgelet_Z;


                        
                        [ Radon_Y_ones,interp_Y_info,sample_points_Y_id,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(ones(size(permute(map_3d_slices_filt3d_depth_,[1 3 2]))));
                        [ Radon_Y_,interp_Y_info,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(permute(map_3d_slices_filt3d_depth_,[1 3 2]));
                        Radon_Y_norm_ = real(Radon_Y_)./real(Radon_Y_ones);
                        Ridgelet_Y = Ridgelet3d_fromRadon_dev3(Radon_Y_norm_,sample_points_Y_id,interp_Y_info,lev_3drid);
%                         dataY = Radon_Y_norm_(1,:);
                        Ridgelet_Y(isnan(Ridgelet_Y(:))) = [];
                        dataY = Ridgelet_Y;
                        
                        [ Radon_X_ones,interp_X_info,sample_points_X_id,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(ones(size(permute(map_3d_slices_filt3d_depth_,[3 2 1]))));
                        [ Radon_X_,interp_X_info,ncx_,ncy_,ncz_] = Radon3d_interp_forwards_dev3(permute(map_3d_slices_filt3d_depth_,[3 2 1]));
                        Radon_X_norm_ = real(Radon_X_)./real(Radon_X_ones);
                        Ridgelet_X = Ridgelet3d_fromRadon_dev3(Radon_X_norm_,sample_points_X_id,interp_X_info,lev_3drid);
%                         dataX = Radon_X_norm_(1,:);
                        Ridgelet_X(isnan(Ridgelet_X(:))) = [];
                        dataX = Ridgelet_X;
                        
%                         data = [dataX,dataY,dataZ];
                        data = [dataX,dataY];
                        
                        
%                         anali2a3_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
%                         anali2a3_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
%                         anali2a3_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
                        anali2a3_depth(1+partition_x+partition_y*partition2d_dept+partition*partition2d_dept*partition2d_dept,:)=[max(data(:)),std(data(:)),max(data(:))/std(data(:)),kurtosis(data(:))];

                        
                        
                        
                    end
                end
        end
        
                dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3_partXY',num2str(partition2d_dept),'_partZ',num2str(partition3rd_dept),'.txt'),anali2a3_depth,'delimiter','\t');


        
%         for slice_depth_id=1:slices_depth
%             
%             %         map_3d_slices_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt2d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
%             map_3d_slices_filt3d_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt3d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
%             
%             
%             if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(3,visual_type))
%                 
%                 %             if ismember(1,visual_in_or_out)
%                 %
%                 %                 fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_depth(:,:,slice_depth_id)); colorbar; axis('image');
%                 %                 xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 %                 ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 %                 set(gca,'FontName','FixedWidth');
%                 %                 set(gca,'FontSize',10);
%                 %                 set(gca,'linewidth',2);
%                 %                 title(strcat('filt2 info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
%                 %
%                 %                 saveas(fig,char(strcat(path_visual_depth_in','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
%                 %
%                 %             end
%                 
%                 if ismember(3,visual_in_or_out)
%                     
%                     fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt3d_depth(:,:,slice_depth_id)); colorbar; axis('image');
%                     xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                     ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                     set(gca,'FontName','FixedWidth');
%                     set(gca,'FontSize',10);
%                     set(gca,'linewidth',2);
%                     title(strcat('filt3d info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
%                     
%                     saveas(fig,char(strcat(path_visual_depth_3dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
%                                 close(fig);
% 
%                 end
%             end
%             
%             C = fdct_wrapping(map_3d_slices_filt3d_depth(:,:,slice_depth_id),0);
%             %     Ct = C;
%             %     Ct2=C;
%             
%             
%             aux_count=1;
%             for s = length(C)-lev:length(C)-1
%                 
%                 As=0;
%                 Bs=0;
%                 Cs=0;
%                 Ds=0;
%                 
%                 size_nt=0;
%                 
%                 
%                 for w = 1:length(C{s})
%                     
%                     ave=mean(abs(C{s}{w}(:))/E{s}{w});
%                     sig=std(abs(C{s}{w}(:))/E{s}{w})^2;
%                     del=skewness(abs(C{s}{w}(:))/E{s}{w})*(sig^(3/2));
%                     rho=kurtosis(abs(C{s}{w}(:))/E{s}{w})*(sig^(2));
%                     
%                     A_=ave;
%                     B_=sig+A_^2;
%                     C_=del+3*B_*A_-2*A_^3;
%                     D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
%                     
%                     size_n=prod(size(C{s}{w}(:)));
%                     size_nt=size_nt+size_n;
%                     
%                     As=A_*size_n+As;
%                     Bs=B_*size_n+Bs;
%                     Cs=C_*size_n+Cs;
%                     Ds=D_*size_n+Ds;
%                     
%                     
%                     
%                 end
%                 %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
%                 AT=As/size_nt;
%                 BT=Bs/size_nt;
%                 CT=Cs/size_nt;
%                 DT=Ds/size_nt;
%                 
%                 sigma_t=BT-AT^2;
%                 delta_t=CT-3*BT*AT+2*AT^3;
%                 rho_t=DT-4*CT*AT+6*BT*AT^2-3*AT^4;
%                 
%                 var_t=sigma_t;
%                 skew_t=delta_t/var_t^(3/2);
%                 kurt_t=rho_t/var_t^2;
%                 
%                 curv2a3_1_depth(aux_count)=AT;
%                 curv2a3_2_depth(aux_count)=var_t;
%                 curv2a3_3_depth(aux_count)=skew_t;
%                 curv2a3_4_depth(aux_count)=kurt_t;
%                 curv2a3_5_depth(aux_count)=rho_t;
%                 
%                 aux_count=aux_count+1;
%             end
%             
%             theta = 0:step_of_degree:180;
%             [R,xp] = radon(map_3d_slices_filt3d_depth(:,:,slice_depth_id),theta);
%             
%             unit=ones(nb);
%             [R_u,xp] = radon(unit,theta);
%             
%             frac_cut=0.5;
%             R_nor=R;
%             R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
%             R_nor(R_u<=nb*frac_cut)=0;
%             
%             %     boudary_removal_factor=2048/nb;
%             
%             n_levels=floor(log2(length(R_nor(:,1))));
%             R_nor_filt=zeros(size(R_nor));
%             for i=1:length(R_nor(1,:))
%                 %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%                 D = wrcoef('d',dc_dwt,levels,'db1',lev_3drid);
%                 %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
%                 %                 D(1:floor((448+200)/boudary_removal_factor))=0;
%                 D(length(xp)-nb*wavel_removal_factor:end)=0;
%                 D(1:nb*wavel_removal_factor)=0;
%                 
%                 R_nor_filt(:,i)=D;
%             end
%             %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
%             %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
%             R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
%             R_nor_filt(1:nb*wavel_removal_factor,:)=[];
%             
%             if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(3,visual_type))&ismember(3,visual_in_or_out)
%                 
%                 %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
%                 fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
%                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
%                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 set(gca,'FontName','FixedWidth');
%                 set(gca,'FontSize',10);
%                 set(gca,'linewidth',2);
%                 title(strcat('ridg filt2d for ',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
%                 
%                 saveas(fig,char(strcat(path_visual_rid_depth_3dfilt','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
%                             close(fig);
% 
%             end
%             
%             this=map_3d_slices_filt3d_depth(:,:,slice_depth_id);
%             
%             anali2a3_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
%             anali2a3_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
%             anali2a3_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
%             anali2a3_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
%             
%             anali2a3_curv_depth(slice_id,1,:)=curv2a3_1_depth;
%             anali2a3_curv_depth(slice_id,2,:)=curv2a3_2_depth;
%             anali2a3_curv_depth(slice_id,3,:)=curv2a3_3_depth;
%             anali2a3_curv_depth(slice_id,4,:)=curv2a3_4_depth;
%             anali2a3_curv_depth(slice_id,5,:)=curv2a3_5_depth;
%         end
        
%         dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3_depth',num2str(sum_depth),'.txt'),anali2a3_depth,'delimiter','\t');
%         dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali2a3_kurt_depth',num2str(sum_depth),'.txt'),anali2a3_curv_depth,'delimiter','\t');
    end
end


%doing the figures for ai analysis

if ismember(4,visual_type)|ismember(5,visual_type)
    for slice_id=1:slices
        
        if ismember(slice_id,snapshot)
            
            if ismember(2,visual_in_or_out)
                
                this=map_3d_slices_filt2d(:,:,slice_id);                
                this=this/(max(this(:)));
                
                thisp=map_3d_slices_filt2d(:,:,mod(slice_id+1-1,slices)+1);
                thisp=thisp/(max(thisp(:)));
                
                thism=map_3d_slices_filt2d(:,:,mod(slice_id-1-1,slices)+1);
                thism=thism/(max(thism(:)));
                
                colorfigure(:,:,1)=thism;
                colorfigure(:,:,2)=this;
                colorfigure(:,:,3)=thisp;
                
                if ismember(4,visual_type)
                
                fig=figure('Visible', 'off');
                set(gcf, 'Position', [0 0 nb-1 nb-1]);
                hold on;
                axes('Position',[0 0 1 1],'Visible','off');
                imagesc([0 nb-1],[0 nb-1],repmat(this,[1 1 3]));
                saveas(fig,char(strcat(path_visual_2dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                close(fig);
                
                end
                
                if ismember(5,visual_type)
                
                fig=figure('Visible', 'off');
                set(gcf, 'Position', [0 0 nb-1 nb-1]);
                hold on;
                axes('Position',[0 0 1 1],'Visible','off');
                imagesc([0 nb-1],[0 nb-1],colorfigure);
                saveas(fig,char(strcat(path_visual_2dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_col_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                close(fig);
                
                end
                
                
            end
            
            if ismember(3,visual_in_or_out)
                
                this=map_3d_slices_filt3d(:,:,slice_id);                
                this=this/(max(this(:)));
                
                thisp=map_3d_slices_filt3d(:,:,mod(slice_id+1-1,slices)+1);
                thisp=thisp/(max(thisp(:)));
                
                thism=map_3d_slices_filt3d(:,:,mod(slice_id-1-1,slices)+1);
                thism=thism/(max(thism(:)));
                
                colorfigure(:,:,1)=thism;
                colorfigure(:,:,2)=this;
                colorfigure(:,:,3)=thisp;
                
                if ismember(4,visual_type)
                
                fig=figure('Visible', 'off');
                set(gcf, 'Position', [0 0 nb-1 nb-1]);
                hold on;
                axes('Position',[0 0 1 1],'Visible','off');
                imagesc([0 nb-1],[0 nb-1],repmat(this,[1 1 3]));
                saveas(fig,char(strcat(path_visual_3dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                close(fig);
                
                end
                
                if ismember(5,visual_type)
                
                fig=figure('Visible', 'off');
                set(gcf, 'Position', [0 0 nb-1 nb-1]);
                hold on;
                axes('Position',[0 0 1 1],'Visible','off');
                imagesc([0 nb-1],[0 nb-1],colorfigure);
                saveas(fig,char(strcat(path_visual_3dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_col_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
                close(fig);
                
                end
                
                
            end
            
        end
    end
end




%doing the figures for ai analysis, depth

if ~ismember(1,sum_depth)
    if ismember(4,visual_type)|ismember(5,visual_type)
        
        slices_depth=slices/sum_depth;
        
        for slice_depth_id=1:slices_depth
            
            if ismember(2,visual_in_or_out)
                
                this=map_3d_slices_filt2d_depth(:,:,slice_depth_id);
                this=this/(max(this(:)));
                
                thisp=map_3d_slices_filt2d_depth(:,:,mod(slice_depth_id+1-1,slices_depth)+1);
                thisp=thisp/(max(thisp(:)));
                
                thism=map_3d_slices_filt2d_depth(:,:,mod(slice_depth_id-1-1,slices_depth)+1);
                thism=thism/(max(thism(:)));
                
                colorfigure(:,:,1)=thism;
                colorfigure(:,:,2)=this;
                colorfigure(:,:,3)=thisp;
                
                if ismember(4,visual_type)
                    
                    fig=figure('Visible', 'off');
                    set(gcf, 'Position', [0 0 nb-1 nb-1]);
                    hold on;
                    axes('Position',[0 0 1 1],'Visible','off');
                    imagesc([0 nb-1],[0 nb-1],repmat(this,[1 1 3]));
                    saveas(fig,char(strcat(path_visual_depth_2dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id),'.png')));
                    close(fig);
                    
                end
                
                if ismember(5,visual_type)
                    
                    fig=figure('Visible', 'off');
                    set(gcf, 'Position', [0 0 nb-1 nb-1]);
                    hold on;
                    axes('Position',[0 0 1 1],'Visible','off');
                    imagesc([0 nb-1],[0 nb-1],colorfigure);
                    saveas(fig,char(strcat(path_visual_depth_2dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_col_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id),'.png')));
                    close(fig);
                    
                end
                
                
            end
            
            if ismember(3,visual_in_or_out)
                
                this=map_3d_slices_filt3d_depth(:,:,slice_depth_id);
                this=this/(max(this(:)));
                
                thisp=map_3d_slices_filt3d_depth(:,:,mod(slice_depth_id+1-1,slices_depth)+1);
                thisp=thisp/(max(thisp(:)));
                
                thism=map_3d_slices_filt3d_depth(:,:,mod(slice_depth_id-1-1,slices_depth)+1);
                thism=thism/(max(thism(:)));
                
                colorfigure(:,:,1)=thism;
                colorfigure(:,:,2)=this;
                colorfigure(:,:,3)=thisp;
                
                if ismember(4,visual_type)
                    
                    fig=figure('Visible', 'off');
                    set(gcf, 'Position', [0 0 nb-1 nb-1]);
                    hold on;
                    axes('Position',[0 0 1 1],'Visible','off');
                    imagesc([0 nb-1],[0 nb-1],repmat(this,[1 1 3]));
                    saveas(fig,char(strcat(path_visual_depth_3dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id),'.png')));
                    close(fig);
                    
                end
                
                if ismember(5,visual_type)
                    
                    fig=figure('Visible', 'off');
                    set(gcf, 'Position', [0 0 nb-1 nb-1]);
                    hold on;
                    axes('Position',[0 0 1 1],'Visible','off');
                    imagesc([0 nb-1],[0 nb-1],colorfigure);
                    saveas(fig,char(strcat(path_visual_depth_3dfilt_fig','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_col_2dproj_3dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id),'.png')));
                    close(fig);
                    
                end
                
                
            end
            
        end
    end
    
end


% 

% path_out=string(strcat(strcat(root_data_2d_anali_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
% string(strcat(root_data_2d_anali_out,spec,aux_path))
% string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
% mkdir(char(strcat(root_data_2d_anali_out,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
% % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
% 
% % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
% dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl','.txt'),anali,'delimiter','\t');
% dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_2a3dcurv_sl','.txt'),anali_2a3d,'delimiter','\t');


cd('../wake_detection/slices_curv_2d/');


end

