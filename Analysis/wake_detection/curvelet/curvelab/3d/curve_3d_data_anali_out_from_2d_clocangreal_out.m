function [  ] = curve_3d_data_anali_out_from_2d_clocangreal_out( root,root_data_2d_in,root_data_2d_out,root_anali_2d_out,root_visual_2d,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,slices,lev_3d,lev_3drig,sigma,step_of_degree,wavel_removal_factor,snapshot,visual_type,visual_in_or_out,partition2d,sum_depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% 
root='/home/asus/Dropbox/extras/storage/graham/ht/';
root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclar_l2lr1na256_data/';
% % root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps256_1024_3dcurv_s5lv3_data/';
root_data_2d_out='';
root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps256_1024_2dclar-l2lr1na256_to_3dcurv-l1lr1_anali/';
root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps256_1024_2dclar-l2lr1na256_to_3dcurv-l1lr1_visual/';
% spec='4Mpc_2048c_1024p_zi63_nowakem';
spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% aux_path='/sample3003/';
aux_path='/sample3007/half_lin_cutoff_half_tot_pert_nvpw/';
aux_path_out='';
filename='3.000xv0.dat';
lenght_factor=1;
resol_factor=1;
pivot=[0,0,0];
rot_angle=[1.5708,0,0];
slices=32;
% lev_2d=2;
% lev_2drig=1;
lev_3d=1;
lev_3drig=1;
sigma=5;
step_of_degree=1*(180/256);
wavel_removal_factor=1/2;
snapshot=[];
% snapshot=[13,28]*(256/32);
% snapshot=[9,29]*(128/32);
visual_type=[1:2]; %if 1, shows the 2d proj; if 2 shows the ridgelet transformation
% visual_in_or_out=[1,2]; %if 1 do visualization of the input, if 2 of the output
 visual_in_or_out=[2];
 partition2d=4;
%  sum_depth=1;
 sum_depth=32;

% filename_read_path=_2dproj_z3_data_sl;

addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

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


F=zeros(nc);
C_zero = fdct_wrapping(F,0);


F2 = ones(new_nc,new_nc,slices);
% F2 = ones(new_nc,new_nc,slices);
X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
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
F2=zeros(new_nc,new_nc,slices/2);

C_zero2 = fdct3d_forward(F2);


map_3d_slices=zeros(nb,nb,slices);
map_3d_slices_filt3d=zeros(nc_anal3d,nc_anal3d,slices);

for slice_id=1:slices
    
    
    
    filename_read_path=strcat(strcat(root_data_2d_in,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/',filename_read,num2str(slice_id),'.bin')
    fid = fopen(filename_read_path);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    %             map = imresize(map,new_nc/nc,'triangle');
    
    
%     map(map<=1)=1;%to remove problem with holes
    
    
    map_3d_slices(:,:,slice_id)=map;
    
    
end



for partition_x=0:partition2d-1
    for partition_y=0:partition2d-1
        
        C_aux=map_3d_slices((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,1:slices);
                       
        
         C = fdct3d_forward(C_aux);
            
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero2{s}{w};
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
            
            for w = 1:length(C{s})
                Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E2{s}{w})
                %                                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                %                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                %                   Ct{s}{w} = C{s}{w};
                %                   curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                
                
                ave=mean(abs(Ct{s}{w}(:))/E2{s}{w});
                    sig=std(abs(Ct{s}{w}(:))/E2{s}{w})^2;
                    del=skewness(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(3/2));
                    rho=kurtosis(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(2));
                    
                    A_=ave;
                    B_=sig+A_^2;
                    C_=del+3*B_*A_-2*A_^3;
                    D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
                    
                    size_n=prod(size(Ct{s}{w}(:)));
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
        map_3d_slices_filt3d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,1:slices) = real(fdct3d_inverse(Ct));
        
    end
end



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

dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali.txt'),anali,'delimiter','\t');




% now we do the depth analysis

if ~ismember(1,sum_depth)
    
    slices_depth=slices/sum_depth;
    map_3d_slices_depth=zeros(nb,nb,slices_depth);
    map_3d_slices_filt3d_depth=zeros(nb,nb,slices_depth);
    
    for slice_depth_id=1:slices_depth
        
        map_3d_slices_depth(:,:,slice_depth_id)=sum(map_3d_slices(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        map_3d_slices_filt3d_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt3d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        
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
        end
%         
%         this=map_3d_slices_filt3d_depth(:,:,slice_depth_id);
%         
%         anali_depth(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
%         
        % theta = 0:180/nc:180;
        
        
        
        theta = 0:step_of_degree:180;
        [R,xp] = radon(map_3d_slices_filt3d_depth(:,:,slice_depth_id),theta);
        
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
            D = wrcoef('d',dc_dwt,levels,'db1',lev_2d);
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
        
        
        if ismember(slice_id,floor(snapshot/sum_depth))&(ismember(2,visual_type))&ismember(2,visual_in_or_out)
            
            %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
            fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
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
        
        anali_depth(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];      
        anali_depth(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
        anali_depth(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
        anali_depth(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
        
        
        
        %     path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
        % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
        
        % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
        
        
        
    end

    dlmwrite(strcat(path_anali_depth_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali.txt'),anali_depth,'delimiter','\t');
    
end



cd('../wake_detection/curvelet/curvelab/3d/');

end

