function [  ] = curve_2d_data_anali_out_from_2d_slices_out_real( root,root_data_2d_in,root_data_2d_out,root_anali_2d_out,root_visual_2d,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,slices,lev,lev_rid,sigma,step_of_degree,wavel_removal_factor,snapshot,visual_type,visual_in_or_out,sum_depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/';
% % root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_1024_2dcabs_s5lv3_data/';
% root_data_2d_out='';
% root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcabs_s5lv2_anali/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcabs_s5lv2_visual/';
% % spec='4Mpc_2048c_1024p_zi63_nowakem';
% spec='4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m';
% % aux_path='/sample3007/';
% aux_path='/sample3007/half_lin_cutoff_half_tot_pert_nvpw/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% pivot=[0,0,0];
% rot_angle=[1.5708,0,0];
% slices=32;
% lev=3;
% lev_rid=1;
% sigma=5;
% step_of_degree=1;
% wavel_removal_factor=1/2;
% % snapshot=[];
% % snapshot=[13,28]*(128/32);
% snapshot=[9,29]*(slices/32);
% visual_type=[1:2]; %if 1, shows the 2d proj; if 2 shows the ridgelet transformation
% visual_in_or_out=[1,2]; %if 1 do visualization of the input, if 2 of the output
% %  sum_depth=1;
% sum_depth=4;


% filename_read_path=_2dproj_z3_data_sl;
% 
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% % % addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

% record the path to data in

% path_data_in=string(strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));


%making path to data out

if ~isempty(root_data_2d_out)
    
    path_data_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/'));
    % string(strcat(root_data_2d_out,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    mkdir(char(strcat(root_data_2d_out,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')));
    
end


%making path to analisis out

path_anali_out=string(strcat(strcat(root_anali_2d_out,spec,aux_path),'anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/'));
% string(strcat(root_data_2d_anali,spec,aux_path))
% string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
mkdir(char(strcat(root_anali_2d_out,spec,aux_path)),char(strcat('anali/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')));

%making path to analisis depth out

if ~ismember(1,sum_depth)
    
    path_anali_depth_out=string(strcat(strcat(root_anali_2d_out,spec,aux_path),'anali_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/'));
    % string(strcat(root_data_2d_anali,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    mkdir(char(strcat(root_anali_2d_out,spec,aux_path)),char(strcat('anali_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')));
    
end

%making path to visualization out

if ~isempty(snapshot)
    
    
    if ismember(1,visual_in_or_out)
        
        path_visual_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')));
        
        
    end
    
    if ismember(2,visual_in_or_out)
        
        path_visual=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/'));
        % string(strcat(root_data_2d_anali,spec,aux_path))
        % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
        mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')));
        
        
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
            
            path_visual_depth_in=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2d/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'_in/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2d/dm/')));
            
            
        end
        
        if ismember(2,visual_in_or_out)
            
            path_visual_depth=string(strcat(strcat(root_visual_2d,spec,aux_path),'visual_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2d/dm/'));
            % string(strcat(root_data_2d_anali,spec,aux_path))
            % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
            mkdir(char(strcat(root_visual_2d,spec,aux_path)),char(strcat('visual_depth_',num2str(sum_depth),'/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2d/dm/')));
            
            
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

filename_read=strcat('_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_sl');

% psi=rot_angle(3);
% theta=rot_angle(2);
% phi=rot_angle(1);

nb=np*resol_factor/lenght_factor;

%
% F=zeros(nc);
% C_zero = fdct_wrapping(F,0);
F = ones(nb);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%X = F * sqrt(prod(size(F)));
%C = fdct_wrapping(X,0);
C = fdct_wrapping(X,0);
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

slices_depth=slices/sum_depth;
slice_depth_count=1;
map_2d_slices_depth=zeros(nb,nb);
map_2d_slices_filt2d_depth=zeros(nb,nb);

for slice_id=1:slices
    
    
    
    filename_read_path=strcat(strcat(root_data_2d_in,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/',filename_read,num2str(slice_id),'.bin')
    fid = fopen(filename_read_path);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    %             map = imresize(map,new_nc/nc,'triangle');
    
    
    map(map<=1)=1;%to remove problem with holes
    
    map_3d_slices=log(map);
    
    
    
    
    
    
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
    
    C = fdct_wrapping(map_3d_slices,0);
    Ct = C;
    for s = 1:length(C)
        for w = 1:length(C{s})
            Ct{s}{w} = C_zero{s}{w};
        end
    end
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        
        %                 thresh=0;
        thresh = sigma + sigma*(s == length(C));
        for w = 1:length(C{s})
            %                     Ct{s}{w} = C{s}{w};
            %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
            %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
%             Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
            %                 Ct{s}{w} = C{s}{w};
            
            Ct{s}{w} = C{s}{w};
            
            
            %              curv(slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
            
            
            
%             ave=mean(abs(Ct{s}{w}(:)));
%             sig=std(abs(Ct{s}{w}(:)))^2;
%             del=skewness(abs(Ct{s}{w}(:)))*(sig^(3/2));
%             rho=kurtosis(abs(Ct{s}{w}(:)))*(sig^(2));

            ave=mean(abs(Ct{s}{w}(:))/E{s}{w});
            sig=std(abs(Ct{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(Ct{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(Ct{s}{w}(:))/E{s}{w})*(sig^(2));

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
    
    
    
    
    
    map_3d_slices_filt2d = real(ifdct_wrapping(Ct,0));
    
    
    %     strcat(root_data_2d_in,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/',filename_read,num2str(slice_id),'.bin')
    
    
    if ismember(slice_id,snapshot)&(ismember(1,visual_type))
        
        if ismember(1,visual_in_or_out)
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices); colorbar; axis('image');
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
            
            fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_3d_slices_filt2d); colorbar; axis('image');
            xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
            set(gca,'FontName','FixedWidth');
            set(gca,'FontSize',10);
            set(gca,'linewidth',2);
            title(strcat('filt2 info:',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
            
            saveas(fig,char(strcat(path_visual','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
            close(fig);
        end
    end
    
%     anali(slice_id,1,:)=[max(map_3d_slices_filt2d(:)),std(map_3d_slices_filt2d(:)),max(map_3d_slices_filt2d(:))/std(map_3d_slices_filt2d(:)),kurtosis(kurtosis(map_3d_slices_filt2d)),kurtosis(map_3d_slices_filt2d(:))];
    
    % theta = 0:180/nc:180;
    
    
   
    
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filt2d,theta);
    
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
        D = wrcoef('d',dc_dwt,levels,'db1',lev_rid);
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
    
    
    if ismember(slice_id,snapshot)&(ismember(1,visual_type))&ismember(2,visual_in_or_out)
        
        %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
        fig=figure('Visible', 'off'); imagesc(0:step_of_degree:180,(size_box/nb)*[1:nb],R_nor_filt);colorbar;
        xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
        ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
        set(gca,'FontName','FixedWidth');
        set(gca,'FontSize',10);
        set(gca,'linewidth',2);
        title(strcat('ridg filt2 for ',' ',aux_path(2:11),';slice =',num2str(slice_id),';G\mu = ',num2str(Gmu)));
        
        saveas(fig,char(strcat(path_visual_rid','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_id),'.png')));
        close(fig);
    end
    
    anali(slice_id,1,:)=[max(map_3d_slices_filt2d(:)),std(map_3d_slices_filt2d(:)),max(map_3d_slices_filt2d(:))/std(map_3d_slices_filt2d(:)),kurtosis(kurtosis(map_3d_slices_filt2d)),kurtosis(map_3d_slices_filt2d(:))];
    anali(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
    anali_curv(slice_id,1,:)=curv_1;
    anali_curv(slice_id,2,:)=curv_2;
    anali_curv(slice_id,3,:)=curv_3;
    anali_curv(slice_id,4,:)=curv_4;
    anali_curv(slice_id,5,:)=curv_5;
    
    %     path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    % string(strcat(root_data_2d_anali,spec,aux_path))
    % string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
    
    % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
    
    if ~isempty(root_data_2d_out)
        
        fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
        fwrite(fileID,map_3d_slices_filt2d, 'float32','l');
        fclose(fileID);
        
    end
    
    
    %now the slice stuff
    
    if ~ismember(1,sum_depth)
        
        map_2d_slices_depth=map_2d_slices_depth+map_3d_slices;
        map_2d_slices_filt2d_depth=map_2d_slices_filt2d_depth+map_3d_slices_filt2d;
        slice_depth_id=floor(slice_id/sum_depth);
        
        if slice_depth_count==sum_depth
            
            if ismember(slice_depth_id,floor(snapshot/sum_depth))&(ismember(1,visual_type))
                
                if ismember(1,visual_in_or_out)
                    
                    fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_2d_slices_depth(:,:)); colorbar; axis('image');
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
                    
                    fig=figure('Visible', 'off'); imagesc((size_box/nb)*[1:nb],(size_box/nb)*[1:nb],map_2d_slices_filt2d_depth(:,:)); colorbar; axis('image');
                    xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                    ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                    set(gca,'FontName','FixedWidth');
                    set(gca,'FontSize',10);
                    set(gca,'linewidth',2);
                    title(strcat('filt2 info:',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                    
                    saveas(fig,char(strcat(path_visual_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                    close(fig);
                end
                
            end
            
            C = fdct_wrapping(map_2d_slices_filt2d_depth,0);
    Ct = C;
    for s = 1:length(C)
        for w = 1:length(C{s})
            Ct{s}{w} = C_zero{s}{w};
        end
    end
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        
        %                 thresh=0;
        for w = 1:length(C{s})
            %                     Ct{s}{w} = C{s}{w};
            %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
            %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
            Ct{s}{w} = C{s}{w};
            %                 Ct{s}{w} = C{s}{w};
            %              curv(slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
            
%             ave=mean(abs(Ct{s}{w}(:)));
%             sig=std(abs(Ct{s}{w}(:)))^2;
%             del=skewness(abs(Ct{s}{w}(:)))*(sig^(3/2));
%             rho=kurtosis(abs(Ct{s}{w}(:)))*(sig^(2));

            ave=mean(abs(Ct{s}{w}(:))/E{s}{w});
            sig=std(abs(Ct{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(Ct{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(Ct{s}{w}(:))/E{s}{w})*(sig^(2));

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
        
        curv_1_depth(aux_count)=AT;
        curv_2_depth(aux_count)=var_t;
        curv_3_depth(aux_count)=skew_t;
        curv_4_depth(aux_count)=kurt_t;
        curv_5_depth(aux_count)=rho_t;
        
        aux_count=aux_count+1;
    end
            
            
            this=map_2d_slices_filt2d_depth(:,:);
            
            
            % theta = 0:180/nc:180;
            
            
            
            theta = 0:step_of_degree:180;
            [R,xp] = radon(map_2d_slices_filt2d_depth(:,:),theta);
            
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
                D = wrcoef('d',dc_dwt,levels,'db1',lev_rid);
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
                title(strcat('ridg filt2 for ',' ',aux_path(2:11),';sliceSum',num2str(sum_depth),' =',num2str(slice_depth_id*sum_depth),';G\mu = ',num2str(Gmu)));
                
                saveas(fig,char(strcat(path_visual_rid_depth','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv&rid_z',num2str(z_glob),'_visual_sl',num2str(slice_depth_id*sum_depth),'.png')));
                close(fig);
            end
            
            anali_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
            anali_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
            anali_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
            
            anali_curv_depth(slice_depth_id,1,:)=curv_1_depth;
            anali_curv_depth(slice_depth_id,2,:)=curv_2_depth;
            anali_curv_depth(slice_depth_id,3,:)=curv_3_depth;
            anali_curv_depth(slice_depth_id,4,:)=curv_4_depth;
            anali_curv_depth(slice_depth_id,5,:)=curv_5_depth;
            
        end
        
        
        
        slice_depth_count=slice_depth_count+1;
        
        if slice_depth_count==sum_depth+1
            
            slice_depth_count=1;
            map_2d_slices_depth=zeros(nb,nb);
            map_2d_slices_filt2d_depth=zeros(nb,nb);
            
        end
        
        
        
        
        
        
    end
    
end

if ~ismember(1,sum_depth)
    dlmwrite(strcat(path_anali_depth_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_depth.txt'),anali_depth,'delimiter','\t');
    dlmwrite(strcat(path_anali_depth_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv_depth.txt'),anali_curv_depth,'delimiter','\t');
end

dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv.txt'),anali_curv,'delimiter','\t');


cd('../wake_detection/curvelet/curvelab/2d/');

end

