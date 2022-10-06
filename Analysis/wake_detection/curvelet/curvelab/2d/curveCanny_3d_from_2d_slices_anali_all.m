function [  ] = curveCanny_3d_from_2d_slices_anali_all(  )
% 
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali/';
% root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali_all/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_visual_all/';
% 
%clearvars;

lenght_factor=1;
resol_factor=0.5;
pivot=[0,0,0];
rot_angle=[1.5708,0,0];
filename='3.000xv0.dat'; GMU_pre=4;
slices=32;
nc=512;
sum_depth=4;
L_box=4;
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_anali/');
root_anali_2d_in=strcat('/home/disraelcunha/Documents/graham/cubep3m/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_anali/');
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_anali_testRid/');
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_anali_testRidAug/');
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar_d2-l1lr1_anali/');
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparcurv-l1lr1_anali/');
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
root='/home/disraelcunha/Documents/graham/cubep3m/ht/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps16_512_2dclar-l2lr1na128_to_3dparcurv-l1lr1_anali/';
% root_anali_2d_in=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclar-l2lr1na256_to_3dparcurv-l1lr2_anali/');
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcr0_l2lr1ap256_anali/';
%root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali_all/';
%root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_visual_all/';






% lev=2;
% sigma=5;
% step_of_degree=1;
% wavel_removal_factor=1/2;
% % snapshot=[];
% % snapshot=[13,28]*(128/32);
% snapshot=[9,29]*(128/32);
% visual_type=[1:2]; %if 1, shows the 2d proj; if 2 shows the ridgelet transformation
%  visual_in_or_out=[1,2]; %if 1 do visualization of the input, if 2 of the output

nowake_samples_ids=[11:60];
wake_samples_ids=[11:60];

% nowake_samples_ids=[1:10];
% wake_samples_ids=[1:10];


% wake_samples_ids=[1:10];


addpath('../../../../processing');
% addpath('../../../../preprocessing');



path_specs_in=strcat(root_anali_2d_in);

specs_nowake=dir(strcat(root_anali_2d_in,'/',num2str(L_box),'*nowake*'));
specs_nowake={specs_nowake.name};
specs_wake=dir(strcat(root_anali_2d_in,'/',num2str(L_box),'*wakeGmu',num2str(GMU_pre),'*'));
specs_wake={specs_wake.name};

specs_list=dir(strcat(path_specs_in,'/',num2str(L_box),'*'));
specs_list={specs_list.name};
specs_list=sort_nat(specs_list);
specs_list=specs_list(1:end);

display(specs_list);

spec_nowake=specs_nowake{1};
path_samples_in=strcat(root_anali_2d_in,spec_nowake);
sample_list_all=dir(strcat(path_samples_in,'/sample*'));
sample_list=sample_list_all(nowake_samples_ids);
sample_list=strcat('/',{sample_list.name},'/');
sample_list_nowake=sort_nat(sample_list)

spec_wake=specs_wake{1};
path_samples_in=strcat(root_anali_2d_in,spec_wake);
sample_list_all=dir(strcat(path_samples_in,'/sample*'));
sample_list=sample_list_all(wake_samples_ids);
% sample_list_short=strcat('/',{sample_list.name},'/');
% sample_list_short=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v0p55');
sample_list_short=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
% sample_list_short=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_wrong/');
% sample_list_short=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert/');
% sample_list_short=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v1/');


% sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw/');
% sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v0p55/');
sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
% sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert/');
% sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_wrong/');
% sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw_v1/');

sample_list_wake=sort_nat(sample_list)
sample_list_wake_short=sort_nat(sample_list_short);

%
% fig_pk=figure;
% set(gcf, 'Position', [0 0 300 800]);
% ax_pk=axes(fig_pk);

%  fig=figure('Visible', 'off')
fig1=figure;
fig2=figure;
fig3=figure;
fig4=figure;
% fig5=figure;

fig1_Canny=figure;
fig2_Canny=figure;
fig3_Canny=figure;
fig4_Canny=figure;

fig1_curv=figure;
fig2_curv=figure;
fig3_curv=figure;
fig4_curv=figure;
fig5_curv=figure;


ax1=axes(fig1);
ax2=axes(fig2);
ax3=axes(fig3);
ax4=axes(fig4);
% ax5=axes(fig5);

ax1_Canny=axes(fig1_Canny);
ax2_Canny=axes(fig2_Canny);
ax3_Canny=axes(fig3_Canny);
ax4_Canny=axes(fig4_Canny);
% ax5=axes(fig5);




ax1_curv=axes(fig1_curv);
ax2_curv=axes(fig2_curv);
ax3_curv=axes(fig3_curv);
ax4_curv=axes(fig4_curv);
ax5_curv=axes(fig5_curv);

if ~ismember(1,sum_depth)
    
    fig1_depth=figure;
    fig2_depth=figure;
    fig3_depth=figure;
    fig4_depth=figure;
    % fig5=figure;
    
    fig1_Canny_depth=figure;
    fig2_Canny_depth=figure;
    fig3_Canny_depth=figure;
    fig4_Canny_depth=figure;
    %
    %     fig1_curv_depth=figure;
    %     fig2_curv_depth=figure;
    %     fig3_curv_depth=figure;
    %     fig4_curv_depth=figure;
    %     fig5_curv_depth=figure;
    
    
    ax1_depth=axes(fig1_depth);
    ax2_depth=axes(fig2_depth);
    ax3_depth=axes(fig3_depth);
    ax4_depth=axes(fig4_depth);
    % ax5=axes(fig5);
    
    ax1_Canny_depth=axes(fig1_Canny_depth);
    ax2_Canny_depth=axes(fig2_Canny_depth);
    ax3_Canny_depth=axes(fig3_Canny_depth);
    ax4_Canny_depth=axes(fig4_Canny_depth);
    
    %
    %     ax1_curv_depth=axes(fig1_curv_depth);
    %     ax2_curv_depth=axes(fig2_curv_depth);
    %     ax3_curv_depth=axes(fig3_curv_depth);
    %     ax4_curv_depth=axes(fig4_curv_depth);
    %     ax5_curv_depth=axes(fig5_curv_depth);
    
end

cd('../../../../preprocessing')
% [xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,sample_list_wake{1} );
display(strcat(root,spec_wake,sample_list_wake_short{1}))
[~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec_wake,sample_list_wake_short{1} );

z_string=char(filename);
z_string=z_string(1:end-7);
z=str2num(z_string);
z_glob=z;



for w_nw=1:2
    % for w_nw=2
    
    if w_nw==1
        % [xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,sample_list_wake{1} );
        [~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec_nowake,sample_list_nowake{1} );
        spec=specs_nowake{1};
        sample_list=sample_list_nowake;
        coul='b';
    else
        display(strcat(root,spec_wake,sample_list_wake{1} ))
        % [xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,sample_list_wake{1} );
        [~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec_wake,sample_list_wake{1} );
        spec=specs_wake{1};
        sample_list=sample_list_wake;
        coul='r';
    end
    
    % if w_nw==1
    %              signal_nw=[];
    %     else
    %             signal_w=[];
    % end
    
    for sample = 1:length(sample_list)
        
        path_in=strcat(strcat(root_anali_2d_in,spec,char(sample_list(sample))),'anali/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')
        
        filename=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali.txt')
        filename_Canny=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_analiCanny.txt')
        filename_curv=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali3_curv.txt')
        
        %         filename=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali.txt')
        
        info = dlmread(filename);
        info_Canny = dlmread(filename_Canny);
        info_curv = dlmread(filename_curv);
        %         lv_sz=prod(size(info_curv))/(slices*5);
        lv_sz=prod(size(info_curv))/(5);
        
        anali(w_nw,sample,:,:,:)=reshape(info,slices,4,5);
        anali_Canny(w_nw,sample,:,:,:)=reshape(info_Canny,slices,4,5);
        anali_curv(w_nw,sample,:,:)=reshape(info_curv,5,lv_sz);
        
        a(:)=max(anali(w_nw,sample,:,1,:),[],3);
        b(:)=max(anali(w_nw,sample,:,2,:),[],3);
        c(:)=max(anali(w_nw,sample,:,3,:),[],3);
        d(:)=max(anali(w_nw,sample,:,4,:),[],3);
        
        a_Canny(:)=max(anali_Canny(w_nw,sample,:,1,:),[],3);
        b_Canny(:)=max(anali_Canny(w_nw,sample,:,2,:),[],3);
        c_Canny(:)=max(anali_Canny(w_nw,sample,:,3,:),[],3);
        d_Canny(:)=max(anali_Canny(w_nw,sample,:,4,:),[],3);
        
        a_curv(:)=(anali_curv(w_nw,sample,1,:));
        b_curv(:)=(anali_curv(w_nw,sample,2,:));
        c_curv(:)=(anali_curv(w_nw,sample,3,:));
        d_curv(:)=(anali_curv(w_nw,sample,4,:));
        e_curv(:)=(anali_curv(w_nw,sample,5,:));
        
        plot1{sample}=   plot(ax1,a,coul);
        plot2{sample}=   plot(ax2,b,coul);
        plot3{sample}=   plot(ax3,c,coul);
        plot4{sample}=   plot(ax4,d,coul);
        
        plot1_Canny{sample}=   plot(ax1_Canny,a_Canny,coul);
        plot2_Canny{sample}=   plot(ax2_Canny,b_Canny,coul);
        plot3_Canny{sample}=   plot(ax3_Canny,c_Canny,coul);
        plot4_Canny{sample}=   plot(ax4_Canny,d_Canny,coul);
        
        
        if (lv_sz==1)
            
            plot1_curv{sample}=   scatter(ax1_curv,1,a_curv,coul);
            plot2_curv{sample}=   scatter(ax2_curv,1,b_curv,coul);
            plot3_curv{sample}=   scatter(ax3_curv,1,c_curv,coul);
            plot4_curv{sample}=   scatter(ax4_curv,1,d_curv,coul);
            plot5_curv{sample}=   scatter(ax5_curv,1,e_curv,coul);
            
        else
            
            plot1_curv{sample}=   plot(ax1_curv,a_curv,coul);
            plot2_curv{sample}=   plot(ax2_curv,b_curv,coul);
            plot3_curv{sample}=   plot(ax3_curv,c_curv,coul);
            plot4_curv{sample}=   plot(ax4_curv,d_curv,coul);
            plot5_curv{sample}=   plot(ax5_curv,e_curv,coul);
            
        end
            
        
%         plot1_curv{sample}=   plot(ax1_curv,a_curv,coul);
%         plot2_curv{sample}=   plot(ax2_curv,b_curv,coul);
%         plot3_curv{sample}=   plot(ax3_curv,c_curv,coul);
%         plot4_curv{sample}=   plot(ax4_curv,d_curv,coul);
%         plot5_curv{sample}=   plot(ax5_curv,e_curv,coul);
        
        clearvars a b c d a_Canny b_Canny c_Canny d_Canny a_curv b_curv c_curv d_curv e_curv
        
        hold(ax1,'on');
        hold(ax2,'on');
        hold(ax3,'on');
        hold(ax4,'on');
        
        hold(ax1_Canny,'on');
        hold(ax2_Canny,'on');
        hold(ax3_Canny,'on');
        hold(ax4_Canny,'on');
        
        hold(ax1_curv,'on');
        hold(ax2_curv,'on');
        hold(ax3_curv,'on');
        hold(ax4_curv,'on');
        hold(ax5_curv,'on');
        
        if ~ismember(1,sum_depth)
            
            path_in_depth=strcat(strcat(root_anali_2d_in,spec,char(sample_list(sample))),'anali_depth_',num2str(sum_depth),'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','3d/dm/')
            
            filename_depth=strcat(path_in_depth,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali3_depth',num2str(sum_depth),'.txt')
            filename_Canny_depth=strcat(path_in_depth,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurvCanny_z',num2str(z_glob),'_anali_depth',num2str(sum_depth),'.txt')
            
            %             filename_curv_depth=strcat(path_in_depth,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv_depth.txt')
            
            info_depth = dlmread(filename_depth);
            info_Canny_depth = dlmread(filename_Canny_depth);
            
            %             info_curv_depth = dlmread(filename_curv_depth);
            %             lv_sz_depth=prod(size(info_curv_depth))*sum_depth/(slices*5);
            %

            
% old setup, from z=0 result            
            anali_depth(w_nw,sample,:,:,:)=reshape(info_depth,slices/sum_depth,4,5);
            anali_Canny_depth(w_nw,sample,:,:,:)=reshape(info_Canny_depth,slices/sum_depth,4,5);

% %ne setup for the ridgelet 3d
%             anali_depth(w_nw,sample,:,:,:)=reshape(info_depth,slices/sum_depth,4);
%             anali_Canny_depth(w_nw,sample,:,:,:)=reshape(info_Canny_depth,slices/sum_depth,4);
            
            %             anali_curv_depth(w_nw,sample,:,:,:)=reshape(info_curv_depth,slices/sum_depth,5,lv_sz_depth);
            
            a_depth(:)=max(anali_depth(w_nw,sample,:,1,:),[],3);
            b_depth(:)=max(anali_depth(w_nw,sample,:,2,:),[],3);
            c_depth(:)=max(anali_depth(w_nw,sample,:,3,:),[],3);
            d_depth(:)=max(anali_depth(w_nw,sample,:,4,:),[],3);
            
            a_Canny_depth(:)=max(anali_Canny_depth(w_nw,sample,:,1,:),[],3);
            b_Canny_depth(:)=max(anali_Canny_depth(w_nw,sample,:,2,:),[],3);
            c_Canny_depth(:)=max(anali_Canny_depth(w_nw,sample,:,3,:),[],3);
            d_Canny_depth(:)=max(anali_Canny_depth(w_nw,sample,:,4,:),[],3);
            
            
            
            %
            %             a_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,1,:),[],3);
            %             b_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,2,:),[],3);
            %             c_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,3,:),[],3);
            %             d_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,4,:),[],3);
            %             e_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,5,:),[],3);
            %
            plot1_depth{sample}=   plot(ax1_depth,a_depth,coul);
            plot2_depth{sample}=   plot(ax2_depth,b_depth,coul);
            plot3_depth{sample}=   plot(ax3_depth,c_depth,coul);
            plot4_depth{sample}=   plot(ax4_depth,d_depth,coul);
            
            plot1_Canny_depth{sample}=   plot(ax1_Canny_depth,a_Canny_depth,coul);
            plot2_Canny_depth{sample}=   plot(ax2_Canny_depth,b_Canny_depth,coul);
            plot3_Canny_depth{sample}=   plot(ax3_Canny_depth,c_Canny_depth,coul);
            plot4_Canny_depth{sample}=   plot(ax4_Canny_depth,d_Canny_depth,coul);
            
            %             plot1_curv_depth{sample}=   plot(ax1_curv_depth,a_curv_depth,coul);
            %             plot2_curv_depth{sample}=   plot(ax2_curv_depth,b_curv_depth,coul);
            %             plot3_curv_depth{sample}=   plot(ax3_curv_depth,c_curv_depth,coul);
            %             plot4_curv_depth{sample}=   plot(ax4_curv_depth,d_curv_depth,coul);
            %             plot5_curv_depth{sample}=   plot(ax5_curv_depth,e_curv_depth,coul);
            %
            clearvars a_depth b_depth c_depth d_depth a_Canny_depth b_Canny_depth c_Canny_depth d_Canny_depth a_curv_depth b_curv_depth c_curv_depth d_curv_depth e_curv_depth
            
            hold(ax1_depth,'on');
            hold(ax2_depth,'on');
            hold(ax3_depth,'on');
            hold(ax4_depth,'on');
            
            hold(ax1_Canny_depth,'on');
            hold(ax2_Canny_depth,'on');
            hold(ax3_Canny_depth,'on');
            hold(ax4_Canny_depth,'on');
            
            %             hold(ax1_curv_depth,'on');
            %             hold(ax2_curv_depth,'on');
            %             hold(ax3_curv_depth,'on');
            %             hold(ax4_curv_depth,'on');
            %             hold(ax5_curv_depth,'on');
            
            
        end
        
    end
    
end

set(ax1, 'YScale', 'log');
title(ax1,'radon of the original map ');
set(ax2, 'YScale', 'log');
title(ax2,'radon of the 2dcurv-filtered map');
set(ax3, 'YScale', 'log');
title(ax3,'1dwavel over radon of the 2dcurv-filtered map');
set(ax4, 'YScale', 'log');
title(ax4,'ridgelet normalized of the 2dcurv-filtered map');


set(ax1_Canny, 'YScale', 'log');
title(ax1_Canny,'radon Canny of the original map ');
set(ax2_Canny, 'YScale', 'log');
title(ax2_Canny,'radon Canny of the 2dcurv-filtered map');
set(ax3_Canny, 'YScale', 'log');
title(ax3_Canny,'1dwavel Canny over radon of the 2dcurv-filtered map');
set(ax4_Canny, 'YScale', 'log');
title(ax4_Canny,'ridgelet Canny normalized of the 2dcurv-filtered map');



set(ax1_curv, 'YScale', 'log');
title(ax1_curv,'average normalised curvelet abs coef fast');
set(ax2_curv, 'YScale', 'log');
title(ax2_curv,'std normalised curvelet abs coef fast');
set(ax3_curv, 'YScale', 'log');
title(ax3_curv,'skewness normalised curvelet abs coef fast');
set(ax4_curv, 'YScale', 'log');
title(ax4_curv,'kurtosis normalised curvelet abs coef fast');
set(ax5_curv, 'YScale', 'log');
title(ax5_curv,'4th moment normalised curvelet abs coef fast');


if ~ismember(1,sum_depth)
    
    set(ax1_depth, 'YScale', 'log');
    title(ax1_depth,strcat('radon of the original map, slice= ',num2str(sum_depth)));
    set(ax2_depth, 'YScale', 'log');
    title(ax2_depth,strcat('radon of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    set(ax3_depth, 'YScale', 'log');
    title(ax3_depth,strcat('1dwavel over radon of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    set(ax4_depth, 'YScale', 'log');
    title(ax4_depth,strcat('ridgelet normalized of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    
    
    set(ax1_Canny_depth, 'YScale', 'log');
    title(ax1_Canny_depth,strcat('radon Canny of the original map, slice= ',num2str(sum_depth)));
    set(ax2_Canny_depth, 'YScale', 'log');
    title(ax2_Canny_depth,strcat('radon Canny of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    set(ax3_Canny_depth, 'YScale', 'log');
    title(ax3_Canny_depth,strcat('1dwavel Canny over radon of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    set(ax4_Canny_depth, 'YScale', 'log');
    title(ax4_Canny_depth,strcat('ridgelet Canny normalized of the 2dcurv-filtered map, slice= ',num2str(sum_depth)));
    
    %
    %     set(ax1_curv_depth, 'YScale', 'log');
    %     title(ax1_curv_depth,strcat('average normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    %     set(ax2_curv_depth, 'YScale', 'log');
    %     title(ax2_curv_depth,strcat('std normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    %     set(ax3_curv_depth, 'YScale', 'log');
    %     title(ax3_curv_depth,strcat('skewness normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    %     set(ax4_curv_depth, 'YScale', 'log');
    %     title(ax4_curv_depth,strcat('kurtosis normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    %     set(ax5_curv_depth, 'YScale', 'log');
    %     title(ax5_curv_depth,strcat('4th moment normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    
end
cd('../wake_detection/curvelet/curvelab/2d/')

% 
% 
% 
% nowake=reshape(permute(anali(1,1:length(sample_list_nowake),:,4,3),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali(2,1:length(sample_list_wake),:,4,3),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_wake),:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[slices,length(sample_list_nowake)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
% 
% % 
% nowake=reshape(permute(anali_Canny(1,1:length(sample_list_nowake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali_Canny(2,1:length(sample_list_wake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_wake),:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[slices,length(sample_list_nowake)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
% 
% 
% nowake=reshape(permute(anali_curv(1,1:length(sample_list_nowake),:,1),[1,3,2,4]),[1,numel(anali_curv(1,1:length(sample_list_nowake),:,1))])
% wake=reshape(permute(anali_curv(2,1:length(sample_list_wake),:,1),[1,3,2,4]),[1,numel(anali_curv(1,1:length(sample_list_wake),:,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[5,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[5,length(sample_list_nowake)])'
% % max_wake_slices_=sort(wake_slices')
% % max_nowake_slices_=sort(nowake_slices')
% % max_wake_slices=max_wake_slices_(end,:)
% % max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(wake_slices)
% mean_nowake=mean(nowake_slices)
% std_wake=std(wake_slices,1)
% std_nowake=std(nowake_slices,1)
% stn_nowake=(nowake_slices-mean_nowake)./std_nowake
% stn_wake=(wake_slices-mean_nowake)./std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)./(std_wake+std_nowake)

% old

nowake=reshape(permute(anali_depth(1,1:length(sample_list_nowake),:,4,3),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_nowake),:,2,1))])
wake=reshape(permute(anali_depth(2,1:length(sample_list_wake),:,4,3),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_wake),:,2,1))])
mean_wake=mean(wake)
mean_nowake=mean(nowake)
std_nowake=std(nowake,1)
stn_nowake=(nowake-mean_nowake)/std_nowake
stn_wake=(wake-mean_nowake)/std_nowake
mean_stn=mean(stn_wake)
std_stn=std(stn_wake,1)
mean_stn-std_stn
wake_slices = reshape(wake,[slices/sum_depth,length(sample_list_wake)])'
nowake_slices = reshape(nowake,[slices/sum_depth,length(sample_list_nowake)])'
max_wake_slices_=sort(wake_slices')
max_nowake_slices_=sort(nowake_slices')
max_wake_slices=sum(max_wake_slices_)
max_nowake_slices=sum(max_nowake_slices_)
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
mean_wake=mean(max_wake_slices)
mean_nowake=mean(max_nowake_slices)
std_wake=std(max_wake_slices,1)
std_nowake=std(max_nowake_slices,1)
stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
stn_wake=(max_wake_slices-mean_nowake)/std_nowake
mean_stn=mean(stn_wake)
std_stn=std(stn_wake,1)
stn=mean_stn-std_stn
significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)

% %new
% 
% nowake=reshape(permute(anali_depth(1,1:length(sample_list_nowake),:,3),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali_depth(2,1:length(sample_list_wake),:,3),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_wake),:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices/sum_depth,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[slices/sum_depth,length(sample_list_nowake)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=sum(max_wake_slices_)
% max_nowake_slices=sum(max_nowake_slices_)
% % max_wake_slices=max_wake_slices_(end,:)
% % max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)

% 
% nowake=reshape(permute(anali_Canny_depth(1,1:length(sample_list_nowake),:,4,1),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali_Canny_depth(2,1:length(sample_list_wake),:,4,1),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_wake),:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices/sum_depth,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[slices/sum_depth,length(sample_list_nowake)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
% 
% nowake=reshape(permute(anali_Canny(1,1:length(sample_list_nowake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali_Canny(2,1:length(sample_list_wake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_wake),:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices,length(sample_list_wake)])'
% nowake_slices = reshape(nowake,[slices,length(sample_list_nowake)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=sum(max_wake_slices_(1:end,:))
% max_nowake_slices=sum(max_nowake_slices_(1:end,:))
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
% sum(max_wake_slices_>max(max_nowake_slices_(:)))
% sum(max_wake_slices_(:)>max(max_nowake_slices_(:)))
% 



end