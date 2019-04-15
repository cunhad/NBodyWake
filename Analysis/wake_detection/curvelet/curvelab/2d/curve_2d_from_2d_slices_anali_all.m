function [  ] = curve_2d_from_2d_slices_anali_all(  )
% 
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali/';
% root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali_all/';
% root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_visual_all/';
% 

root='/home/asus/Dropbox/extras/storage/graham/ht/';
root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dclar_l2lr1na256_anali/';
% root_anali_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcr0_l2lr1ap256_anali/';
%root_anali_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_anali_all/';
%root_visual_2d='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2dcurv_s5lv2_visual_all/';

filename='3.000xv0.dat';
lenght_factor=1;
resol_factor=1;
pivot=[0,0,0];
rot_angle=[1.5708,0,0];
slices=32;
sum_depth=2;
% lev=2;
% sigma=5;
% step_of_degree=1;
% wavel_removal_factor=1/2;
% % snapshot=[];
% % snapshot=[13,28]*(128/32);
% snapshot=[9,29]*(128/32);
% visual_type=[1:2]; %if 1, shows the 2d proj; if 2 shows the ridgelet transformation
%  visual_in_or_out=[1,2]; %if 1 do visualization of the input, if 2 of the output

addpath('../../../../processing');
% addpath('../../../../preprocessing');



path_specs_in=strcat(root_anali_2d_in);

specs_nowake=dir(strcat(root_anali_2d_in,'/*nowake*'));
specs_nowake={specs_nowake.name};
specs_wake=dir(strcat(root_anali_2d_in,'/*wakeGmu*'));
specs_wake={specs_wake.name};

specs_list=dir(strcat(path_specs_in,'/*'));
specs_list={specs_list.name};
specs_list=sort_nat(specs_list);
specs_list=specs_list(3:end);

display(specs_list);

spec_nowake=specs_nowake{1};
path_samples_in=strcat(root_anali_2d_in,spec_nowake);
sample_list=dir(strcat(path_samples_in,'/sample*'));
sample_list=strcat('/',{sample_list.name},'/');
sample_list_nowake=sort_nat(sample_list)

spec_wake=specs_wake{1};
path_samples_in=strcat(root_anali_2d_in,spec_wake);
sample_list=dir(strcat(path_samples_in,'/sample*'));
sample_list_short=strcat('/',{sample_list.name},'/');
sample_list=strcat('/',{sample_list.name},'/half_lin_cutoff_half_tot_pert_nvpw/');
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
    
    fig1_curv_depth=figure;
    fig2_curv_depth=figure;
    fig3_curv_depth=figure;
    fig4_curv_depth=figure;
    fig5_curv_depth=figure;
    
    
    ax1_depth=axes(fig1_depth);
    ax2_depth=axes(fig2_depth);
    ax3_depth=axes(fig3_depth);
    ax4_depth=axes(fig4_depth);
    % ax5=axes(fig5);
    
    
    ax1_curv_depth=axes(fig1_curv_depth);
    ax2_curv_depth=axes(fig2_curv_depth);
    ax3_curv_depth=axes(fig3_curv_depth);
    ax4_curv_depth=axes(fig4_curv_depth);
    ax5_curv_depth=axes(fig5_curv_depth);
    
end

cd('../../../../preprocessing')
% [xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,sample_list_wake{1} );
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
        
        path_in=strcat(strcat(root_anali_2d_in,spec,char(sample_list(sample))),'anali/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')
        
        filename=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali.txt')
        filename_curv=strcat(path_in,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv.txt')
        
        info = dlmread(filename);
        info_curv = dlmread(filename_curv);
        lv_sz=prod(size(info_curv))/(slices*5);
        
        anali(w_nw,sample,:,:,:)=reshape(info,slices,4,5);
        anali_curv(w_nw,sample,:,:,:)=reshape(info_curv,slices,5,lv_sz);

        a(:)=max(anali(w_nw,sample,:,1,:),[],3);
        b(:)=max(anali(w_nw,sample,:,2,:),[],3);
        c(:)=max(anali(w_nw,sample,:,3,:),[],3);
        d(:)=max(anali(w_nw,sample,:,4,:),[],3);
        
        a_curv(:)=max(anali_curv(w_nw,sample,:,1,:),[],3);
        b_curv(:)=max(anali_curv(w_nw,sample,:,2,:),[],3);
        c_curv(:)=max(anali_curv(w_nw,sample,:,3,:),[],3);
        d_curv(:)=max(anali_curv(w_nw,sample,:,4,:),[],3);
        e_curv(:)=max(anali_curv(w_nw,sample,:,5,:),[],3);

        plot1{sample}=   plot(ax1,a,coul);
        plot2{sample}=   plot(ax2,b,coul);
        plot3{sample}=   plot(ax3,c,coul);
        plot4{sample}=   plot(ax4,d,coul);
        
        plot1_curv{sample}=   plot(ax1_curv,a_curv,coul);
        plot2_curv{sample}=   plot(ax2_curv,b_curv,coul);
        plot3_curv{sample}=   plot(ax3_curv,c_curv,coul);
        plot4_curv{sample}=   plot(ax4_curv,d_curv,coul);
        plot5_curv{sample}=   plot(ax5_curv,e_curv,coul);
        
        clearvars a b c d a_curv b_curv c_curv d_curv e_curv
        
        hold(ax1,'on');
        hold(ax2,'on');
        hold(ax3,'on');
        hold(ax4,'on');
        
        hold(ax1_curv,'on');
        hold(ax2_curv,'on');
        hold(ax3_curv,'on');
        hold(ax4_curv,'on');
        hold(ax5_curv,'on');
        
        if ~ismember(1,sum_depth)
            
            path_in_depth=strcat(strcat(root_anali_2d_in,spec,char(sample_list(sample))),'anali_depth_',num2str(sum_depth),'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2)),'-',num2str(rot_angle(3))),'ra','/','2dproj/dm/')
            
            filename_depth=strcat(path_in_depth,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_depth.txt')
            filename_curv_depth=strcat(path_in_depth,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv_depth.txt')
            
            info_depth = dlmread(filename_depth);
            info_curv_depth = dlmread(filename_curv_depth);
            lv_sz_depth=prod(size(info_curv_depth))*sum_depth/(slices*5);
            
            anali_depth(w_nw,sample,:,:,:)=reshape(info_depth,slices/sum_depth,4,5);
            anali_curv_depth(w_nw,sample,:,:,:)=reshape(info_curv_depth,slices/sum_depth,5,lv_sz_depth);
            
            a_depth(:)=max(anali_depth(w_nw,sample,:,1,:),[],3);
            b_depth(:)=max(anali_depth(w_nw,sample,:,2,:),[],3);
            c_depth(:)=max(anali_depth(w_nw,sample,:,3,:),[],3);
            d_depth(:)=max(anali_depth(w_nw,sample,:,4,:),[],3);
            
            a_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,1,:),[],3);
            b_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,2,:),[],3);
            c_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,3,:),[],3);
            d_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,4,:),[],3);
            e_curv_depth(:)=max(anali_curv_depth(w_nw,sample,:,5,:),[],3);
            
            plot1_depth{sample}=   plot(ax1_depth,a_depth,coul);
            plot2_depth{sample}=   plot(ax2_depth,b_depth,coul);
            plot3_depth{sample}=   plot(ax3_depth,c_depth,coul);
            plot4_depth{sample}=   plot(ax4_depth,d_depth,coul);
            
            plot1_curv_depth{sample}=   plot(ax1_curv_depth,a_curv_depth,coul);
            plot2_curv_depth{sample}=   plot(ax2_curv_depth,b_curv_depth,coul);
            plot3_curv_depth{sample}=   plot(ax3_curv_depth,c_curv_depth,coul);
            plot4_curv_depth{sample}=   plot(ax4_curv_depth,d_curv_depth,coul);
            plot5_curv_depth{sample}=   plot(ax5_curv_depth,e_curv_depth,coul);
            
            clearvars a_depth b_depth c_depth d_depth a_curv_depth b_curv_depth c_curv_depth d_curv_depth e_curv_depth
            
            hold(ax1_depth,'on');
            hold(ax2_depth,'on');
            hold(ax3_depth,'on');
            hold(ax4_depth,'on');
            
            hold(ax1_curv_depth,'on');
            hold(ax2_curv_depth,'on');
            hold(ax3_curv_depth,'on');
            hold(ax4_curv_depth,'on');
            hold(ax5_curv_depth,'on');
            
            
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
    
    set(ax1_curv_depth, 'YScale', 'log');
    title(ax1_curv_depth,strcat('average normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    set(ax2_curv_depth, 'YScale', 'log');
    title(ax2_curv_depth,strcat('std normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    set(ax3_curv_depth, 'YScale', 'log');
    title(ax3_curv_depth,strcat('skewness normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    set(ax4_curv_depth, 'YScale', 'log');
    title(ax4_curv_depth,strcat('kurtosis normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    set(ax5_curv_depth, 'YScale', 'log');
    title(ax5_curv_depth,strcat('4th moment normalised curvelet abs coef fast, slice= ',num2str(sum_depth)));
    
end

% nowake=reshape(permute(anali(1,1:length(sample_list_nowake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_nowake),:,2,1))])
% wake=reshape(permute(anali(2,1:length(sample_list_wake),:,4,1),[1,3,2,4,5]),[1,numel(anali(1,1:length(sample_list_wake),:,2,1))])
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

nowake=reshape(permute(anali_depth(1,1:length(sample_list_nowake),:,4,1),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_nowake),:,2,1))])
wake=reshape(permute(anali_depth(2,1:length(sample_list_wake),:,4,1),[1,3,2,4,5]),[1,numel(anali_depth(1,1:length(sample_list_wake),:,2,1))])
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
max_wake_slices=max_wake_slices_(end,:)
max_nowake_slices=max_nowake_slices_(end,:)
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

cd('../wake_detection/curvelet/curvelab/2d/')


end