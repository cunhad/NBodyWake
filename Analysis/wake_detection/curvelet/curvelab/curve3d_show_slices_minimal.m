function [ anali_sum2,sample_id_range_nw,sample_id_range_w ] = curve3d_show_slices_minimal(  )


%example:
% nowake(:)=anali(1,:,3,4);
%wake(:)=anali(2,:,3,4);

% addpath('../../processing');
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct3d'));

addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));

filename='_2dproj_z3_data_sl';
nc=1024;
new_nc=512;
trsh=20;
cut=1;
lev=1;
sigma = 5;
slices=32;
anal_lev=1;
size_mpc=4;
step_of_degree=1/8;
wavel_removal_factor=1/2;
sample_id_range_nw=[1:10];
sample_id_range_w=[1:10];
% sample_id_range_nw=[4,7];
% sample_id_range_w=[3,7];

% sample_id_range_nw=[2,7];
% sample_id_range_w=[2,7];


display_slice_nw = cell(1,length(sample_id_range_nw));
display_slice_w = cell(1,length(sample_id_range_w));

% display_slice_nw{find(sample_id_range_nw==4)}=[15];
% display_slice_nw{find(sample_id_range_nw==7)}=[13];
% display_slice_w{find(sample_id_range_w==3)}=[29];
% display_slice_w{find(sample_id_range_w==7)}=[28];

% display_slice_nw{find(sample_id_range_nw==2)}=[22];
% display_slice_nw{find(sample_id_range_nw==7)}=[28];
% display_slice_w{find(sample_id_range_w==2)}=[10];
% display_slice_w{find(sample_id_range_w==7)}=[29];

% display_slice_nw={[],[]};
% display_slice_w={[],[]};

% sample_id_range=[3, 7];
% sample_id_range=[1 : length(sample_list_nowake)];


specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw');
% sample_list_wake=sort_nat(sample_list_wake)


%
% F=zeros(nc);
% C_zero = fdct_wrapping(F,0);
F = ones(new_nc);
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

F=zeros(new_nc);
C_zero = fdct_wrapping(F,0);


% F2 = ones(new_nc,new_nc,slices/2);
F2 = ones(new_nc,new_nc,slices);
X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
C2 = fdct3d_forward(X2);
E2 = cell(size(C2));
for s=1:length(C2)
    E{s} = cell(size(C2{s}));
    for w=1:length(C2{s})
        A2 = C2{s}{w};
        E2{s}{w} = sqrt(sum(sum(A2.*conj(A2))) / prod(size(A2)));
    end
end



% F2=zeros(new_nc,new_nc,slices/2);
F2=zeros(new_nc,new_nc,slices);
C_zero2 = fdct3d_forward(F2);



% fig_curv=figure;
% fig_test1=figure;
% fig_test2=figure;
% fig_test3=figure;
% fig_test4=figure;

% fig_test1_m=figure;
% fig_test2_m=figure;
% fig_test3_m=figure;
% fig_test4_m=figure;

% fig_test1_max=figure;
% fig_test2_max=figure;
% fig_test3_max=figure;
% fig_test4_max=figure;

fig_test1_max2=figure;
fig_test2_max2=figure;
fig_test3_max2=figure;
fig_test4_max2=figure;
fig_test5_max2=figure;

% ax_curv=axes(fig_curv);
% ax_test1=axes(fig_test1);
% ax_test2=axes(fig_test2);
% ax_test3=axes(fig_test3);
% ax_test4=axes(fig_test4);
% 
% ax_test1_m=axes(fig_test1_m);
% ax_test2_m=axes(fig_test2_m);
% ax_test3_m=axes(fig_test3_m);
% ax_test4_m=axes(fig_test4_m);
% 
% ax_test1_max=axes(fig_test1_max);
% ax_test2_max=axes(fig_test2_max);
% ax_test3_max=axes(fig_test3_max);
% ax_test4_max=axes(fig_test4_max);

ax_test1_max2=axes(fig_test1_max2);
ax_test2_max2=axes(fig_test2_max2);
ax_test3_max2=axes(fig_test3_max2);
ax_test4_max2=axes(fig_test4_max2);
ax_test5_max2=axes(fig_test5_max2);

% sample_id_range=[1 : length(sample_list_nowake)];

for w_nw=1:2
% for w_nw=2
    
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
        ch='_7';
        coul='b';
        sample_id_range=sample_id_range_nw;
        display_slice=display_slice_nw;
    else
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
        ch='_4';
        coul='r';
        sample_id_range=sample_id_range_w;
        display_slice=display_slice_w;
    end
    
    
     for sample = sample_id_range
%     for sample = 1:length(sample_id_range)
%     for sample = 1:2
        
        map_3d_slices=zeros(new_nc,new_nc,slices);
        map_3d_slices_filt2d=zeros(new_nc,new_nc,slices);
        
        for slice_id=1:slices
            
            sample_id=(slices*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/',ch,filename,num2str(slice_id),'.bin')
            fid = fopen(filename_nowake);
            map = fread(fid,[nc nc], 'float32','l') ;
            fclose(fid);            
             map = imresize(map,new_nc/nc,'triangle');





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
            map(map<=1)=1;%to remove problem with holes
            map_3d_slices(:,:,slice_id)=log(map);
            %
            %             map(map<=1)=1;%to remove problem with holes
            %             map=map+1;
            %             test=log(log(map));
            
%             if false
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/new_nc)*[1:new_nc],(size_mpc/new_nc)*[1:new_nc],map_3d_slices(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('log for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            
            C = fdct_wrapping(log(map),0);
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero{s}{w};
                end
            end
            
            aux_count=1;
            for s = length(C)-lev:length(C)-1
                %                 thresh=0;
                thresh = sigma + sigma*(s == length(C));
                for w = 1:length(C{s})
                    %                     Ct{s}{w} = C{s}{w};
                    %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                    %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                    %                 Ct{s}{w} = C{s}{w};
                    curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                end
                curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                aux_count=aux_count+1;
            end
            
            
            
            
            
            BW2 = real(ifdct_wrapping(Ct,0));
            
            map_3d_slices_filt2d(:,:,slice_id) =real(ifdct_wrapping(Ct,0));
            
            
            %             if false
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/new_nc)*[1:new_nc],(size_mpc/new_nc)*[1:new_nc],map_3d_slices_filt2d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2 for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            
            anali(w_nw,sample,slice_id,1,:)=[max(BW2(:)),std(BW2(:)),max(BW2(:))/std(BW2(:)),kurtosis(kurtosis(BW2)),kurtosis(BW2(:))];
            
            % theta = 0:180/nc:180;
            theta = 0:step_of_degree:180;
            [R,xp] = radon(BW2,theta);
            
            unit=ones(new_nc);
            [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(R_u>new_nc*frac_cut)=R_nor(R_u>new_nc*frac_cut)./R_u(R_u>new_nc*frac_cut);
            R_nor(R_u<=new_nc*frac_cut)=0;
            
            boudary_removal_factor=2048/new_nc;
            
            n_levels=floor(log2(length(R_nor(:,1))));
            R_nor_filt=zeros(size(R_nor));
            for i=1:length(R_nor(1,:))
                %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
                [dc_dwt,levels] = wavedec(R(:,i),n_levels,'db1');
                D = wrcoef('d',dc_dwt,levels,'db1',lev);
                %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
                %                 D(1:floor((448+200)/boudary_removal_factor))=0;
                D(length(xp)-new_nc*wavel_removal_factor:end)=0;
                D(1:new_nc*wavel_removal_factor)=0;
                
                R_nor_filt(:,i)=D;
            end
            %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
            %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
            R_nor_filt(length(xp)-new_nc*wavel_removal_factor:end,:)=[];
            R_nor_filt(1:new_nc*wavel_removal_factor,:)=[];
            
            
            
            anali(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
            anali(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
            
%             
%             aa(:)=curv2(w_nw,sample,slice_id,:);
%             
%             a(:)=anali(w_nw,sample,slice_id,1,:);
%             b(:)=anali(w_nw,sample,slice_id,2,:);
%             c(:)=anali(w_nw,sample,slice_id,3,:);
%             d(:)=anali(w_nw,sample,slice_id,4,:);
%             
%             aa
%             a
%             b
%             c
%             d
%             
%             
%             test_curv_lv{sample_id}=   plot(ax_curv,aa,coul);
%             test1_lv{sample_id}=   plot(ax_test1,a,coul);
%             test2_lv{sample_id}=   plot(ax_test2,b,coul);
%             test3_lv{sample_id}=   plot(ax_test3,c,coul);
%             test4_lv{sample_id}=   plot(ax_test4,d,coul);
%             
%             clearvars a b c d aa
%             
%             %         clearvars R R_nor R_nor_filt kurt2;
%             
%             hold(ax_curv,'on');
%             hold(ax_test1,'on');
%             hold(ax_test2,'on');
%             hold(ax_test3,'on');
%             hold(ax_test4,'on');
            
        end
%         
%         for partition=0:1
%             
%             C_aux=map_3d_slices(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2);
%             
%             C = fdct3d_forward(C_aux);
%             Ct = C;
%             for s = 1:length(C)
%                 for w = 1:length(C{s})
%                     Ct{s}{w} = C_zero2{s}{w};
%                 end
%             end
%             
%             aux_count=1;
%             for s = 1:length(C)
%                 %                 thresh=0;
%                 thresh = sigma + sigma*(s == length(C));
%                 for w = 1:length(C{s})
%                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) >  0*E2{s}{w});
%                     %                     Ct{s}{w} = C{s}{w};
%                     %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
%                     %                 Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
%                     %                 Ct{s}{w} = C{s}{w};
%                     %                     curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
%                 end
%                 %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
%                 aux_count=aux_count+1;
%             end
%             
%             map_3d_slices_filt3d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2) = real(fdct3d_inverse(Ct));
%             
%         end
%         
%         for slice_id=1:slices
%             
%             this=map_3d_slices_filt3d(:,:,slice_id);
%             
%             theta = 0:180;
%             [R,xp] = radon(map_3d_slices_filt3d(:,:,slice_id),theta);
%             
%             unit=ones(nc);
%             [R_u,xp] = radon(unit,theta);
%             
%             frac_cut=0.5;
%             R_nor=R;
%             R_nor(R_u>nc*frac_cut)=R_nor(R_u>nc*frac_cut)./R_u(R_u>nc*frac_cut);
%             R_nor(R_u<=nc*frac_cut)=0;
%             
%             boudary_removal_factor=2048/nc;
%             
%             n_levels=floor(log2(length(R_nor(:,1))));
%             R_nor_filt=zeros(size(R_nor));
%             for i=1:length(R_nor(1,:))
%                 %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%                 [dc_dwt,levels] = wavedec(R(:,i),n_levels,'db1');
%                 D = wrcoef('d',dc_dwt,levels,'db1',lev);
%                 %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
%                 %                 D(1:floor((448+200)/boudary_removal_factor))=0;
%                 D(1237-256:end)=0;
%                 D(1:217+256)=0;
%                 
%                 R_nor_filt(:,i)=D;
%             end
%             %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
%             %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
%             R_nor_filt(1237-256:end,:)=[];
%             R_nor_filt(1:217+256,:)=[];
%             
%             anali_sum(w_nw,sample,slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
%             anali_sum(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
%             anali_sum(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
%             anali_sum(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
%             
%             a_max(:)=anali_sum(w_nw,sample,slice_id,1,:);
%             b_max(:)=anali_sum(w_nw,sample,slice_id,2,:);
%             c_max(:)=anali_sum(w_nw,sample,slice_id,3,:);
%             d_max(:)=anali_sum(w_nw,sample,slice_id,4,:);
%         end
%         
%         a(:)=max(anali(w_nw,sample,:,1,:),[],3);
%         b(:)=max(anali(w_nw,sample,:,2,:),[],3);
%         c(:)=max(anali(w_nw,sample,:,3,:),[],3);
%         d(:)=max(anali(w_nw,sample,:,4,:),[],3);
%         
%         a
%         b
%         c
%         d
%         
%         
%         test1m_lv{sample_id}=   plot(ax_test1_m,a,coul);
%         test2m_lv{sample_id}=   plot(ax_test2_m,b,coul);
%         test3m_lv{sample_id}=   plot(ax_test3_m,c,coul);
%         test4m_lv{sample_id}=   plot(ax_test4_m,d,coul);
%         
%         test1max_lv{sample_id}=   plot(ax_test1_max,a_max,coul);
%         test2max_lv{sample_id}=   plot(ax_test2_max,b_max,coul);
%         test3max_lv{sample_id}=   plot(ax_test3_max,c_max,coul);
%         test4max_lv{sample_id}=   plot(ax_test4_max,d_max,coul);
%         
        
        
        
%         for partition=0:1
        for partition=0
            
%             C_aux=map_3d_slices_filt2d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2);
            C_aux=map_3d_slices_filt2d(:,:,:);
            
            C = fdct3d_forward(C_aux);
            
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero2{s}{w};
                end
            end
            
            aux_count=1;
            for s = 2:3
                %                 thresh=0;
                thresh = sigma + sigma*(s == length(C));
                for w = 1:length(C{s})
                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E2{s}{w})
                    %                                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w};
                    %                   curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                end
                %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                aux_count=aux_count+1;
            end
            
%             map_3d_slices_filt2a3d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2) = real(fdct3d_inverse(Ct));
            map_3d_slices_filt2a3d(:,:,:) = real(fdct3d_inverse(Ct));            

        end
        
        for slice_id=1:slices
            
            
            this=map_3d_slices_filt2a3d(:,:,slice_id);

            
%             if false            
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/new_nc)*[1:new_nc],(size_mpc/new_nc)*[1:new_nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            
            theta = 0:step_of_degree:180;
            [R,xp] = radon(map_3d_slices_filt2a3d(:,:,slice_id),theta);
            
            unit=ones(new_nc);
            [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(R_u>new_nc*frac_cut)=R_nor(R_u>new_nc*frac_cut)./R_u(R_u>new_nc*frac_cut);
            R_nor(R_u<=new_nc*frac_cut)=0;
            
            boudary_removal_factor=2048/new_nc;
            
            n_levels=floor(log2(length(R_nor(:,1))));
            R_nor_filt=zeros(size(R_nor));
            for i=1:length(R_nor(1,:))
                [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
                 D = wrcoef('d',dc_dwt,levels,'db1',lev);
%                A = wrcoef('a',dc_dwt,levels,'db1',3);
                %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
                %                 D(1:floor((448+200)/boudary_removal_factor))=0;
%                 D(1237-256:end)=0;
%                 D(1:217+256)=0;
%                 A(length(xp)-new_nc*wavel_removal_factor:end)=0;
%                 A(1:new_nc*wavel_removal_factor)=0;
%                ridg_coeff=[];
 %               ridg_coeff=R(:,i)-A;
  %              ridg_coeff=diff(ridg_coeff);
   %             ridg_coeff(end+1,:)=0;
    %            ridg_coeff=(abs(ridg_coeff));
% 
%                 R_nor_filt(:,i)=ridg_coeff;
                 R_nor_filt(:,i)=D;
            end
            %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
            %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
%             R_nor_filt(1237-256:end,:)=[];
%             R_nor_filt(1:217+256,:)=[];

            
            
            R_nor_filt(length(xp)-new_nc*wavel_removal_factor:end,:)=[];
            R_nor_filt(1:new_nc*wavel_removal_factor,:)=[];
            
%             n_levels2=floor(log2(length(R_nor_filt(1,:))));
%             R_nor_filt2=zeros(size(R_nor_filt));
%             for i=1:length(R_nor_filt(:,1))
%                 %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
%                 [dc_dwt,levels] = wavedec(R_nor_filt(i,:),n_levels2,'db1');
%                 D = wrcoef('d',dc_dwt,levels,'db1',lev);
%                 %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
%                 %                 D(1:floor((448+200)/boudary_removal_factor))=0;
% %                 D(length(xp)-new_nc*wavel_removal_factor:end)=0;
% %                 D(1:new_nc*wavel_removal_factor)=0;
%                 
%                 R_nor_filt2(i,:)=D;
%             end
%             
            %2d wavelet filt
            
            
            orig_img=R_nor_filt;            
%             orig_img=R_nor_filt2;
%             
%             orig_img(length(xp)-new_nc*wavel_removal_factor:end,:)=[];
%             orig_img(1:new_nc*wavel_removal_factor,:)=[];
            
%             [c_2dwav,s_2dwav]=wavedec2(orig_img,2,'haar');
%             
% %             [H1,V1,D1] = detcoef2('all',c_2dwav,s_2dwav,1);
% %             A1 = appcoef2(c_2dwav,s_2dwav,'haar',1); 
% %             V1img = wcodemat(V1,255,'mat',1);
% %             H1img = wcodemat(H1,255,'mat',1);
% %             D1img = wcodemat(D1,255,'mat',1);
% %             A1img = wcodemat(A1,255,'mat',1);
% 
%             A1_2dwav = appcoef2(c_2dwav,s_2dwav,'haar',1); 
%             R_2dwav_filt=R_nor_filt-A1_2dwav;

            [cA1,cH1,cV1,cD1]=dwt2(orig_img,'haar');
            %cA1(:)=0;
            [cA2,cH2,cV2,cD2]=dwt2(cA1,'haar');
            cA2(:)=0;
            cA1=idwt2(cA2,cH2,cV2,cD2,'haar',size(cA1));
            R_2dwav_filt=idwt2(cA1,cH1,cV1,cD1,'haar',size(R_nor_filt));
            
            
%             if false            
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
%                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
                figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R_nor_filt);colorbar;
                xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
                figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R_2dwav_filt);colorbar;
                xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            
            anali_sum2(w_nw,sample,slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
            anali_sum2(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali_sum2(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
%             anali_sum2(w_nw,sample,slice_id,4,:)=[max(R_nor_filt2(:)),std(R_nor_filt2(:)),max(R_nor_filt2(:))/std(R_nor_filt2(:)),kurtosis(kurtosis(R_nor_filt2)),kurtosis(R_nor_filt2(:))];
            anali_sum2(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
            anali_sum2(w_nw,sample,slice_id,5,:)=[max(R_2dwav_filt(:)),std(R_2dwav_filt(:)),max(R_2dwav_filt(:))/std(R_2dwav_filt(:)),kurtosis(kurtosis(R_2dwav_filt)),kurtosis(R_2dwav_filt(:))];
%             a_max2(:)=anali_sum2(w_nw,sample,slice_id,1,:);
%             b_max2(:)=anali_sum2(w_nw,sample,slice_id,2,:);
%             c_max2(:)=anali_sum2(w_nw,sample,slice_id,3,:);
%             d_max2(:)=anali_sum2(w_nw,sample,slice_id,4,:);
%             
%             test1max2_lv{sample_id}=   plot(ax_test1_max2,a_max2,coul);
%             test2max2_lv{sample_id}=   plot(ax_test2_max2,b_max2,coul);
%             test3max2_lv{sample_id}=   plot(ax_test3_max2,c_max2,coul);
%             test4max2_lv{sample_id}=   plot(ax_test4_max2,d_max2,coul);
        end
        
        a_max2(:)=max(anali_sum2(w_nw,sample,:,1,:),[],3);
            b_max2(:)=max(anali_sum2(w_nw,sample,:,2,:),[],3);
            c_max2(:)=max(anali_sum2(w_nw,sample,:,3,:),[],3);
            d_max2(:)=max(anali_sum2(w_nw,sample,:,4,:),[],3);
            e_max2(:)=max(anali_sum2(w_nw,sample,:,5,:),[],3);
            
            test1max2_lv{sample_id}=   plot(ax_test1_max2,a_max2,coul);
            test2max2_lv{sample_id}=   plot(ax_test2_max2,b_max2,coul);
            test3max2_lv{sample_id}=   plot(ax_test3_max2,c_max2,coul);
            test4max2_lv{sample_id}=   plot(ax_test4_max2,d_max2,coul);
            test5max2_lv{sample_id}=   plot(ax_test5_max2,e_max2,coul);
        
        clearvars a b c d aa a_max b_max c_max d_max2 a_max2 b_max2 c_max2 d_max2 e_max2
        
        %         clearvars R R_nor R_nor_filt kurt2;
        
%         hold(ax_test1_m,'on');
%         hold(ax_test2_m,'on');
%         hold(ax_test3_m,'on');
%         hold(ax_test4_m,'on');
%         
%         hold(ax_test1_max,'on');
%         hold(ax_test2_max,'on');
%         hold(ax_test3_max,'on');
%         hold(ax_test4_max,'on');
%         
        hold(ax_test1_max2,'on');
        hold(ax_test2_max2,'on');
        hold(ax_test3_max2,'on');
        hold(ax_test4_max2,'on');
        hold(ax_test5_max2,'on');
        
    end
    
    
    
end
% 
% set(ax_curv, 'YScale', 'log');
% title(ax_curv,'curvelet kurtosis')
% set(ax_test1, 'YScale', 'log');
% title(ax_test1,'all map curvelet');
% set(ax_test2, 'YScale', 'log');
% title(ax_test2,'all ridgelet');
% set(ax_test3, 'YScale', 'log');
% title(ax_test3,'all ridgelet normalized');
% set(ax_test4, 'YScale', 'log');
% title(ax_test4,'all ridgelet normalized with wavelet');
% 
% 
% set(ax_test1_m, 'YScale', 'log');
% title(ax_test1_m,'map curvelet sum max');
% set(ax_test2_m, 'YScale', 'log');
% title(ax_test2_m,'ridgelet from max');
% set(ax_test3_m, 'YScale', 'log');
% title(ax_test3_m,'ridgelet normalized from max');
% set(ax_test4_m, 'YScale', 'log');
% title(ax_test4_m,'ridgelet normalized with wavelet from max');
% 
% 
% set(ax_test1_max, 'YScale', 'log');
% title(ax_test1_max,'map curvelet from curv');
% set(ax_test2_max, 'YScale', 'log');
% title(ax_test2_max,'ridgelet from curv');
% set(ax_test3_max, 'YScale', 'log');
% title(ax_test3_max,'ridgelet normalized from curv');
% set(ax_test4_max, 'YScale', 'log');
% title(ax_test4_max,'ridgelet normalized with wavelet from curv');

set(ax_test1_max2, 'YScale', 'log');
title(ax_test1_max2,'map curvelet from 2dcurv plus curv');
set(ax_test2_max2, 'YScale', 'log');
title(ax_test2_max2,'ridgelet from sum 2dcurv curv');
set(ax_test3_max2, 'YScale', 'log');
title(ax_test3_max2,'ridgelet normalized from 2dcurv plus curv');
set(ax_test4_max2, 'YScale', 'log');
title(ax_test4_max2,'ridgelet normalized with wavelet from 2dcurv plus curv');
set(ax_test5_max2, 'YScale', 'log');
title(ax_test5_max2,'ridgelet normalized with 2dwavelet from 2dcurv plus curv');
end




% [ anali_sum2,sample_id_range_nw,sample_id_range_w ] = curve3d_show_slices_minimal(  )
% 
% 
% nowake=reshape(permute(anali_sum2(1,sample_id_range_nw,:,4,1),[1,3,2,4,5]),[1,numel(anali_sum2(1,sample_id_range_nw,:,2,1))])
% wake=reshape(permute(anali_sum2(2,sample_id_range_w,:,4,1),[1,3,2,4,5]),[1,numel(anali_sum2(1,sample_id_range_w,:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[32,length(sample_id_range_w)])'
% nowake_slices = reshape(nowake,[32,length(sample_id_range_nw)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn