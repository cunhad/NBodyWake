function [ anali,anali_curv_mom,anali3,anali3s,anali3_curv_mom,sample_id_range_nw,sample_id_range_w ] = curve3d_show_3D_minimal_nor_locangreala_CanRid_dev1(  )


%example:
% nowake(:)=anali(1,:,3,4);
%wake(:)=anali(2,:,3,4);

addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct3d'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));
addpath(genpath('/home/asus/Dropbox/Disrael/Doutorado/Research/NBodyWake/production/Analysis/wake_detection/ridgelet/ppft3_nomex/ppft3'))
run('/home/asus/Dropbox/Disrael/Doutorado/Research/NBodyWake/production/Analysis/wake_detection/ridgelet/ppft3_nomex/ppft3/initpath.m')


filename='_2dproj_z3_data_sl';
aux_path='/half_lin_cutoff_half_tot_pert_nvpw_v0p6'
nc=256;
rf=0.25;  %rf=nc/1024
lev_2d=1;
lev_2drig=1;
lev_3d=1;
lev_3drig=1;
slices=256;
size_mpc=4;
partition2d=2;
partition3rd=2;
step_of_degree=1*(180/256);
wavel_removal_factor=1/2;

% sample_id_range_nw=[1:10];
% sample_id_range_w=[1:10];
% sample_id_range_nw=[4,7];
% sample_id_range_w=[3,7];
% sample_id_range_nw=[4,7];
% sample_id_range_w=[4,7];
% 
% % sample_id_range=[3, 7];
% % sample_id_range=[1 : length(sample_list_nowake)];

sample_id_range_nw=[1];
sample_id_range_w=[8];

% sample_id_range_nw=[1:10];
% sample_id_range_w=[1:10];

display_slice_nw = cell(1,length(sample_id_range_nw));
display_slice_w = cell(1,length(sample_id_range_w));

% display_slice_nw{find(sample_id_range_nw==4)}=[15];
% display_slice_nw{find(sample_id_range_nw==7)}=[13];
% display_slice_w{find(sample_id_range_w==3)}=[29];
% display_slice_w{find(sample_id_range_w==7)}=[28];

% display_slice_nw{find(sample_id_range_nw==3)}=[22];
% display_slice_nw{find(sample_id_range_nw==4)}=[16];
% display_slice_w{find(sample_id_range_w==3)}=[22];
% display_slice_w{find(sample_id_range_w==4)}=[16];
% display_slice_nw{find(sample_id_range_nw==10)}=[5];
% display_slice_w{find(sample_id_range_w==10)}=[5];

% display_slice_nw{find(sample_id_range_nw==3)}=[22];
% display_slice_nw{find(sample_id_range_nw==4)}=[24];
% display_slice_w{find(sample_id_range_w==3)}=[22];
% display_slice_w{find(sample_id_range_w==4)}=[24];
% display_slice_nw{find(sample_id_range_nw==10)}=[5];
% display_slice_w{find(sample_id_range_w==10)}=[5];

% display_slice_nw{find(sample_id_range_nw==1)}=[8];
% display_slice_w{find(sample_id_range_w==1)}=[8];

% display_slice_nw={[],[]};
% display_slice_w={[],[]};


specs_path_list_nowake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'/4Mpc_2048c_1024p_zi63_nowakem')
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m')
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,aux_path);
% sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p55');
% sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_wrong');
% sample_list_wake=sort_nat(sample_list_wake)


%
% F=zeros(nc);
% C_zero = fdct_wrapping(F,0);
F = ones(nc);
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


F=zeros(nc);
C_zero = fdct_wrapping(F,0);

nc_anal3d=nc/partition2d;

% F2 = ones(nc_anal3d,nc_anal3d,slices/partition3rd);
% % F2 = ones(new_nc,new_nc,slices);
% X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
X2 = zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
X2(1+nc_anal3d/2,1+nc_anal3d/2,1+slices/(2*partition3rd))=(nc_anal3d*nc_anal3d*slices/partition3rd)^(1/3);
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
% F2=zeros(new_nc,new_nc,slices);
C_zero2 = fdct3d_forward(F2);


%radon

% nb_m=max(nb,slices)
% nb_m=min(nb,slices)
nb_m=sqrt(nb*slices);
Z2 = ones(nb_m,nb_m,nb_m);
Rad_norm = radon3(Z2);




% fig_curv=figure;
% fig_test1=figure;
fig_test2=figure;
fig_test3=figure;
fig_test4=figure;

% fig_test1_m=figure;
% fig_test2_m=figure;
% fig_test3_m=figure;
% fig_test4_m=figure;

% fig_test1_max=figure;
% fig_test2_max=figure;
% fig_test3_max=figure;
% fig_test4_max=figure;

% fig_test1_max2=figure;
% fig_test2_max2=figure;
% fig_test3_max2=figure;
% fig_test4_max2=figure;
% fig_test5_max2=figure;

% fig_test1_curv=figure;
% fig_test2_curv=figure;
% fig_test3_curv=figure;


fig_test1_curv_mom=figure;
fig_test2_curv_mom=figure;
fig_test3_curv_mom=figure;
fig_test4_curv_mom=figure;
fig_test5_curv_mom=figure;

fig3_test1=figure;
fig3_test2=figure;
fig3_test3=figure;
fig3_test4=figure;

fig3_test1_curv_mom=figure;
fig3_test2_curv_mom=figure;
fig3_test3_curv_mom=figure;
fig3_test4_curv_mom=figure;
fig3_test5_curv_mom=figure;

fig3_test1s=figure;
fig3_test2s=figure;
fig3_test3s=figure;
fig3_test4s=figure;

% ax_curv=axes(fig_curv);
% ax_test1=axes(fig_test1);
ax_test2=axes(fig_test2);
ax_test3=axes(fig_test3);
ax_test4=axes(fig_test4);

% ax_test1_m=axes(fig_test1_m);
% ax_test2_m=axes(fig_test2_m);
% ax_test3_m=axes(fig_test3_m);
% ax_test4_m=axes(fig_test4_m);
% 
% ax_test1_max=axes(fig_test1_max);
% ax_test2_max=axes(fig_test2_max);
% ax_test3_max=axes(fig_test3_max);
% ax_test4_max=axes(fig_test4_max);

% ax_test1_max2=axes(fig_test1_max2);
% ax_test2_max2=axes(fig_test2_max2);
% ax_test3_max2=axes(fig_test3_max2);
% ax_test4_max2=axes(fig_test4_max2);
% ax_test5_max2=axes(fig_test5_max2);

% ax_test1_curv=axes(fig_test1_curv);
% ax_test2_curv=axes(fig_test2_curv);
% ax_test3_curv=axes(fig_test3_curv);


ax_test1_curv_mom=axes(fig_test1_curv_mom);
ax_test2_curv_mom=axes(fig_test2_curv_mom);
ax_test3_curv_mom=axes(fig_test3_curv_mom);
ax_test4_curv_mom=axes(fig_test4_curv_mom);
ax_test5_curv_mom=axes(fig_test5_curv_mom);

ax3_test1=axes(fig3_test1);
ax3_test2=axes(fig3_test2);
ax3_test3=axes(fig3_test3);
ax3_test4=axes(fig3_test4);

ax3_test1s=axes(fig3_test1s);
ax3_test2s=axes(fig3_test2s);
ax3_test3s=axes(fig3_test3s);
ax3_test4s=axes(fig3_test4s);

ax3_test1_curv_mom=axes(fig3_test1_curv_mom);
ax3_test2_curv_mom=axes(fig3_test2_curv_mom);
ax3_test3_curv_mom=axes(fig3_test3_curv_mom);
ax3_test4_curv_mom=axes(fig3_test4_curv_mom);
ax3_test5_curv_mom=axes(fig3_test5_curv_mom);



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
%     for sample = 1

        
        map_3d_slices_pre=zeros(nc,nc,slices);
        map_3d_slices=zeros(nc,nc,slices);
        map_3d_slices_filt2d=zeros(nc,nc,slices);
        map_3d_slices_filt2a3d=zeros(nc,nc,slices);
        
        for slice_id=1:slices
            
            %             sample_id=(slices*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_',num2str(rf),'rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/',ch,filename,num2str(slice_id),'.bin')
            fid = fopen(filename_nowake);
            map = fread(fid,[nc nc], 'float32','l') ;
            fclose(fid);
            
            %             map(map<=1)=1;%to remove problem with holes
%                         map_3d_slices(:,:,slice_id)=log(map);

%                 map_3d_slices(:,:,slice_id)=log(map+1);
            map_3d_slices_pre(:,:,slice_id)=map;
            
            
        end
        
        for slice_id=1:slices

            
            sample_id=(slices*(sample-1))+slice_id;
            map=map_3d_slices_pre(:,:,slice_id);
            dc=(map-mean(map_3d_slices_pre(:)))/mean(map_3d_slices_pre(:));
%             dc=(map-mean(map(:)))/mean(map(:));
            %             map_3d_slices(:,:,slice_id)=(atan(dc+1));
            map_3d_slices(:,:,slice_id)=(atan((dc+1)*16));
            %             map_3d_slices(:,:,slice_id)=(atan(map));
            %             map_3d_slices(:,:,slice_id)=(atan(map+1));



%             if false
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],map_3d_slices(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('log for sample ',num2str(sample),' slice ',num2str(slice_id)));                
            end
            
%             C = fdct_wrapping(log(map),0);
%             Ct = C;
%             for s = 1:length(C)
%                 for w = 1:length(C{s})
%                     Ct{s}{w} = C_zero{s}{w};
%                 end
%             end
%             
%             aux_count=1;
%             for s = length(C)-lev:length(C)-1
%                 %                 thresh=0;
%                 thresh = sigma + sigma*(s == length(C));
%                 for w = 1:length(C{s})
%                     %                     Ct{s}{w} = C{s}{w};
%                     %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
%                     %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
%                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
%                     %                 Ct{s}{w} = C{s}{w};
% %                     curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
%                     
%                 end
% %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
% %                 curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}(:)));
%                 aux_count=aux_count+1;
%                 
%             end

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
            for s = length(C)-lev_2d:length(C)-1
                %                 thresh=0;
                %                 thresh = Sigma + Sigma*(s == length(C_));
                %                 collect=[];
                %                 collect_=[];
                
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
                    %                                         Ct{s}{w} = C_{s}{w};
                    %                                         Ct{s}{w} = C_{s}{w}.* ((C_{s}{w}) > thresh*E{s}{w});
                    %                                        Ct{s}{w} = C_{s}{w}.* ((C_{s}{w}) > 0*E{s}{w});
                    %                     Ct{s}{w} = C_{s}{w}.* ((C_{s}{w}) > -thresh*E{s}{w}&(C_{s}{w}) < thresh*E{s}{w});
                    %                                     Ct{s}{w} = C_{s}{w};
                    %                     curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                    %                     collect=[ collect ;abs(Ct{s}{w}(:))];
                    %                     collect_=[ collect_ ;abs(Ct{s}{w}(:))/E{s}{w}];
                    
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
                %                 curv(w_nw,sample,slice_id,aux_count)=kurtosis(collect);
                %
                %                   curv(aux_count)=kurtosis(collect);
                %                   curv_(aux_count)=kurtosis(collect_);
                %                   curv_m4(aux_count)=kurtosis(collect)*std(collect)^4;
                %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                %                 curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}(:)));
                
                
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
                
                curv_mom1(aux_count)=AT;
                curv_mom2(aux_count)=var_t;
                curv_mom3(aux_count)=skew_t;
                curv_mom4(aux_count)=kurt_t;
                curv_mom5(aux_count)=rho_t;

                aux_count=aux_count+1;
                
                
                
            end
          
            map_3d_slices_filt2d(:,:,slice_id) =real(ifdct_wrapping(Ct2,0));
            
            
            %             if false
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],map_3d_slices_filt2d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2 for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            

            theta = 0:step_of_degree:180;
            [R,xp] = radon(map_3d_slices_filt2d(:,:,slice_id),theta);
            
            unit=ones(nc);
            [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(R_u>nc*frac_cut)=R_nor(R_u>nc*frac_cut)./R_u(R_u>nc*frac_cut);
            R_nor(R_u<=nc*frac_cut)=0;
            
            
            
            
            boudary_removal_factor=2048/nc;
            
            n_levels=floor(log2(length(R_nor(:,1))));
            R_nor_filt=zeros(size(R_nor));
            for i=1:length(R_nor(1,:))
                %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
                [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
                D = wrcoef('d',dc_dwt,levels,'db1',lev_2drig);
                %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
                %                 D(1:floor((448+200)/boudary_removal_factor))=0;
                D(length(xp)-nc*wavel_removal_factor:end)=0;
                D(1:nc*wavel_removal_factor)=0;
                
                R_nor_filt(:,i)=D;
            end
            %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
            %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
            R_nor_filt(length(xp)-nc*wavel_removal_factor:end,:)=[];
            R_nor_filt(1:nc*wavel_removal_factor,:)=[];
            
% 
%             if ismember(slice_id,display_slice{find(sample_id_range==sample)})
%                 
%                 figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R);colorbar;
%                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
%                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 set(gca,'FontName','FixedWidth');
%                 set(gca,'FontSize',16);
%                 set(gca,'linewidth',2);
%                 title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
%                 
%                 figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R_nor);colorbar;
%                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
%                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 set(gca,'FontName','FixedWidth');
%                 set(gca,'FontSize',16);
%                 set(gca,'linewidth',2);
%                 title(strcat('wfilt norm filt2 for sample ',num2str(sample),' slice ',num2str(slice_id)));
%                 
%                 
%                 

            
%             anali(w_nw,sample,slice_id,1,:)=[max(BW2(:)),std(BW2(:)),max(BW2(:))/std(BW2(:)),kurtosis(kurtosis(BW2)),kurtosis(BW2(:))];
%             anali(w_nw,sample,slice_id,1,:)=curv;
            anali(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
            anali(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
%             anali(w_nw,sample,slice_id,5,:)=curv_;
  
%             anali_curv(w_nw,sample,slice_id,1,:)=curv;
%             anali_curv(w_nw,sample,slice_id,2,:)=curv_;
%             anali_curv(w_nw,sample,slice_id,3,:)=curv_m4;

            
            anali_curv_mom(w_nw,sample,slice_id,1,:)=curv_mom1;
            anali_curv_mom(w_nw,sample,slice_id,2,:)=curv_mom2;
            anali_curv_mom(w_nw,sample,slice_id,3,:)=curv_mom3;
            anali_curv_mom(w_nw,sample,slice_id,4,:)=curv_mom4;
            anali_curv_mom(w_nw,sample,slice_id,5,:)=curv_mom5;


%             aa(:)=curv2(w_nw,sample,slice_id,:);
%             
%             a(:)=anali(w_nw,sample,slice_id,1,:);
%             b(:)=anali(w_nw,sample,slice_id,2,:);
%             c(:)=anali(w_nw,sample,slice_id,3,:);
%             d(:)=anali(w_nw,sample,slice_id,4,:);
%             
%             a_curv(:)=anali_curv(w_nw,sample,slice_id,1,:);
%             b_curv(:)=anali_curv(w_nw,sample,slice_id,2,:);
      
%            
% %             test1_lv{sample_id}=   plot(ax_test1,a,coul);
%             test2_lv{sample_id}=   plot(ax_test2,b,coul);
%             test3_lv{sample_id}=   plot(ax_test3,c,coul);
%             test4_lv{sample_id}=   plot(ax_test4,d,coul);
%             
%             test1_curv{sample_id}=   plot(ax_test1_curv,a_curv,coul);
%             test2_curv{sample_id}=   plot(ax_test2_curv,b_curv,coul);
%             
%             clearvars a b c d aa a_curv b_curv
% %             
% %             %         clearvars R R_nor R_nor_filt kurt2;
% %             
% %             hold(ax_curv,'on');
% %             hold(ax_test1,'on');
%             hold(ax_test2,'on');
%             hold(ax_test3,'on');
%             hold(ax_test4,'on');
% 
%             hold(ax_test1_curv,'on');
%             hold(ax_test2_curv,'on');

            
        end
        
        
            b(:)=max(anali(w_nw,sample,:,2,:),[],3);
            c(:)=max(anali(w_nw,sample,:,3,:),[],3);
            d(:)=max(anali(w_nw,sample,:,4,:),[],3);
            
%             a_curv(:)=max(anali_curv(w_nw,sample,:,1,:),[],3);
%             b_curv(:)=max(anali_curv(w_nw,sample,:,2,:),[],3);
%             c_curv(:)=max(anali_curv(w_nw,sample,:,3,:),[],3);
            
            a_curv_mom(:)=max(anali_curv_mom(w_nw,sample,:,1,:),[],3);
            b_curv_mom(:)=max(anali_curv_mom(w_nw,sample,:,2,:),[],3);
            c_curv_mom(:)=max(anali_curv_mom(w_nw,sample,:,3,:),[],3);
            d_curv_mom(:)=max(anali_curv_mom(w_nw,sample,:,4,:),[],3);
            e_curv_mom(:)=max(anali_curv_mom(w_nw,sample,:,5,:),[],3);
       
%             test1_lv{sample_id}=   plot(ax_test1,a,coul);
            test2_lv{sample_id}=   plot(ax_test2,b,coul);
            test3_lv{sample_id}=   plot(ax_test3,c,coul);
            test4_lv{sample_id}=   plot(ax_test4,d,coul);
            
            %             test1_curv{sample_id}=   plot(ax_test1_curv,a_curv,coul);
            %             test2_curv{sample_id}=   plot(ax_test2_curv,b_curv,coul);
            %             test3_curv{sample_id}=   plot(ax_test3_curv,c_curv,coul);
            
            if (size(a_curv_mom)==1)
                
                test1_curv_mom{sample_id}=   scatter(ax_test1_curv_mom,1,a_curv_mom,coul);
                test2_curv_mom{sample_id}=   scatter(ax_test2_curv_mom,1,b_curv_mom,coul);
                test3_curv_mom{sample_id}=   scatter(ax_test3_curv_mom,1,c_curv_mom,coul);
                test4_curv_mom{sample_id}=   scatter(ax_test4_curv_mom,1,d_curv_mom,coul);
                test5_curv_mom{sample_id}=   scatter(ax_test5_curv_mom,1,e_curv_mom,coul);
                
            else
                
                test1_curv_mom{sample_id}=   plot(ax_test1_curv_mom,a_curv_mom,coul);
                test2_curv_mom{sample_id}=   plot(ax_test2_curv_mom,b_curv_mom,coul);
                test3_curv_mom{sample_id}=   plot(ax_test3_curv_mom,c_curv_mom,coul);
                test4_curv_mom{sample_id}=   plot(ax_test4_curv_mom,d_curv_mom,coul);
                test5_curv_mom{sample_id}=   plot(ax_test5_curv_mom,e_curv_mom,coul);
                
            end
            clearvars b c d a_curv b_curv c_curv a_curv_mom b_curv_mom c_curv_mom d_curv_mom e_curv_mom
%             

%             hold(ax_test1,'on');
            hold(ax_test2,'on');
            hold(ax_test3,'on');
            hold(ax_test4,'on');
% 
%             hold(ax_test1_curv,'on');
%             hold(ax_test2_curv,'on');
%             hold(ax_test3_curv,'on');
%             
            hold(ax_test1_curv_mom,'on');
            hold(ax_test2_curv_mom,'on');
            hold(ax_test3_curv_mom,'on');
            hold(ax_test4_curv_mom,'on');            
            hold(ax_test5_curv_mom,'on');            
        
        
            
            
            
            
            for partition=0:partition3rd-1
                for partition_x=0:partition2d-1
                    for partition_y=0:partition2d-1
                        
                        
                        %             C_aux=map_3d_slices_filt2d(:,:,:);
                        %             C_aux=map_3d_slices((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
                        Caux=map_3d_slices_filt2d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
                        
                        
                        C = fdct3d_forward(Caux);
                        
                        %                         Ct = C;
                        %                         for s = 1:length(C)
                        %                             for w = 1:length(C{s})
                        %                                 Ct{s}{w} = C_zero2{s}{w};
                        %                             end
                        %                         end
                        
                        Ct = C;
                        Ct2=C;
                        for s = 1:length(C)
                            for w = 1:length(C{s})
                                Ct{s}{w} = C_zero2{s}{w};
                                Ct2{s}{w} = C_zero2{s}{w};
                            end
                        end
                        
                        aux_count=1;
                        for s = length(C)-lev_3d:length(C)-1
                            %                 thresh=0;
                            %                             thresh = Sigma + Sigma*(s == length(C));
                            
                            As=0;
                            Bs=0;
                            Cs=0;
                            Ds=0;
                            
                            size_nt=0;
                            
                            
                            %                             for w = 1:length(C{s})
                            %                                 %                                         Ct{s}{w} = C_{s}{w};
                            %                                 %                                         Ct{s}{w} = C_{s}{w}.* ((C_{s}{w}) > thresh*E2{s}{w});
                            %                                 Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*C{s}{w});
                            %                                 %                     Ct{s}{w} = C_{s}{w}.* ((C_{s}{w}) > -thresh*E2{s}{w}&(C_{s}{w}) < thresh*E{2s}{w});
                            %                                 %                                     Ct{s}{w} = C{s}{w};
                            %                                 %                     curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                            %                                 %                     collect=[ collect ;abs(Ct{s}{w}(:))];
                            %                                 %                     collect_=[ collect_ ;abs(Ct{s}{w}(:))/E{s}{w}];
                            %
                            %                                 ave=mean(abs(Ct{s}{w}(:))/E2{s}{w});
                            %                                 sig=std(abs(Ct{s}{w}(:))/E2{s}{w})^2;
                            %                                 del=skewness(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(3/2));
                            %                                 rho=kurtosis(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(2));
                            %
                            %                                 A_=ave;
                            %                                 B_=sig+A_^2;
                            %                                 C_=del+3*B_*A_-2*A_^3;
                            %                                 D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
                            %
                            %                                 size_n=prod(size(Ct{s}{w}(:)));
                            %                                 size_nt=size_nt+size_n;
                            %
                            %                                 As=A_*size_n+As;
                            %                                 Bs=B_*size_n+Bs;
                            %                                 Cs=C_*size_n+Cs;
                            %                                 Ds=D_*size_n+Ds;
                            %
                            %
                            %
                            %                             end
                            %                             %                 curv(w_nw,sample,slice_id,aux_count)=kurtosis(collect);
                            %
                            %                             %                   curv(aux_count)=kurtosis(collect);
                            %                             %                   curv_(aux_count)=kurtosis(collect_);
                            %                             %                   curv_m4(aux_count)=kurtosis(collect)*std(collect)^4;
                            %                             % %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                            %                             % %                 curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}(:)));
                            
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
                            
                            curv3_mom1(aux_count)=AT;
                            curv3_mom2(aux_count)=var_t;
                            curv3_mom3(aux_count)=skew_t;
                            curv3_mom4(aux_count)=kurt_t;
                            curv3_mom5(aux_count)=rho_t;
                            
                            aux_count=aux_count+1;
                            
                            
                            
                        end
                        
                        %                         map_3d_slices_filt2a3d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct));
                        
                        map_3d_slices_filt2a3d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct2));
%                         map_3d_slices_filt2a3d((partition_x)*nc/partition2d+1:(partition_x)*nc/partition2d+nc/partition2d,(partition_y)*nc/partition2d+1:(partition_y)*nc/partition2d+nc/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = edge(real(fdct3d_inverse(Ct2)),'Canny',0.5);
                        
%                         [BW,thresOut] = edge(this,'Canny',0.5);
                        
                        anali3_curv_mom(w_nw,sample,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,5,:)=curv3_mom5;                        
                        anali3_curv_mom(w_nw,sample,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,1,:)=curv3_mom1;
                        anali3_curv_mom(w_nw,sample,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,2,:)=curv3_mom2;
                        anali3_curv_mom(w_nw,sample,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,3,:)=curv3_mom3;
                        anali3_curv_mom(w_nw,sample,1+partition_x+partition_y*partition2d+partition*partition2d*partition3rd,4,:)=curv3_mom4;
                        
                    end
                end
            end
        
%          anali3_curv_mom(w_nw,sample,slice_id,1,:)=curv3_mom1;
%             anali3_curv_mom(w_nw,sample,slice_id,2,:)=curv3_mom2;
%             anali3_curv_mom(w_nw,sample,slice_id,3,:)=curv3_mom3;
%             anali3_curv_mom(w_nw,sample,slice_id,4,:)=curv3_mom4;
%             anali3_curv_mom(w_nw,sample,slice_id,5,:)=curv3_mom5;
            
        
        for slice_id=1:slices
            
            map_3d_slices_filt2a3d(:,:,slice_id)= edge(map_3d_slices_filt2a3d(:,:,slice_id),'Canny',0.5);
            this=map_3d_slices_filt2a3d(:,:,slice_id);

            BW=this;
%             BW = edge(this);
%             [BW,thresOut] = edge(this,'Canny',0.5);
%             if false            
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],BW); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
            end
            
            theta = 0:step_of_degree:180;
%             [R,xp] = radon(map_3d_slices_filt2a3d(:,:,slice_id),theta);
            [R,xp] = radon(BW,theta);

            

            unit=ones(nc);
            [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(R_u>nc*frac_cut)=R_nor(R_u>nc*frac_cut)./R_u(R_u>nc*frac_cut);
            R_nor(R_u<=nc*frac_cut)=0;
            
            boudary_removal_factor=2048/nc;
            
            n_levels=floor(log2(length(R_nor(:,1))));
            R_nor_filt=zeros(size(R_nor));
            for i=1:length(R_nor(1,:))
                [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
                 D = wrcoef('d',dc_dwt,levels,'db1',lev_2drig);
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

            
            
            R_nor_filt(length(xp)-nc*wavel_removal_factor:end,:)=[];
            R_nor_filt(1:nc*wavel_removal_factor,:)=[];
            
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
%             %2d wavelet filt
%             
%             
%             orig_img=R_nor_filt;            
% %             orig_img=R_nor_filt2;
% %             
% %             orig_img(length(xp)-new_nc*wavel_removal_factor:end,:)=[];
% %             orig_img(1:new_nc*wavel_removal_factor,:)=[];
%             
% %             [c_2dwav,s_2dwav]=wavedec2(orig_img,2,'haar');
% %             
% % %             [H1,V1,D1] = detcoef2('all',c_2dwav,s_2dwav,1);
% % %             A1 = appcoef2(c_2dwav,s_2dwav,'haar',1); 
% % %             V1img = wcodemat(V1,255,'mat',1);
% % %             H1img = wcodemat(H1,255,'mat',1);
% % %             D1img = wcodemat(D1,255,'mat',1);
% % %             A1img = wcodemat(A1,255,'mat',1);
% % 
% %             A1_2dwav = appcoef2(c_2dwav,s_2dwav,'haar',1); 
% %             R_2dwav_filt=R_nor_filt-A1_2dwav;
% 
%             [cA1,cH1,cV1,cD1]=dwt2(orig_img,'haar');
%             %cA1(:)=0;
%             [cA2,cH2,cV2,cD2]=dwt2(cA1,'haar');
%             cA2(:)=0;
%             cA1=idwt2(cA2,cH2,cV2,cD2,'haar',size(cA1));
%             R_2dwav_filt=idwt2(cA1,cH1,cV1,cD1,'haar',size(R_nor_filt));
            
            
%             if false            
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
%                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
                figure; imagesc(0:step_of_degree:180,(size_mpc/nc)*[1:nc],R_nor_filt);colorbar;
                xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
%                 figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R_2dwav_filt);colorbar;
%                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
%                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 set(gca,'FontName','FixedWidth');
%                 set(gca,'FontSize',16);
%                 set(gca,'linewidth',2);
%                 title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
            end
            
            anali3(w_nw,sample,slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
            anali3(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali3(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
%             anali_sum2(w_nw,sample,slice_id,4,:)=[max(R_nor_filt2(:)),std(R_nor_filt2(:)),max(R_nor_filt2(:))/std(R_nor_filt2(:)),kurtosis(kurtosis(R_nor_filt2)),kurtosis(R_nor_filt2(:))];
            anali3(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
%             anali_sum2(w_nw,sample,slice_id,5,:)=[max(R_2dwav_filt(:)),std(R_2dwav_filt(:)),max(R_2dwav_filt(:))/std(R_2dwav_filt(:)),kurtosis(kurtosis(R_2dwav_filt)),kurtosis(R_2dwav_filt(:))];
%             a_max2(:)=anali_sum2(w_nw,sample,slice_id,1,:);
%             b_max2(:)=anali_sum2(w_nw,sample,slice_id,2,:);
%             c_max2(:)=anali_sum2(w_nw,sample,slice_id,3,:);
%             d_max2(:)=anali_sum2(w_nw,sample,slice_id,4,:);
%             
%             test1max2_lv{sample_id}=   plot(ax_test1_max2,a_max2,coul);
%             test2max2_lv{sample_id}=   plot(ax_test2_max2,b_max2,coul);
%             test3max2_lv{sample_id}=   plot(ax_test3_max2,c_max2,coul);
%             test4max2_lv{sample_id}=   plot(ax_test4_max2,d_max2,coul);


%             anali3_curv_mom(w_nw,sample,slice_id,1,:)=curv3_mom1;
%             anali3_curv_mom(w_nw,sample,slice_id,2,:)=curv3_mom2;
%             anali3_curv_mom(w_nw,sample,slice_id,3,:)=curv3_mom3;
%             anali3_curv_mom(w_nw,sample,slice_id,4,:)=curv3_mom4;
%             anali3_curv_mom(w_nw,sample,slice_id,5,:)=curv3_mom5;
%             
            
        end
        
            a3(:)=max(anali3(w_nw,sample,:,1,:),[],3);
            b3(:)=max(anali3(w_nw,sample,:,2,:),[],3);
            c3(:)=max(anali3(w_nw,sample,:,3,:),[],3);
            d3(:)=max(anali3(w_nw,sample,:,4,:),[],3);
%             e_max2(:)=max(anali_sum2(w_nw,sample,:,5,:),[],3);

            
            

            a3_curv_mom(:)=max(anali3_curv_mom(w_nw,sample,:,1,:),[],[3]);
            b3_curv_mom(:)=max(anali3_curv_mom(w_nw,sample,:,2,:),[],[3]);
            c3_curv_mom(:)=max(anali3_curv_mom(w_nw,sample,:,3,:),[],[3]);
            d3_curv_mom(:)=max(anali3_curv_mom(w_nw,sample,:,4,:),[],[3]);
            e3_curv_mom(:)=max(anali3_curv_mom(w_nw,sample,:,5,:),[],[3]);

            
            test1max3_lv{sample_id}=   plot(ax3_test1,a3,coul);
            test2max3_lv{sample_id}=   plot(ax3_test2,b3,coul);
            test3max3_lv{sample_id}=   plot(ax3_test3,c3,coul);
            test4max3_lv{sample_id}=   plot(ax3_test4,d3,coul);
            %             test5max2_lv{sample_id}=   plot(ax_test5_max2,e_max2,coul);
            
            if (size(a3_curv_mom)==1)
                
                test1_curv3_mom{sample_id}=   scatter(ax3_test1_curv_mom,1,a3_curv_mom,coul);
                test2_curv3_mom{sample_id}=   scatter(ax3_test2_curv_mom,1,b3_curv_mom,coul);
                test3_curv3_mom{sample_id}=   scatter(ax3_test3_curv_mom,1,c3_curv_mom,coul);
                test4_curv3_mom{sample_id}=   scatter(ax3_test4_curv_mom,1,d3_curv_mom,coul);
                test5_curv3_mom{sample_id}=   scatter(ax3_test5_curv_mom,1,e3_curv_mom,coul);
                
            else
                
                test1_curv3_mom{sample_id}=   plot(ax3_test1_curv_mom,a3_curv_mom,coul);
                test2_curv3_mom{sample_id}=   plot(ax3_test2_curv_mom,b3_curv_mom,coul);
                test3_curv3_mom{sample_id}=   plot(ax3_test3_curv_mom,c3_curv_mom,coul);
                test4_curv3_mom{sample_id}=   plot(ax3_test4_curv_mom,d3_curv_mom,coul);
                test5_curv3_mom{sample_id}=   plot(ax3_test5_curv_mom,e3_curv_mom,coul);
                
            end
            
         clearvars a3 b3 c3 d3 a3_curv_mom b3_curv_mom c3_curv_mom d3_curv_mom e3_curv_mom
%             clearvars a b c d aa a_max b_max c_max
        
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
        hold(ax3_test1,'on');
        hold(ax3_test2,'on');
        hold(ax3_test3,'on');
        hold(ax3_test4,'on');
%         hold(ax_test5_max2,'on');
%         
%         hold(ax_test1_max2,'on');
%         hold(ax_test2_max2,'on');

            hold(ax3_test1_curv_mom,'on');
            hold(ax3_test2_curv_mom,'on');
            hold(ax3_test3_curv_mom,'on');
            hold(ax3_test4_curv_mom,'on');            
            hold(ax3_test5_curv_mom,'on');       
            
            
            
            
            
            

            
            tst=map_3d_slices_filt2a3d;
%             tst(tst<0)=0;
% %             tst=abs(tst);
            this=sum(tst,3);
% 
            BW=this;
%             BW = edge(this);
%             [BW,thresOut] = edge(this,'Canny',0.5);
%             if false            
            if ismember(slice_id,display_slice{find(sample_id_range==sample)})
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
                figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/nc)*[1:nc],BW); colorbar; axis('image');
                xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
                set(gca,'FontName','FixedWidth');
                set(gca,'FontSize',16);
                set(gca,'linewidth',2);
                title(strcat('filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
                
            end
            
%             theta = 0:step_of_degree:180;
%             [R,xp] = radon(map_3d_slices_filt2a3d(:,:,slice_id),theta);
%             [R,xp] = radon(BW,theta);
        
            

            R = radon3(                                                                                                                                                                                                                                                                                                                                         (map_3d_slices_filt2a3d,[nb_m nb_m nb_m]));
%             
%             R_nor_=R./Rad_norm;
%             R_nor_(isinf(R_nor)|isnan(R_nor))=0;
%             

%             unit=ones(nc);
%             [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(Rad_norm>nb_m*nb_m*frac_cut)=R_nor(Rad_norm>nb_m*nb_m*frac_cut)./Rad_norm(Rad_norm>nb_m*nb_m*frac_cut);
            R_nor(Rad_norm<=nb_m*nb_m*frac_cut)=0;
            
%             R_nor=R_nor_(:,1+end/4:-1+3*end/4,1+end/4:-1+3*end/4,1+end/4:-1+3*end/4);

            n_levels=floor(log2(max(size(squeeze(R_nor(1,:,:,:))))));
            
            WDEC = wavedec3(squeeze(R_nor(1,:,:,:)),n_levels,'db1') ;
            R_nor_filt_x = waverec3(WDEC,'d',lev_2drig);
            WDEC = wavedec3(squeeze(R_nor(2,:,:,:)),1,'db1') ;
            R_nor_filt_y = waverec3(WDEC,'d',lev_2drig);
            WDEC = wavedec3(squeeze(R_nor(3,:,:,:)),1,'db1') ;
            R_nor_filt_z = waverec3(WDEC,'d',lev_2drig);
            
            R_nor_filt=[R_nor_filt_x;R_nor_filt_y;R_nor_filt_z];
            
% %             boudary_removal_factor=2048/nc;
% %             
% % %             n_levels=floor(log2(length(R_nor(1,:,1,1))));
% % %             R_nor_filt=zeros(size(R_nor));
% % %             for i=1:length(R_nor(1,:))
% % %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
% % %                  D = wrcoef('d',dc_dwt,levels,'db1',lev_2drig);
% % %                  R_nor_filt(:,i)=D;
% % %             end
% %                 
% %             R_nor_filt(length(xp)-nc*wavel_removal_factor:end,:)=[];
% %             R_nor_filt(1:nc*wavel_removal_factor,:)=[];
% 
% %             if false            
%             if ismember(slice_id,display_slice{find(sample_id_range==sample)})
%                 
% %                 figure; imagesc((size_mpc/nc)*[1:nc],(size_mpc/1024)*[1:nc],map_3d_slices_filt2a3d(:,:,slice_id)); colorbar; axis('image');
%                 figure; imagesc(0:step_of_degree:180,(size_mpc/nc)*[1:nc],R_nor_filt);colorbar;
%                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
%                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 set(gca,'FontName','FixedWidth');
%                 set(gca,'FontSize',16);
%                 set(gca,'linewidth',2);
%                 title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
%                 
% %                 figure; imagesc(0:step_of_degree:180,(size_mpc/new_nc)*[1:new_nc],R_2dwav_filt);colorbar;
% %                 xlabel('$\theta (degrees)$', 'interpreter', 'latex', 'fontsize', 20);
% %                 ylabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
% %                 set(gca,'FontName','FixedWidth');
% %                 set(gca,'FontSize',16);
% %                 set(gca,'linewidth',2);
% %                 title(strcat('wfilt ridg filt2a3 for sample ',num2str(sample),' slice ',num2str(slice_id)));
%             end
            
            anali3s(w_nw,sample,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(this(:))];
            anali3s(w_nw,sample,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(R(:))];
            anali3s(w_nw,sample,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(R_nor(:))];
            anali3s(w_nw,sample,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(R_nor_filt(:))];

%             
            
        
        
            a3s(:)=(anali3s(w_nw,sample,1,:));
            b3s(:)=(anali3s(w_nw,sample,2,:));
            c3s(:)=(anali3s(w_nw,sample,3,:));
            d3s(:)=(anali3s(w_nw,sample,4,:));
            
            test1max3_lvs{sample_id}=   plot(ax3_test1s,a3s,coul);
            test2max3_lvs{sample_id}=   plot(ax3_test2s,b3s,coul);
            test3max3_lvs{sample_id}=   plot(ax3_test3s,c3s,coul);
            test4max3_lvs{sample_id}=   plot(ax3_test4s,d3s,coul);
            %             test5max2_lv{sample_id}=   plot(ax_test5_max2,e_max2,coul);
            
           
            
         clearvars a3s b3s c3s d3s
         
        hold(ax3_test1s,'on');
        hold(ax3_test2s,'on');
        hold(ax3_test3s,'on');
        hold(ax3_test4s,'on');
     




    end
    
    
    
end
% 
% set(ax_curv, 'YScale', 'log');
% title(ax_curv,'curvelet kurtosis')
% set(ax_test1, 'YScale', 'log');
% title(ax_test1,'all map curvelet');
set(ax_test2, 'YScale', 'log');
title(ax_test2,'ridgelet');
set(ax_test3, 'YScale', 'log');
title(ax_test3,'ridgelet normalized');
set(ax_test4, 'YScale', 'log');
title(ax_test4,'ridgelet normalized with wavelet');

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
% 
% set(ax_test1_max2, 'YScale', 'log');
% title(ax_test1_max2,'map curvelet from 2dcurv plus curv');
% set(ax_test2_max2, 'YScale', 'log');
% title(ax_test2_max2,'ridgelet from sum 2dcurv curv');
% set(ax_test3_max2, 'YScale', 'log');
% title(ax_test3_max2,'ridgelet normalized from 2dcurv plus curv');
% set(ax_test4_max2, 'YScale', 'log');
% title(ax_test4_max2,'ridgelet normalized with wavelet from 2dcurv plus curv');
% set(ax_test5_max2, 'YScale', 'log');
% title(ax_test5_max2,'ridgelet normalized with 2dwavelet from 2dcurv plus curv');

% set(ax_test1_curv, 'YScale', 'log');
% title(ax_test1_curv,'kurtosis normalised');
% set(ax_test2_curv, 'YScale', 'log');
% title(ax_test2_curv,'kurtosis curvelet');
% set(ax_test3_curv, 'YScale', 'log');
% title(ax_test3_curv,'4th moment curvelet');

set(ax_test1_curv_mom, 'YScale', 'log');
title(ax_test1_curv_mom,'average normalised fast');
set(ax_test2_curv_mom, 'YScale', 'log');
title(ax_test2_curv_mom,'std curvelet fast');
set(ax_test3_curv_mom, 'YScale', 'log');
title(ax_test3_curv_mom,'skewness curvelet fast');
set(ax_test4_curv_mom, 'YScale', 'log');
title(ax_test4_curv_mom,'kurtosis curvelet fast');
set(ax_test5_curv_mom, 'YScale', 'log');
title(ax_test5_curv_mom,'4th moment curvelet fast');

set(ax3_test1, 'YScale', 'log');
title(ax3_test1,'radon from 2 and 3d filtered');
set(ax3_test2, 'YScale', 'log');
title(ax3_test2,'ridgelet from 2 and 3');
set(ax3_test3, 'YScale', 'log');
title(ax3_test3,'ridgelet normalized from 2 and 3');
set(ax3_test4, 'YScale', 'log');
title(ax3_test4,'ridgelet normalized with wavelet from 2 and 3');

set(ax3_test1s, 'YScale', 'log');
title(ax3_test1s,'radon from 2 and 3d filtered sum');
set(ax3_test2s, 'YScale', 'log');
title(ax3_test2s,'ridgelet from 2 and 3 sum');
set(ax3_test3s, 'YScale', 'log');
title(ax3_test3s,'ridgelet normalized from 2 and 3 sum');
set(ax3_test4s, 'YScale', 'log');
title(ax3_test4s,'ridgelet normalized with wavelet from 2 and 3 sum');

set(ax3_test1_curv_mom, 'YScale', 'log');
title(ax3_test1_curv_mom,'average normalised fast from 2 and 3d filtered');
set(ax3_test2_curv_mom, 'YScale', 'log');
title(ax3_test2_curv_mom,'std curvelet fast from 2 and 3d filtered');
set(ax3_test3_curv_mom, 'YScale', 'log');
title(ax3_test3_curv_mom,'skewness curvelet fast from 2 and 3d filtered');
set(ax3_test4_curv_mom, 'YScale', 'log');
title(ax3_test4_curv_mom,'kurtosis curvelet fast from 2 and 3d filtered');
set(ax3_test5_curv_mom, 'YScale', 'log');
title(ax3_test5_curv_mom,'4th moment curvelet fast from 2 and 3d filtered');

end



% 
% [ anali,sample_id_range_nw,sample_id_range_w ] = curve2d_show_slices_minimal_nor_test_absreal(  )
% 
% 
% nowake=reshape(permute(anali(1,sample_id_range_nw,:,4,1),[1,3,2,4,5]),[1,numel(anali(1,sample_id_range_nw,:,2,1))])
% wake=reshape(permute(anali(2,sample_id_range_w,:,4,1),[1,3,2,4,5]),[1,numel(anali(1,sample_id_range_w,:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices,length(sample_id_range_w)])'
% nowake_slices = reshape(nowake,[slices,length(sample_id_range_nw)])'
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



% %             d3(:)=max(anali3(w_nw,sample,:,4,:),[],3);
% nowake=reshape(permute(anali3(1,sample_id_range_nw,:,4,1),[1,3,2,4,5]),[1,numel(anali3(1,sample_id_range_nw,:,2,1))])
% wake=reshape(permute(anali3(2,sample_id_range_w,:,4,1),[1,3,2,4,5]),[1,numel(anali3(1,sample_id_range_w,:,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% std_wake=std(wake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% wake_slices = reshape(wake,[slices,length(sample_id_range_w)])'
% nowake_slices = reshape(nowake,[slices,length(sample_id_range_nw)])'
% max_wake_slices_=sort(wake_slices')
% max_nowake_slices_=sort(nowake_slices')
% max_wake_slices=max_wake_slices_(end,:)
% max_nowake_slices=max_nowake_slices_(end,:)
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_nowake=std(max_nowake_slices,1)
% std_wake=std(max_wake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
% 


% nowake=reshape(anali3s(1,sample_id_range_nw,4,1),[1,numel(anali3(1,sample_id_range_nw,2,1))])
% wake=reshape(anali3s(2,sample_id_range_w,4,1),[1,numel(anali3(1,sample_id_range_w,2,1))])
% mean_wake=mean(wake)
% mean_nowake=mean(nowake)
% std_nowake=std(nowake,1)
% std_wake=std(wake,1)
% stn_nowake=(nowake-mean_nowake)/std_nowake
% stn_wake=(wake-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)


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
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)
