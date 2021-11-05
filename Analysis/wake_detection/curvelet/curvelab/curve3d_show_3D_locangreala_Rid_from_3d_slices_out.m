function [ anali,anali_curv_mom,anali3,anali3s,anali3_curv_mom,sample_id_range_nw,sample_id_range_w ] = curve3d_show_3D_locangreala_Rid_from_3d_slices_out(  )



clearvars;

%example:
% nowake(:)=anali(1,:,3,4);
%wake(:)=anali(2,:,3,4);

% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct3d'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/asus/Dropbox/Disrael/Doutorado/Research/NBodyWake/production/Analysis/wake_detection/ridgelet/ppft3_nomex/ppft3'))
% run('/home/asus/Dropbox/Disrael/Doutorado/Research/NBodyWake/production/Analysis/wake_detection/ridgelet/ppft3_nomex/ppft3/initpath.m')

% addpath(genpath('../../ridgelet/Ridgelet3d_interp_forwards_dev3.m'));
addpath(genpath('../../ridgelet/'));


filename='_2dproj_3dcurv_z3_data_sl';
aux_path='/half_lin_cutoff_half_tot_pert_nvpw_v0p6'
nc=512;
rf=0.5;  %rf=nc/1024
lev_2d=1;
lev_2drig=1;
lev_3d=1;
lev_3drig=1;
slices=32;
size_mpc=4;
partition2d=1;
partition3rd=1;
step_of_degree=1*(180/256);
wavel_removal_factor=1/2;

sample_id_range_nw=[1:2];
sample_id_range_w=[1:2];
% sample_id_range_nw=[4,7];
% sample_id_range_w=[3,7];
% sample_id_range_nw=[4,7];
% sample_id_range_w=[4,7];
% 
% % sample_id_range=[3, 7];
% % sample_id_range=[1 : length(sample_list_nowake)];

% sample_id_range_nw=[1];
% sample_id_range_w=[8];

% sample_id_range_nw=[1:10];
% sample_id_range_w=[1:10];

% sample_id_range_nw=[1];
% sample_id_range_w=[1];


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


% specs_path_list_nowake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_data_test','/4Mpc_2048c_1024p_zi63_nowakem')
specs_path_list_nowake=strcat('/home/disraelcunha/Documents/graham/cubep3m/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_data_test','/4Mpc_2048c_1024p_zi63_nowakem')
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

% specs_path_list_wake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_data_test','/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m')
specs_path_list_wake=strcat('/home/disraelcunha/Documents/graham/cubep3m/ht/data_cps',num2str(slices),'_',num2str(nc),'_2dclara-l1lr1na1024_to_3dparclar-l1lr1_data_test','/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m')
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,aux_path);
% sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p55');
% sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_wrong');
% sample_list_wake=sort_nat(sample_list_wake)


% %
% % F=zeros(nc);
% % C_zero = fdct_wrapping(F,0);
% F = ones(nc);
% X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
% %X = F * sqrt(prod(size(F)));
% %C = fdct_wrapping(X,0);
% C = fdct_wrapping(X,0);
% %C = fdct_wrapping(F,0);
% E = cell(size(C));
% for s=1:length(C)
%     E{s} = cell(size(C{s}));
%     for w=1:length(C{s})
%         A = C{s}{w};
%         E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
%     end
% end
% 
% 
% F=zeros(nc);
% C_zero = fdct_wrapping(F,0);

% nc_anal3d=nc/partition2d;
% 
% % F2 = ones(nc_anal3d,nc_anal3d,slices/partition3rd);
% % % F2 = ones(new_nc,new_nc,slices);
% % X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
% X2 = zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
% X2(1+nc_anal3d/2,1+nc_anal3d/2,1+slices/(2*partition3rd))=(nc_anal3d*nc_anal3d*slices/partition3rd)^(1/3);
% C2 = fdct3d_forward(X2);
% E2 = cell(size(C2));
% for s=1:length(C2)
%     E2{s} = cell(size(C2{s}));
%     for w=1:length(C2{s})
%         A2 = C2{s}{w};
%         E2{s}{w} = sqrt(sum(sum(sum(A2.*conj(A2)))) / prod(size(A2)));
%     end
% end
% 
% 
% 
% F2=zeros(nc_anal3d,nc_anal3d,slices/partition3rd);
% % F2=zeros(new_nc,new_nc,slices);
% C_zero2 = fdct3d_forward(F2);
% 
% 
% %radon
% 
% % nb_m=max(nb,slices)
% % nb_m=min(nb,slices)
% nb_m=sqrt(nb*slices);
% Z2 = ones(nb_m,nb_m,nb_m);
% Rad_norm = radon3(Z2);




% % fig_curv=figure;
% % fig_test1=figure;
% fig_test2=figure;
% fig_test3=figure;
% fig_test4=figure;
% 
% % fig_test1_m=figure;
% % fig_test2_m=figure;
% % fig_test3_m=figure;
% % fig_test4_m=figure;
% 
% % fig_test1_max=figure;
% % fig_test2_max=figure;
% % fig_test3_max=figure;
% % fig_test4_max=figure;
% 
% % fig_test1_max2=figure;
% % fig_test2_max2=figure;
% % fig_test3_max2=figure;
% % fig_test4_max2=figure;
% % fig_test5_max2=figure;
% 
% % fig_test1_curv=figure;
% % fig_test2_curv=figure;
% % fig_test3_curv=figure;
% 
% 
% fig_test1_curv_mom=figure;
% fig_test2_curv_mom=figure;
% fig_test3_curv_mom=figure;
% fig_test4_curv_mom=figure;
% fig_test5_curv_mom=figure;
% 
% fig3_test1=figure;
% fig3_test2=figure;
% fig3_test3=figure;
% fig3_test4=figure;
% 
% fig3_test1_curv_mom=figure;
% fig3_test2_curv_mom=figure;
% fig3_test3_curv_mom=figure;
% fig3_test4_curv_mom=figure;
% fig3_test5_curv_mom=figure;
% 
% % fig3_test1s=figure;
% % fig3_test2s=figure;
% % fig3_test3s=figure;
fig3_test4s=figure;


fig_hist = figure;
% 
% % ax_curv=axes(fig_curv);
% % ax_test1=axes(fig_test1);
% ax_test2=axes(fig_test2);
% ax_test3=axes(fig_test3);
% ax_test4=axes(fig_test4);
% 
% % ax_test1_m=axes(fig_test1_m);
% % ax_test2_m=axes(fig_test2_m);
% % ax_test3_m=axes(fig_test3_m);
% % ax_test4_m=axes(fig_test4_m);
% % 
% % ax_test1_max=axes(fig_test1_max);
% % ax_test2_max=axes(fig_test2_max);
% % ax_test3_max=axes(fig_test3_max);
% % ax_test4_max=axes(fig_test4_max);
% 
% % ax_test1_max2=axes(fig_test1_max2);
% % ax_test2_max2=axes(fig_test2_max2);
% % ax_test3_max2=axes(fig_test3_max2);
% % ax_test4_max2=axes(fig_test4_max2);
% % ax_test5_max2=axes(fig_test5_max2);
% 
% % ax_test1_curv=axes(fig_test1_curv);
% % ax_test2_curv=axes(fig_test2_curv);
% % ax_test3_curv=axes(fig_test3_curv);
% 
% 
% ax_test1_curv_mom=axes(fig_test1_curv_mom);
% ax_test2_curv_mom=axes(fig_test2_curv_mom);
% ax_test3_curv_mom=axes(fig_test3_curv_mom);
% ax_test4_curv_mom=axes(fig_test4_curv_mom);
% ax_test5_curv_mom=axes(fig_test5_curv_mom);
% 
% ax3_test1=axes(fig3_test1);
% ax3_test2=axes(fig3_test2);
% ax3_test3=axes(fig3_test3);
% ax3_test4=axes(fig3_test4);
% 
% % ax3_test1s=axes(fig3_test1s);
% % ax3_test2s=axes(fig3_test2s);
% % ax3_test3s=axes(fig3_test3s);
ax3_test4s=axes(fig3_test4s);
ax_hist=axes(fig_hist);
% 
% ax3_test1_curv_mom=axes(fig3_test1_curv_mom);
% ax3_test2_curv_mom=axes(fig3_test2_curv_mom);
% ax3_test3_curv_mom=axes(fig3_test3_curv_mom);
% ax3_test4_curv_mom=axes(fig3_test4_curv_mom);
% ax3_test5_curv_mom=axes(fig3_test5_curv_mom);



for w_nw=1:2
% for w_nw=2
% for w_nw=1
    
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

        
%         map_3d_slices_pre=zeros(nc,nc,slices);
        map_3d_slices=zeros(nc,nc,slices);
        map_3d_slices_filt2d=zeros(nc,nc,slices);
        map_3d_slices_filt2a3d=zeros(nc,nc,slices);
        
        for slice_id=1:slices
            
            %             sample_id=(slices*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_',num2str(rf),'rf_0-0-0pv_1.5708-0-0ra/3d/dm/',ch,filename,num2str(slice_id),'.bin')
            fid = fopen(filename_nowake);
            map = fread(fid,[nc nc], 'float32','l') ;
            fclose(fid);
            
            %             map(map<=1)=1;%to remove problem with holes
%                         map_3d_slices(:,:,slice_id)=log(map);

%                 map_3d_slices(:,:,slice_id)=log(map+1);
            map_3d_slices(:,:,slice_id)=map;
            
            
        end
        


[ Radon_Z_ones,interp_Z_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(ones(size(map_3d_slices)));
[ Radon_Z_,interp_Z_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(map_3d_slices);
Radon_Z_norm_ = real(Radon_Z_)./real(Radon_Z_ones);
dataZ = Radon_Z_norm_(1,:);

[ Radon_Y_ones,interp_Y_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(ones(size(permute(map_3d_slices,[1 3 2]))));
[ Radon_Y_,interp_Y_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(permute(map_3d_slices,[1 3 2]));
Radon_Y_norm_ = real(Radon_Y_)./real(Radon_Y_ones);
dataY = Radon_Y_norm_(1,:);

[ Radon_X_ones,interp_X_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(ones(size(permute(map_3d_slices,[3 2 1]))));
[ Radon_X_,interp_X_info,ncx_,ncy_,ncz_] = Ridgelet3d_interp_forwards_dev3(permute(map_3d_slices,[3 2 1]));
Radon_X_norm_ = real(Radon_X_)./real(Radon_X_ones);
dataX = Radon_X_norm_(1,:);

data = [dataX,dataY,dataZ];
% data = [dataX,dataY];
         
%             anali3s(w_nw,sample,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(this(:))];
%             anali3s(w_nw,sample,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(R(:))];
%             anali3s(w_nw,sample,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(R_nor(:))];
            anali3s(w_nw,sample,4,:)=[max(data(:)),std(data(:)),max(data(:))/std(data(:)),kurtosis(data(:))];

%             
            
        
        
%             a3s(:)=(anali3s(w_nw,sample,1,:));
%             b3s(:)=(anali3s(w_nw,sample,2,:));
%             c3s(:)=(anali3s(w_nw,sample,3,:));
            d3s(:)=(anali3s(w_nw,sample,4,:));
            
%             test1max3_lvs{sample_id}=   plot(ax3_test1s,a3s,coul);
%             test2max3_lvs{sample_id}=   plot(ax3_test2s,b3s,coul);
%             test3max3_lvs{sample_id}=   plot(ax3_test3s,c3s,coul);
            test4max3_lvs{sample}=   plot(ax3_test4s,d3s,coul);
            %             test5max2_lv{sample_id}=   plot(ax_test5_max2,e_max2,coul);
            
           
            
         clearvars a3s b3s c3s d3s
         
%         hold(ax3_test1s,'on');
%         hold(ax3_test2s,'on');
%         hold(ax3_test3s,'on');
        hold(ax3_test4s,'on');
        
        
        if (w_nw==1&&sample==1)
            figure;
            h1 = histogram(data);
            set(gca, 'YScale', 'log')
            Hist{w_nw,sample} = h1.BinCounts;
            BinWidth1 = h1.BinWidth;
            BinEdgesHist{w_nw,sample} = h1.BinEdges;
            
%             plot(ax_hist,BinEdgesHist{1,1}(1:end-1),Hist{1,1},coul);
%             hold(ax_hist,'on');

            
            clearvars h1
        else
            
            figure;
            h = histogram(data,'BinWidth',BinWidth1);    
            set(gca, 'YScale', 'log')
            Hist{w_nw,sample} = h.BinCounts;
            BinEdgesHist{w_nw,sample} = h.BinEdges;      
            
%             plot(ax_hist,BinEdgesHist{w_nw,sample}(1:end-1),Hist{w_nw,sample},coul);
%             hold(ax_hist,'on');

            
            clearvars h
            
        end

% Hist{w_nw,sample} = h.BinCounts;
% 
% {w_nw,sample} 

% figure; 
% h{w_nw,sample} = histogram(data);
% % h2 = histogram(data,'BinWidth',h{1,1}.BinWidth);
% set(gca, 'YScale', 'log')



    end
    
    
    
end




dom_min=min(BinEdgesHist{1});
for i = 2:prod(size(BinEdgesHist))
    dom_min_=min(BinEdgesHist{i})
    if (dom_min_<dom_min) dom_min = dom_min_;end
    
end

dom_max=max(BinEdgesHist{1});
for i = 2:prod(size(BinEdgesHist))
    dom_max_=max(BinEdgesHist{i})
    if (dom_max_<dom_max) dom_max = dom_max_;end
    
end

domain = dom_min_:BinWidth1:dom_max;

[a b] = size(BinEdgesHist)

for i = 1:b    
    Hist_comDom_nw(i,:) = interp1(BinEdgesHist{1,i}(1:end-1)',Hist{1,i}',domain','nearest','extrap');
end

for i = 1:b    
    Hist_comDom_w(i,:) = interp1(BinEdgesHist{2,i}(1:end-1)',Hist{2,i}',domain','nearest','extrap');
end



figure;
plot(domain,sum(Hist_comDom_nw))
hold on
plot(domain,sum(Hist_comDom_w))
set(gca, 'YScale', 'log')


% figure;
% hist(BinEdgesHist{1,1}(1:end-1), Hist{1,1})
% set(gca, 'YScale', 'log')
% 
% 
% f1 = 

% figure;
% plot(BinEdgesHist{1,1}(1:end-1),Hist{1,1})
% set(ax_hist, 'YScale', 'log')
% 
% hold(ax_hist,'on');
% plot(ax_hist,BinEdgesHist{1,1}(1:end-1),Hist{1,1},coul);
% plot(ax_hist,BinEdgesHist{2,1}(1:end-1),Hist{2,1},coul);


% 
% % set(ax_curv, 'YScale', 'log');
% % title(ax_curv,'curvelet kurtosis')
% % set(ax_test1, 'YScale', 'log');
% % title(ax_test1,'all map curvelet');
% set(ax_test2, 'YScale', 'log');
% title(ax_test2,'ridgelet');
% set(ax_test3, 'YScale', 'log');
% title(ax_test3,'ridgelet normalized');
% set(ax_test4, 'YScale', 'log');
% title(ax_test4,'ridgelet normalized with wavelet');
% 
% % 
% % set(ax_test1_m, 'YScale', 'log');
% % title(ax_test1_m,'map curvelet sum max');
% % set(ax_test2_m, 'YScale', 'log');
% % title(ax_test2_m,'ridgelet from max');
% % set(ax_test3_m, 'YScale', 'log');
% % title(ax_test3_m,'ridgelet normalized from max');
% % set(ax_test4_m, 'YScale', 'log');
% % title(ax_test4_m,'ridgelet normalized with wavelet from max');
% % 
% % 
% % set(ax_test1_max, 'YScale', 'log');
% % title(ax_test1_max,'map curvelet from curv');
% % set(ax_test2_max, 'YScale', 'log');
% % title(ax_test2_max,'ridgelet from curv');
% % set(ax_test3_max, 'YScale', 'log');
% % title(ax_test3_max,'ridgelet normalized from curv');
% % set(ax_test4_max, 'YScale', 'log');
% % title(ax_test4_max,'ridgelet normalized with wavelet from curv');
% % 
% % set(ax_test1_max2, 'YScale', 'log');
% % title(ax_test1_max2,'map curvelet from 2dcurv plus curv');
% % set(ax_test2_max2, 'YScale', 'log');
% % title(ax_test2_max2,'ridgelet from sum 2dcurv curv');
% % set(ax_test3_max2, 'YScale', 'log');
% % title(ax_test3_max2,'ridgelet normalized from 2dcurv plus curv');
% % set(ax_test4_max2, 'YScale', 'log');
% % title(ax_test4_max2,'ridgelet normalized with wavelet from 2dcurv plus curv');
% % set(ax_test5_max2, 'YScale', 'log');
% % title(ax_test5_max2,'ridgelet normalized with 2dwavelet from 2dcurv plus curv');
% 
% % set(ax_test1_curv, 'YScale', 'log');
% % title(ax_test1_curv,'kurtosis normalised');
% % set(ax_test2_curv, 'YScale', 'log');
% % title(ax_test2_curv,'kurtosis curvelet');
% % set(ax_test3_curv, 'YScale', 'log');
% % title(ax_test3_curv,'4th moment curvelet');
% 
% set(ax_test1_curv_mom, 'YScale', 'log');
% title(ax_test1_curv_mom,'average normalised fast');
% set(ax_test2_curv_mom, 'YScale', 'log');
% title(ax_test2_curv_mom,'std curvelet fast');
% set(ax_test3_curv_mom, 'YScale', 'log');
% title(ax_test3_curv_mom,'skewness curvelet fast');
% set(ax_test4_curv_mom, 'YScale', 'log');
% title(ax_test4_curv_mom,'kurtosis curvelet fast');
% set(ax_test5_curv_mom, 'YScale', 'log');
% title(ax_test5_curv_mom,'4th moment curvelet fast');
% 
% set(ax3_test1, 'YScale', 'log');
% title(ax3_test1,'radon from 2 and 3d filtered');
% set(ax3_test2, 'YScale', 'log');
% title(ax3_test2,'ridgelet from 2 and 3');
% set(ax3_test3, 'YScale', 'log');
% title(ax3_test3,'ridgelet normalized from 2 and 3');
% set(ax3_test4, 'YScale', 'log');
% title(ax3_test4,'ridgelet normalized with wavelet from 2 and 3');
% 
% set(ax3_test1s, 'YScale', 'log');
% title(ax3_test1s,'radon from 2 and 3d filtered sum');
% set(ax3_test2s, 'YScale', 'log');
% title(ax3_test2s,'ridgelet from 2 and 3 sum');
% set(ax3_test3s, 'YScale', 'log');
% title(ax3_test3s,'ridgelet normalized from 2 and 3 sum');
set(ax3_test4s, 'YScale', 'log');
title(ax3_test4s,'ridgelet normalized with wavelet from 2 and 3 sum');

% set(ax3_test1_curv_mom, 'YScale', 'log');
% title(ax3_test1_curv_mom,'average normalised fast from 2 and 3d filtered');
% set(ax3_test2_curv_mom, 'YScale', 'log');
% title(ax3_test2_curv_mom,'std curvelet fast from 2 and 3d filtered');
% set(ax3_test3_curv_mom, 'YScale', 'log');
% title(ax3_test3_curv_mom,'skewness curvelet fast from 2 and 3d filtered');
% set(ax3_test4_curv_mom, 'YScale', 'log');
% title(ax3_test4_curv_mom,'kurtosis curvelet fast from 2 and 3d filtered');
% set(ax3_test5_curv_mom, 'YScale', 'log');
% title(ax3_test5_curv_mom,'4th moment curvelet fast from 2 and 3d filtered');

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


nowake=reshape(anali3s(1,sample_id_range_nw,4,1),[1,numel(anali3s(1,sample_id_range_nw,2,1))])
wake=reshape(anali3s(2,sample_id_range_w,4,1),[1,numel(anali3s(1,sample_id_range_w,2,1))])
mean_wake=mean(wake)
mean_nowake=mean(nowake)
std_nowake=std(nowake,1)
std_wake=std(wake,1)
stn_nowake=(nowake-mean_nowake)/std_nowake
stn_wake=(wake-mean_nowake)/std_nowake
mean_stn=mean(stn_wake)
mean_stn-std_stn
significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)


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
