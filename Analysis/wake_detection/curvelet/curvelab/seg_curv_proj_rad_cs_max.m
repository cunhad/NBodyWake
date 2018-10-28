function [ anali ] = seg_curv_proj_rad_cs_max(  )


%example:
% nowake(:)=anali(1,:,3,4);
%wake(:)=anali(2,:,3,4);

% addpath('../../processing');
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));

filename='_2dproj_z3_data_sl';
nc=1024;
trsh=20;
cut=1;
lev=3;
slice=8;
anal_lev=1;

specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps8_1012_nopro/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps8_1012_nopro/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw');
% sample_list_wake=sort_nat(sample_list_wake)


%
% F=zeros(nc);
% C_zero = fdct_wrapping(F,0);
F = ones(nc);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%X = F * sqrt(prod(size(F)));
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

fig_curv=figure;
fig_test1=figure;
fig_test2=figure;
fig_test3=figure;
fig_test4=figure;

fig_test1_m=figure;
fig_test2_m=figure;
fig_test3_m=figure;
fig_test4_m=figure;


ax_curv=axes(fig_curv);
ax_test1=axes(fig_test1);
ax_test2=axes(fig_test2);
ax_test3=axes(fig_test3);
ax_test4=axes(fig_test4);

ax_test1_m=axes(fig_test1_m);
ax_test2_m=axes(fig_test2_m);
ax_test3_m=axes(fig_test3_m);
ax_test4_m=axes(fig_test4_m);

sample_id_range=[1 : length(sample_list_nowake)];

for w_nw=1:2
    
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
        ch='_7';
        coul='b';
    else
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
        ch='_4';
        coul='r';
    end
    
    
    
    for sample = 1:length(sample_id_range)
        
        for slice_id=1:slice
            
            sample_id=(slice*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/',ch,filename,num2str(slice_id),'.bin')
            fid = fopen(filename_nowake);
            %         scalefactor = fread(fid, [1 1], 'float32','l') ;
            map = fread(fid,[nc nc], 'float32','l') ;
            fclose(fid);
            dc=(map-mean(map(:)))/mean(map(:));
            %         dc=map;
            %         dc(dc>cut)=cut;
            
            %         nc_red=nc/red;
            %         conv_=ones(red);
            %         dc_red = conv2(dc,conv_,'valid');
            %         dc = dc_red(1:red:end,1:red:end)/(red*red);
            %
            
            
            
            dc_cut=dc;
            dc_cut(dc_cut>cut)=cut;
            %         dc_cut(dc_cut>cut)=-1;
            %         dc_cut = edge(dc_cut,'canny');
            %         dc=double(dc_cut);
            
            thresh = multithresh(dc_cut,trsh);
            seg_I = imquantize(dc_cut,thresh);
            %         test=dc_cut;
            %         test=seg_I;
            test=log(seg_I);
            
            
            
            n = length(test);
            C = fdct_wrapping(test,0);
            F=zeros(n);
            C_zero = fdct_wrapping(F,0);
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero{s}{w};
                end
            end
            aux_count=1;
            for s = length(C)-lev:length(C)
                thresh=0;
                for w = 1:length(C{s})
                    Ct{s}{w} = C{s}{w};
                    curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                end
                curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                aux_count=aux_count+1;
            end
            
            
            
            BW2 = real(ifdct_wrapping(Ct,0));
            
            
            
            anali(w_nw,sample,slice_id,1,:)=[max(BW2(:)),std(BW2(:)),max(BW2(:))/std(BW2(:)),kurtosis(kurtosis(BW2)),kurtosis(BW2(:))];
            
            % theta = 0:180/nc:180;
            theta = 0:180;
            [R,xp] = radon(BW2,theta);
            
            unit=ones(nc);
            [R_u,xp] = radon(unit,theta);
            
            frac_cut=0.5;
            R_nor=R;
            R_nor(R_u>nc*frac_cut)=R_nor(R_u>nc*frac_cut)./R_u(R_u>nc*frac_cut);
            R_nor(R_u<=nc*frac_cut)=0;
            
            n_levels=floor(log2(length(R_nor(:,1))));
            R_nor_filt=zeros(size(R_nor));
            for i=1:length(R(1,:))
                [dc_dwt,levels] = wavedec(R(:,i),n_levels,'db1');
                D = wrcoef('d',dc_dwt,levels,'db1',lev);
                D(2465-200:end)=0;
                D(1:448+200)=0;
                R_nor_filt(:,i)=D;
            end
            R_nor_filt(2465-200:end,:)=[];
            R_nor_filt(1:448+200,:)=[];
            
            
            
            anali(w_nw,sample,slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
            anali(w_nw,sample,slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
            anali(w_nw,sample,slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
            
            
            aa(:)=curv2(w_nw,sample,slice_id,:);
            
            a(:)=anali(w_nw,sample,slice_id,1,:);
            b(:)=anali(w_nw,sample,slice_id,2,:);
            c(:)=anali(w_nw,sample,slice_id,3,:);
            d(:)=anali(w_nw,sample,slice_id,4,:);
            
            aa
            a
            b
            c
            d
            
            
            test_curv_lv{sample_id}=   plot(ax_curv,aa,coul);
            test1_lv{sample_id}=   plot(ax_test1,a,coul);
            test2_lv{sample_id}=   plot(ax_test2,b,coul);
            test3_lv{sample_id}=   plot(ax_test3,c,coul);
            test4_lv{sample_id}=   plot(ax_test4,d,coul);
            
            clearvars a b c d aa
            
            %         clearvars R R_nor R_nor_filt kurt2;
            
            hold(ax_curv,'on');
            hold(ax_test1,'on');
            hold(ax_test2,'on');
            hold(ax_test3,'on');
            hold(ax_test4,'on');
            
        end
        
        a(:)=max(anali(w_nw,sample,:,1,:),[],3);
        b(:)=max(anali(w_nw,sample,:,2,:),[],3);
        c(:)=max(anali(w_nw,sample,:,3,:),[],3);
        d(:)=max(anali(w_nw,sample,:,4,:),[],3);
        
        a
        b
        c
        d
        
        
        test1m_lv{sample_id}=   plot(ax_test1_m,a,coul);
        test2m_lv{sample_id}=   plot(ax_test2_m,b,coul);
        test3m_lv{sample_id}=   plot(ax_test3_m,c,coul);
        test4m_lv{sample_id}=   plot(ax_test4_m,d,coul);
        
        clearvars a b c d aa
        
        %         clearvars R R_nor R_nor_filt kurt2;
        
        hold(ax_test1_m,'on');
        hold(ax_test2_m,'on');
        hold(ax_test3_m,'on');
        hold(ax_test4_m,'on');
        
    end
    
    
    
end

set(ax_curv, 'YScale', 'log');
set(ax_test1, 'YScale', 'log');
set(ax_test2, 'YScale', 'log');
set(ax_test3, 'YScale', 'log');
set(ax_test4, 'YScale', 'log');

set(ax_test1_m, 'YScale', 'log');
set(ax_test2_m, 'YScale', 'log');
set(ax_test3_m, 'YScale', 'log');
set(ax_test4_m, 'YScale', 'log');

end

