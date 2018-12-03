function [ map_3d_slices ] = curve3d_slices(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% (obs) on asus only works with 1024*1024*16

addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct3d'));

addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));

filename='_2dproj_z3_data_sl';
nc=1024;
lev=2;
sigma = 5;
slice=32;

F = ones(nc,nc,slice/2);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
C = fdct3d_forward(X);
E = cell(size(C));
for s=1:length(C)
    E{s} = cell(size(C{s}));
    for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
    end
end

specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};

specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw');


sample_id_range=[1 : length(sample_list_nowake)];

for w_nw=2
% for w_nw=1:2
    
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
    
    
    
%     for sample = 1:length(sample_id_range)
      for sample = 1:1
        
        
        map_3d_slices=zeros(nc,nc,slice);
        
        for slice_id=1:slice
            
            sample_id=(slice*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/',ch,filename,num2str(slice_id),'.bin')
            fid = fopen(filename_nowake);
            %         scalefactor = fread(fid, [1 1], 'float32','l') ;
            map = fread(fid,[nc nc], 'float32','l') ;
            fclose(fid);
            
            map(map<=1)=1;%to remove problem with holes
            map_3d_slices(:,:,slice_id)=log(map);
            
            
        end
        
        C = fdct3d_forward(map_3d_slices(:,:,1:slices/2));
        
        F=zeros(nc,nc,slice/2);
        C_zero = fdct3d_forward(F);
        Ct = C;
        for s = 1:length(C)
            for w = 1:length(C{s})
                Ct{s}{w} = C_zero{s}{w};
            end
        end
        
        aux_count=1;
        for s = 1:length(C)
            %                 thresh=0;
            thresh = sigma + sigma*(s == length(C));
            for w = 1:length(C{s})
%                                      Ct{s}{w} = C{s}{w};
                                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
%                 Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                %                 Ct{s}{w} = C{s}{w};
                %                     curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
            end
            %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
            aux_count=aux_count+1;
        end
        
        
        
        map_3d_slices_filt = real(fdct3d_inverse(Ct));
        
        
      end
    
    
    
end

% 
% figure; imagesc([2/1024:4/1024:4],[2/1024:4/1024:4],map_3d_slices(:,:,12)); colorbar; axis('image');
% figure; imagesc([2/1024:4/1024:4],[2/1024:4/1024:4],map_3d_slices_filt(:,:,12)); colorbar; axis('image');
