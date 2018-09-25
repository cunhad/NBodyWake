
% addpath('../../processing');
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));

filename='3.000proj_xz.dat';
nc=2048;
cut=10;
% anal_lev=2;
clims = [-1 cut]

specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
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



fig_test=figure;
% set(gcf, 'Position', [0 0 300 800]);

ax_test=axes(fig_test);

sample_id_range=[1 : length(sample_list_nowake)];

for w_nw=1:2
    
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
        coul='b';
    else
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
        coul='r';
    end
    
    
    for sample = 1:length(sample_id_range)
        
        filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/',filename)
        fid = fopen(filename_nowake);
        scalefactor = fread(fid, [1 1], 'float32','l') ;
        map = fread(fid,[nc nc], 'float32','l') ;
        fclose(fid);
        dc=(map-mean(map(:)))/mean(map(:));
%         dc=map;
%         dc(dc<cut)=cut;
        C = fdct_wrapping(dc,0);
        kurt2(1,1)=1;
        for s=2:length(C)-1
            for w=1:length(C{s})
                %clearvars kurt;s=5; for w=1:length(C{s})
                test(w)=kurtosis(abs(C{s}{w}(:)));
%                 test(w)=kurtosis(abs(C{s}{w}(:)))/E{s}{w};
                
%                 test(w)=mean(abs(C{s}{w}(:)))*kurtosis(abs(C{s}{w}(:)))/std(abs(C{s}{w}(:)));
                %kurt(w)=std(C{s}{w}(:));
            end
%             figure; plot(test)
            mean(test)
            std(test)
            skewness(test)
            kurt2(1,s)=kurtosis(test)
%             kurt2(1,s)=std(test)*kurtosis(test)/mean(test)
            clearvars test;
        end
        kurt2(1)=[];
        test_lv{sample}=   plot(ax_test,[2:length(C)-1],kurt2,coul);
        hold(ax_test,'on');
        
        
    end
    
end


