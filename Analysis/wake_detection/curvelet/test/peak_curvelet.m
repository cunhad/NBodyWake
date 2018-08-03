function [  ] = peak_curvelet( root,root_out,spec,aux_path,aux_path_out,filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%(example)  [proj1d_angles] = box_statistics_dm_data_out_cubic_fast_ap('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/curvelet/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','','10.000xv0.dat');
%(example)  [proj1d_angles] = box_statistics_dm_data_out_cubic_fast_ap('/home/asus/Dropbox/extras/storage/guillimin/', '/home/asus/Dropbox/extras/storage/guillimin/curvelet/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','10.000xv0.dat');

addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));

% [root,root_out,spec,aux_path,aux_path_out,filename,1,1,[0,0,0],[0,0],[1,2]]

cd ../../../2dproj/

[ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_dm_data_out(root,root_out,spec,aux_path,aux_path_out,filename,1,1,[0,0,0],[0,0],[1,2])



n = size(count_sum,1);
sigma = 10;        
is_real = 1;

disp('Compute all thresholds');
% F = ones(n);
F=zeros(n);
% X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
tic, C_zero = fdct_wrapping(F,0); toc;

disp('Compute all thresholds');
F = ones(n);
% F=zeros(n);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
tic, C = fdct_wrapping(X,0); toc;

% Compute norm of curvelets (exact)
E = cell(size(C));
for s=1:length(C)
  E{s} = cell(size(C{s}));
  for w=1:length(C{s})
    A = C{s}{w};
    E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
  end
end

C = fdct_wrapping(count_sum,0);

% Apply thresholding
Ct = C;
for s = 1:length(C)-2
%   thresh = 3*sigma + sigma*(s == length(C));
  for w = 1:length(C{s})
%      Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
     Ct{s}{w} = C_zero{s}{w};
  end
end

for s = length(C)-1:length(C)
%   thresh = 3*sigma + sigma*(s == length(C));
  thresh = 3*sigma;  
  for w = 1:length(C{s})
     Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
%      Ct{s}{w} = C_zero{s}{w};
  end
end

% for s = length(C)
%   thresh = 3*sigma + sigma*(s == length(C));
%   for w = 1:length(C{s})
% %      Ct{s}{w} = C{s}{w}.* (abs(C{s}{w}) > thresh*E{s}{w});
%      Ct{s}{w} = C_zero{s}{w};
%   end
% end

disp(' ');
disp('Take inverse transform of thresholded data: ifdct_wrapping');
tic; restored_img = real(ifdct_wrapping(Ct,1)); toc;

subplot(1,2,1); imagesc(count_sum); colormap gray; axis('image');
subplot(1,2,2); imagesc(restored_img); colormap gray; axis('image');


end