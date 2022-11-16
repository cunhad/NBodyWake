function slice_2d = read_slices_bin_slices( filename)
% function slice_2d = read_slices_binAll( filename)
% function [slice_2d,userdata,done] = read_slices_binAll(filename,userdata)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% a = userdata
% display(userdata)

% fid = fopen(filename);
% slice_2d = fread(fid,[512 512], 'float32','l') ;
% % data=transpose(data);
% fclose(fid);


% 
% 
% 
% slice_aux = split(filename,'.');
% slice_aux2 = split(slice_aux{end-1},'All');
% slice = str2num(slice_aux2{end});
% 
% 


% 
% 
% %read 3d
% 
% fid = fopen(filename);
% slice_3d = fread(fid,512*512*32, 'float32','l') ;  %Optimize this, so only the needed data is loaded
% % data=transpose(data);
% fclose(fid);
% slice_3d = reshape(slice_3d,512,512,32);
% slice_2d = sum(slice_3d,3);





% %read slice from 323d
% 
% fid = fopen(filename);
% slice_3d = fread(fid,512*512*32, 'float32','l') ;  %Optimize this, so only the needed data is loaded
% % data=transpose(data);
% fclose(fid);
% slice_3d = reshape(slice_3d,512,512,32);
% 
% slice_2d = {slice_3d(:,:,1);slice_3d(:,:,2)}





% reads part of the file

nc = 512;



% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_sl32All3.bin'

slice_aux = split(filename,'/');
slice_aux2 = split(slice_aux{end},'All');
slice_aux3 = split(slice_aux2{1},'sl');
slices = str2num(slice_aux3{end});
slice_aux4 = split(slice_aux2{end},'.');
sliceID = str2num(slice_aux4{1});

filename_rec = join({slice_aux3{1};'slAll.bin'},'');

filename_in = join({slice_aux{1:end-1},filename_rec{1}},'/');
filename_in = filename_in{:};




skip = 4*(nc*nc)*(sliceID-1);     %each number has 4 `bytes to skip

fid = fopen(filename_in);
% slice_2d = fread(fid,nc*nc, 'float32','l') ; 
fseek(fid,skip,'bof');
slice_2d = fread(fid,[512 512], 'float32','l') ; 
% slice_2d = reshape(slice_2d,512,512);
fclose(fid);


end

