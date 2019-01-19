function slice_2d = log_read_slices_bin( filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



fid = fopen(filename);
slice_2d = fread(fid,[1024 1024], 'float32','l') ;
% data=transpose(data);
slice_2d(slice_2d<=1)=1;%to remove problem with holes
            slice_2d=log(slice_2d);
            maxim=max(max(slice_2d));
            slice_2d=slice_2d/maxim;
fclose(fid);


end

