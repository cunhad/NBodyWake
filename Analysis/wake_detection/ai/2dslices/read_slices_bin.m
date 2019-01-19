function slice_2d = read_slices_bin( filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



fid = fopen(filename);
slice_2d = fread(fid,[1024 1024], 'float32','l') ;
% data=transpose(data);
fclose(fid);


end

