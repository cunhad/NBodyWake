function data = read_xv( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fid = fopen(filename);
fread(fid, [12 1], 'float32','l') ;
data=fread(fid, [3 Inf], 'float32','l');
data=transpose(data);
% data=data(3,:);
fclose(fid);

end

