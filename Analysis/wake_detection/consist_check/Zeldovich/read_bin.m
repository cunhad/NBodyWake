function data = read_bin( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fid = fopen(filename);
data=fread(fid, [3 Inf], 'float32','l');
data=transpose(data);
% data=data(3,:);
fclose(fid);

end

