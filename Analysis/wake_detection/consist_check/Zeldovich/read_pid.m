function data = read_pid( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fid = fopen(filename);
fread(fid, 6, 'int64','l') ;
data=fread(fid, [1 Inf], 'int64','l');
data=transpose(data);
% data=data(3,:);
fclose(fid);

end

