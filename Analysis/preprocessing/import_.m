function [ dat ] = import_( path_and_name,formatSpec,n_col)

%for matlab R2016b
%import data

fid=fopen(path_and_name);
%formatSpec = '%f %f';
size = [n_col Inf];
dat=fscanf(fid,formatSpec,size);
dat=transpose(dat);

end
