function [ xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info( path_files,spec_files,aux_path_files)
%   This function takes the path to the xv output of CUBEP3M and saves a list
%   with the filenames

%(example)[ xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info('/home/asus/Dropbox/extras/storage/graham/small_res/','64Mpc_96c_48p_zi255_nowakem','/sample1001/');

path_files_read=strcat(path_files,'data_list_out/',spec_files,aux_path_files);


cd('../processing');

%read xv list

fid = fopen(strcat(path_files_read,'xv_files_list.txt'),'r');
xv_files_list=textscan(fid,'%s');
xv_files_list=xv_files_list{1};
xv_files_list=transpose(xv_files_list);
fclose(fid);


%read redshift list


fid = fopen(strcat(path_files_read,'z_list.txt'),'r');
redshift_list=textscan(fid,'%s');
redshift_list=redshift_list{1};
redshift_list=transpose(redshift_list);
fclose(fid);

%read node list


fid = fopen(strcat(path_files_read,'nodes_list.txt'),'r');
nodes_list=textscan(fid,'%s');
nodes_list=nodes_list{1};
nodes_list=transpose(nodes_list);
fclose(fid);

 cd('../preprocessing');

 %reads the specifications and extract the information on variables
spec_arr = strsplit(spec_files,'_');

%extract the box size

size_box = spec_arr(1);
size_box = char(size_box);
size_box = size_box(1:end-3);
size_box = str2num(size_box);

%extract the number of cells per dimension

nc = spec_arr(2);
nc = char(nc);
nc = nc(1:end-1);
nc = str2num(nc);

%extract the number of particle per dimension

np = spec_arr(3);
np = char(np);
np = np(1:end-1);
np = str2num(np);

 
%extract the initial redshift of the simulation
  
 zi = spec_arr(4);
 zi = char(zi);
 zi = zi(3:end);
 zi = str2num(zi);
 
 
 % extracts the informations of the wake if there is one
 
 wake_spec = spec_arr(5);
 wake_spec = char(wake_spec);
 if wake_spec(1)=='n'
    wake_or_no_wake='no wake';
    multiplicity_of_files=wake_spec(end); 
    Gmu=0;
    ziw=0;
 end
 if wake_spec(1)=='w'
     wake_or_no_wake='wake';
     wake_spec2=strsplit(wake_spec,{'u','t10m','zi'},'CollapseDelimiters',true);
     Gmu=str2num(char(wake_spec2(2)))*10^(-str2num(char(wake_spec2(3))));
     ziw=char(wake_spec2(4));
     ziw=str2num(ziw(1:end-1));
     multiplicity_of_files=char(wake_spec(end));
 end
 




end