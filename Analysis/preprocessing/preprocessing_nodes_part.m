function [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes_part( path,spec,aux_path,filename,part,part_id)
%   This function takes the phase space output from CUBEP3M and returns
%   some relevant information plus the global positions of all particles at the
%   corresponding redshift and node volume


%(example) [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes_part( '/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','0.000xv0.dat',8,1);



%   Detailed explanation:

%reads the specifications and extract the information on variables
spec_arr = strsplit(spec,'_');

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
 
  %extracts the redshift of the file to be analised
 
 z=char(filename);
 z=str2num(z(1:end-7));
 
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
 
 %extract the information of the multiplicity of the files to be analysed
 
 if multiplicity_of_files=='s'
     path_file_in=strcat(path,spec,aux_path,filename);
 end
 if multiplicity_of_files=='m'
     path_file_in=strcat(path,spec,aux_path,filename);
 end
 
fid = fopen(path_file_in);
directory = dir(path_file_in);
particles=(directory.bytes-48)/24;
headear = fread(fid, [12 1], 'float32','l') ;
for i=1:part_id
data=fread(fid, [6 particles/part], 'float32','l');
end
fclose(fid);
Pos = data(1:3,:);

%Pos=transpose(Pos);
%Pos=mod(Pos,nc);
 
%in this part we will get the position of the wake taking into acount the
%node structure

 node=char(filename);
 node=str2num(node(strfind(filename, 'xv')+2:strfind(filename,'.dat')-1));
 [ nodes_list redshift_list ] = preprocessing_many_nodes(path,spec,aux_path);
 
 number_node_dim=nthroot(numel(nodes_list), 3);
 k_node=floor(node/number_node_dim^2);
 res=mod(node,number_node_dim^2);
 j_node=floor(res/number_node_dim);
 i_node=mod(res,number_node_dim);
 
 
 Pos(1,:)=Pos(1,:)+(nc/number_node_dim)*i_node;
 Pos(2,:)=Pos(2,:)+(nc/number_node_dim)*j_node;
 Pos(3,:)=Pos(3,:)+(nc/number_node_dim)*k_node;
 
 
%  XM=nc;
%  Xm=0;
%  YM=nc;
%  Ym=0;
%  ZM=nc;
%  Zm=0;
%  ncx=nc;
%  ncy=nc;
%  ncz=nc;
%  halfx=(XM+Xm)/2;
%  halfy=(YM+Ym)/2;
%  halfz=(ZM+Zm)/2;
%  limxinf=halfx-percentage_analysed*ncx/2;
%  limxsup=halfx+percentage_analysed*ncx/2;
%  limyinf=halfy-percentage_analysed*ncy/2;
%  limysup=halfy+percentage_analysed*ncy/2;
%  limzinf=halfz-percentage_analysed*ncz/2;
%  limzsup=halfz+percentage_analysed*ncz/2;
%  conditionsx=Pos(:,1)<=limxinf|Pos(:,1)>=limxsup;
%  conditionsy=Pos(:,2)<=limyinf|Pos(:,2)>=limysup;
%  conditionsz=Pos(:,3)<=limzinf|Pos(:,3)>=limzsup;
%  conditions=conditionsx|conditionsy|conditionsz;
%  Pos(conditions,:)=[];


end

