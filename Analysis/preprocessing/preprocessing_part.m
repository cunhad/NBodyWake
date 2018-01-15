function [ size_box, nc, np, zi, wake_or_no_wake ,multiplicity_of_files ,Gmu ,ziw ,z, path_file_in, Pos  ] = preprocessing_part( path,spec,aux_path,filename,part,part_id)
%   This function takes the phase space output from CUBEP3M and returns
%   some relevant information plus the global positions of all particles at the
%   corresponding redshift and node volume


%(example) [ size_box, nc, np, zi, wake_or_no_wake ,multiplicity_of_files ,Gmu ,ziw ,z, path_file_in, Pos ] = preprocessing_part( '/home/asus/Dropbox/extras/storage/graham/small_res/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','0.000xv0.dat',8,1);



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
 
 z_string=char(filename);
 z_string=z_string(1:end-7);
 z=str2num(z_string);
 
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
 
%look at the chunks
 
files_node_list=dir(strcat(path,spec,aux_path,num2str(z_string),'xv*','.dat'));
files_node_list={files_node_list.name};

particle_number_per_file=ones(1,length(files_node_list)+1);
particle_sum_per_file=zeros(1,length(files_node_list)+1);
particle_sum=0;
for node_ID=1:length(files_node_list)
    directory = dir(strcat(path,spec,aux_path,files_node_list{node_ID}));
    particle_number_per_file(node_ID+1)=(directory.bytes-48)/24;
    particle_sum = (directory.bytes-48)/24 + particle_sum;
    particle_sum_per_file(node_ID+1)=particle_sum;
end

total_number_of_particles=sum(particle_number_per_file)-1;
size_of_particle_chunk=ceil(total_number_of_particles/part);

particle_sum_per_chunk=ones(1,part+1);
for part_ID=1:part
    particle_sum_per_chunk(part_ID+1)=size_of_particle_chunk*part_ID;
end
particle_sum_per_chunk(part+1)=total_number_of_particles;

part_chunck_start=particle_sum_per_chunk(part_id);
part_chunck_end=particle_sum_per_chunk(part_id+1);


for node_ID=1:length(files_node_list)
    if ((part_chunck_start>particle_sum_per_file(node_ID))&&(part_chunck_start<=particle_sum_per_file(node_ID+1)))
        part_files_start_ID=node_ID;
    end
    if ((part_chunck_end>particle_sum_per_file(node_ID))&&(part_chunck_end<=particle_sum_per_file(node_ID+1)))
        part_files_end_ID=node_ID;
    end
end

% data=[];
% for node_ID=part_files_start_ID:part_files_end_ID
%     fid = fopen(strcat(path,spec,aux_path,files_node_list{node_ID}));
%     fread(fid, [12 1], 'float32','l') ;
%     if (part_chunck_start>particle_sum_per_file(node_ID))
%         for i=1:part_chunck_start-particle_sum_per_file(node_ID)-1
%             fread(fid, [6 1], 'float32','l');                            
%         end
%     end
%     if (part_chunck_end<=particle_sum_per_file(node_ID+1))
%         data=[data fread(fid, [6 size_of_particle_chunk], 'float32','l')];
%     else 
%         data=[data fread(fid, [6 Inf], 'float32','l')];
%     end   
% end

[ nodes_list redshift_list ] = preprocessing_many_nodes(path,spec,aux_path);
number_node_dim=nthroot(numel(nodes_list), 3);

data=[];
for node_ID=part_files_start_ID:part_files_end_ID
    fid = fopen(strcat(path,spec,aux_path,files_node_list{node_ID}));
    fread(fid, [12 1], 'float32','l') ;
    if (part_chunck_start>particle_sum_per_file(node_ID))
        for i=1:part_chunck_start-particle_sum_per_file(node_ID)-1
            fread(fid, [6 1], 'float32','l');                            
        end
    end
    if (part_chunck_end<=particle_sum_per_file(node_ID+1))
        data_file=fread(fid, [6 size_of_particle_chunk], 'float32','l');
    else 
        data_file=fread(fid, [6 Inf], 'float32','l');
    end
    node=node_ID-1;
    k_node=floor(node/number_node_dim^2);
    res=mod(node,number_node_dim^2);
    j_node=floor(res/number_node_dim);
    i_node=mod(res,number_node_dim);
    data_file(1,:)=data_file(1,:)+(nc/number_node_dim)*i_node;
    data_file(2,:)=data_file(2,:)+(nc/number_node_dim)*j_node;
    data_file(3,:)=data_file(3,:)+(nc/number_node_dim)*k_node;
    data=[data data_file];
end


Pos = data(1:3,:);

%Pos=transpose(Pos);
%Pos=mod(Pos,nc);
 
%in this part we will get the position of the wake taking into acount the
%node structure


%  k_node=floor(node/number_node_dim^2);
%  res=mod(node,number_node_dim^2);
%  j_node=floor(res/number_node_dim);
%  i_node=mod(res,number_node_dim);
%  
%  
%  Pos(1,:)=Pos(1,:)+(nc/number_node_dim)*i_node;
%  Pos(2,:)=Pos(2,:)+(nc/number_node_dim)*j_node;
%  Pos(3,:)=Pos(3,:)+(nc/number_node_dim)*k_node;
 
 

end

