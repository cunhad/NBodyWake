function [  ] = preprocessing_nodes_to_csv( path,spec,aux_path,filename)
%   This function takes the phase space output from CUBEP3M and returns
%   some relevant information plus the global positions of all particles at the
%   corresponding redshift and node volume

%(example) [  ] = preprocessing_nodes_to_csv('/home/asus/Dropbox/extras/storage/graham/small_res/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/','10.000xv0.dat')


[ ~,~,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info( path,spec,aux_path)

files_node_list=dir(strcat(path,spec,aux_path,filename(1:end-5),'*','.dat'));
files_node_list={files_node_list.name};

% display(files_node_list)

% fid_o = fopen(strcat(path,spec,aux_path,filename(1:end-7),'positions.csv'),'w');

tic

parfor file=1:length(files_node_list)
% for file=1:1
    
    [ i_node,j_node,k_node,number_node_dim,~ ] = preprocessing_filename_info( path,spec,aux_path,files_node_list{file})
    
    fid = fopen(strcat(path,spec,aux_path,files_node_list{file}));
    directory = dir(strcat(path,spec,aux_path,files_node_list{file}));
    particles=(directory.bytes-48)/24;
    fseek(fid,48,'bof');
%     for part=1:particles
    pos=fread(fid, [6 particles], 'float32','l');
    pos(1,:)=pos(1,:)+(nc/number_node_dim)*i_node;
    pos(2,:)=pos(2,:)+(nc/number_node_dim)*j_node;
    pos(3,:)=pos(3,:)+(nc/number_node_dim)*k_node;
    dlmwrite(strcat(path,spec,aux_path,filename(1:end-7),'positions.csv'),pos(1:3,:)','-append','newline','unix','precision','%.10f') ;
%     end
    fclose(fid);
end

toc


%  
%  
% fid = fopen(path_file_in);
% directory = dir(path_file_in);
% particles=(directory.bytes-48)/24;
% headear = fread(fid, [12 1], 'float32','l') ;
% data=fread(fid, [6 particles], 'float32','l');
% Pos = data(1:3,:);
% %Pos=transpose(Pos);
% %Pos=mod(Pos,nc);
%  
% %in this part we will get the position of the wake taking into acount the
% %node structure
% 
%  node=char(filename);
%  node=str2num(node(strfind(filename, 'xv')+2:strfind(filename,'.dat')-1));
%  [ nodes_list redshift_list ] = preprocessing_many_nodes(path,spec,aux_path);
%  
%  number_node_dim=nthroot(numel(nodes_list), 3);
%  k_node=floor(node/number_node_dim^2);
%  res=mod(node,number_node_dim^2);
%  j_node=floor(res/number_node_dim);
%  i_node=mod(res,number_node_dim);
%  
% %  display(numel(nodes_list));
% %  display(number_node_dim);
% %  display(node);
% %  display(i_node);
% %  display(j_node);
% %  display(k_node);
%  
%  Pos(1,:)=Pos(1,:)+(nc/number_node_dim)*i_node;
%  Pos(2,:)=Pos(2,:)+(nc/number_node_dim)*j_node;
%  Pos(3,:)=Pos(3,:)+(nc/number_node_dim)*k_node;
%  
%  
% %  XM=nc;
% %  Xm=0;
% %  YM=nc;
% %  Ym=0;
% %  ZM=nc;
% %  Zm=0;
% %  ncx=nc;
% %  ncy=nc;
% %  ncz=nc;
% %  halfx=(XM+Xm)/2;
% %  halfy=(YM+Ym)/2;
% %  halfz=(ZM+Zm)/2;
% %  limxinf=halfx-percentage_analysed*ncx/2;
% %  limxsup=halfx+percentage_analysed*ncx/2;
% %  limyinf=halfy-percentage_analysed*ncy/2;
% %  limysup=halfy+percentage_analysed*ncy/2;
% %  limzinf=halfz-percentage_analysed*ncz/2;
% %  limzsup=halfz+percentage_analysed*ncz/2;
% %  conditionsx=Pos(:,1)<=limxinf|Pos(:,1)>=limxsup;
% %  conditionsy=Pos(:,2)<=limyinf|Pos(:,2)>=limysup;
% %  conditionsz=Pos(:,3)<=limzinf|Pos(:,3)>=limzsup;
% %  conditions=conditionsx|conditionsy|conditionsz;
% %  Pos(conditions,:)=[];
% 

end

