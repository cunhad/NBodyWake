function [ count_sum] = proj2d_cic_dm_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,data_stream)
%Computes the 2d projections aconding to the input specifications and stores (and returns) the resulting data

%   (example) [ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_cic_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',1,1,[0,0,0],[0,0],[1,2]);
%   (example) [ count_sum] = proj2d_cic_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,[0,0,0],[0,0,0],[1,2]);
%   (example) [ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_cic_dm_data_out('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/data/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','10.000xv0.dat',1,1,[0,0,0],[0,0],[1,2]);


% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% if specified, data will be stored in  root_out+spec+aux_path+aux_path_out

% filename is the output file from the nbody simulation

% lenght_factor = the analysis cube will have a lateral size given by the
% lateral size of the simulation cube divided by this number

% resol_factor= the bin will hte the particle bin size divided by this
%number

% pivot = a 3d array containing the translation wrt to the center of the
% cube (in grid cell units)

%rot_angle = 2d aray containing the theta and phy spherical angles pointing
%to the direction where the new z axis will be rotated


% data_stream=[0,1,2]
% if data_stream = 0, no output
% if data_stream = 1, binaries generated
% if data_stream = 2, text files generated

cd('../preprocessing');

% [ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );


% [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
[ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,length(nodes_list),1);

% display(z)
z_glob=z;

% cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
% cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
% cell_bins1d_y(end)=[];
% cell_bins1d_z(end)=[];
count_sum=zeros((np*resol_factor/lenght_factor),(np*resol_factor/lenght_factor));


psi=rot_angle(1);
theta=rot_angle(2);
phi=rot_angle(3);


%     cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
%     cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];

% display(cell_bins1d_y);
% display(cell_bins1d_z);

% display(length(nodes_list));

nb=np*resol_factor/lenght_factor;

% for node = 1 : 1
% for node = 1 : 16
for node = 1 : length(nodes_list)
    
%     cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
%     cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
        
    path_in=strcat(root,spec,aux_path);
    file_name = dir(strcat(path_in,num2str(z_glob),'*xv',char(nodes_list(node)),'.dat'));
    filename=file_name.name;
    
    display(filename);
    
    
    [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
%     [ ~,~,~,~,~,~,~,~,~,~,Pos ] = preprocessing_part(root,spec,aux_path,filename,length(nodes_list),node);

% display(length(Pos));   

     Pos=mod(Pos,nc);

    Pos=Pos*(np*resol_factor)/(nc);

    
    Pos(1,:)=Pos(1,:)-(np*resol_factor/2)-pivot(1)*(np*resol_factor)/(nc);
    Pos(2,:)=Pos(2,:)-(np*resol_factor/2)-pivot(2)*(np*resol_factor)/(nc);
    Pos(3,:)=Pos(3,:)-(np*resol_factor/2)-pivot(3)*(np*resol_factor)/(nc);
    
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    Rx = [ 1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
    
    Pos=Rz*Pos;
    Pos=Ry*Pos;
    Pos=Rx*Pos;
    
    
%         Pos(:,conditions)=[];
    
    Pos(1,:)=Pos(1,:)+(1/(2*lenght_factor))*np*resol_factor;
    Pos(2,:)=Pos(2,:)+(1/(2*lenght_factor))*np*resol_factor;
    Pos(3,:)=Pos(3,:)+(1/(2*lenght_factor))*np*resol_factor;
    
    
    
    
    
    Pos=mod(Pos,np*resol_factor);
    
        lim= (1/(lenght_factor))*np*resol_factor;
    conditionsx=Pos(1,:)<0|Pos(1,:)>lim;
    conditionsy=Pos(2,:)<0|Pos(2,:)>lim;
    conditionsz=Pos(3,:)<0|Pos(3,:)>lim;
    conditions=conditionsx|conditionsy|conditionsz;
    
    Pos(:,conditions)=[];
    
%     Pos(:,conditionsz)=[];

Pos=transpose(Pos);

count=zeros(nb,nb);

display(size(Pos,1))

if ~isempty(length(Pos))
    
    for particle=1:size(Pos,1)
        
        
        
        x=Pos(particle,1)-0.5;
        y=Pos(particle,2)-0.5;
        i1=floor(x)+1;
        j1=floor(y)+1;
        i2=i1+1;
        j2=j1+1;
        dx1=i1-x;
        dy1=j1-y;
        dx2=1-dx1;
        dy2=1-dy1;
        
        %         if (i1>=0 && i2< (np*resol_factor/lenght_factor)&& j1>=0 && j2< (np*resol_factor/lenght_factor))
        count(mod(i1,nb)+1,mod(j1,nb)+1)=count(mod(i1,nb)+1,mod(j1,nb)+1)+dx1*dy1;
        count(mod(i2,nb)+1,mod(j1,nb)+1)=count(mod(i2,nb)+1,mod(j1,nb)+1)+dx2*dy1;
        count(mod(i1,nb)+1,mod(j2,nb)+1)=count(mod(i1,nb)+1,mod(j2,nb)+1)+dx1*dy2;
        count(mod(i2,nb)+1,mod(j2,nb)+1)=count(mod(i2,nb)+1,mod(j2,nb)+1)+dx2*dy2;
        %         else
        %             display([i1,i2,j1,j2])
        %         end
        
    end
end
%     [count edges mid loc] = histcn(Pos,1,cell_bins1d_y,cell_bins1d_z);
%     count=count(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
%     %     average=mean2(count);
%     %     count=(count-average)/average;
%     count=squeeze(count);
% 
% % [count,Xedges,Yedges] = histcounts2(Pos(:,2),Pos(:,3),cell_bins1d_y,cell_bins1d_z);
% 
% % h = histogram2(Pos(2,:),Pos(3,:),cell_bins1d_y,cell_bins1d_z);


    %    cell_bins1d(end)=[];
    count_sum=count_sum+count;
    Pos=[];
%     clearvars count
end

% average=mean2(count_sum);
% count_sum=(count_sum-average)/average;

if ~ismember(0,data_stream)
    
    mkdir(root_out);
    mkdir(root_out,strcat(spec,aux_path));
    
    path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');
    mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/'));
    
    
    mkdir(path_out,'dc/');
    
   
    if ismember(1,data_stream)        
        fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data.bin'),'w');
        fwrite(fileID,count_sum, 'float32','l');
        fclose(fileID);
    end

    if ismember(2,data_stream)
        dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data.txt'),count_sum,'delimiter','\t');
    end
    
end


cd('../2dproj');

end



