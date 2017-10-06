function [ cell_bins1d_y cell_bins1d_z  count_sum Pos_halos mass_halos radius_halos count_sum_h_number count_sum_h_mass count_sum_dmdc] = proj2d_halos_data_out( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle)
%Computes the 2d projections aconding to the input specifications and stores (and returns) the resulting data

%   (example) [ cell_bins1d_y cell_bins1d_z  count_sum Pos_halos mass_halos radius_halos count_sum_h_number count_sum_h_mass count_sum_dmdc] = proj2d_halos_data_out( '/home/asus/Dropbox/extras/storage/guillimin/old/', '/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','','31.000halo0.dat',1,1,[0,0,0],[0,0])
%   (example) [ cell_bins1d_y cell_bins1d_z  count_sum Pos_halos mass_halos radius_halos count_sum_h_number count_sum_h_mass count_sum_dmdc] = proj2d_halos_data_out( '/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','','0.000halo0.dat',1,1,[0,0,0],[0,0])



cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos_h mass Radiusd halos] = preprocessing_halo_nodes( root,spec,aux_path,filename);

%to particle mass

mass=mass*(np/nc)^3;

cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
cell_bins1d_y(end)=[];
cell_bins1d_z(end)=[];

count_sum_h_number=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));
count_sum_h_mass=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));
count_sum=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));


theta=rot_angle(1);
phi=rot_angle(2);

Pos_halos=[];
mass_halos=[];
radius_halos=[];

path_in=strcat(root,spec,aux_path);


%for node = 1 : 1
for node = 1 : length(nodes_list)
    
    cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
    cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
    
    file_name = dir(strcat(path_in,num2str(z),'*halo',char(nodes_list(node)),'.dat'));
    filename=file_name.name;
    
    [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos_h mass Radiusd halos] = preprocessing_halo_nodes( root,spec,aux_path,filename);
    filename_pos=filename;
    [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,strcat(filename_pos(1:strfind(filename,'halo')-1),'xv',filename_pos(strfind(filename,'halo')+4:end)));
    
    %to particle mass
    
    mass=mass*(np/nc)^3;
    
    Pos_h=mod(Pos_h,nc);    
    Pos_h(1,:)=Pos_h(1,:)-(nc/2)-pivot(1);
    Pos_h(2,:)=Pos_h(2,:)-(nc/2)-pivot(2);
    Pos_h(3,:)=Pos_h(3,:)-(nc/2)-pivot(3);
    
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    Pos_h=Rz*Pos_h;
    Pos_h=Ry*Pos_h;
    %Pos=Rz*Pos;
    
    liminf=-(1/(2*lenght_factor))*nc;
    limsup= (1/(2*lenght_factor))*nc;
    conditionsx=Pos_h(1,:)<=liminf|Pos_h(1,:)>=limsup;
    conditionsy=Pos_h(2,:)<=liminf|Pos_h(2,:)>=limsup;
    conditionsz=Pos_h(3,:)<=liminf|Pos_h(3,:)>=limsup;
    conditions=conditionsx|conditionsy|conditionsz;
    Pos_h(:,conditions)=[];
    mass(conditions)=[];
    Radiusd(conditions)=[];
    
    Pos_h(1,:)=Pos_h(1,:)+(nc/2)+pivot(1);
    Pos_h(2,:)=Pos_h(2,:)+(nc/2)+pivot(2);
    Pos_h(3,:)=Pos_h(3,:)+(nc/2)+pivot(3);
        
    Pos_halos=[ Pos_halos , Pos_h];  
    mass_halos=[mass_halos , mass];
    radius_halos=[radius_halos , Radiusd];
        
    if (~isempty(Pos_h))
        
        
        [count_h_n edges mid loc] = histcn(transpose(Pos_h),1,cell_bins1d_y,cell_bins1d_z);
        count_h_n=count_h_n(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
        %     average=mean2(count);
        %     count=(count-average)/average;
        count_h_n=squeeze(count_h_n);
        
        %    cell_bins1d(end)=[];
        count_sum_h_number=count_sum_h_number+count_h_n;
        
        
        [count_h_m edges mid loc] = histcn(transpose(Pos_h),1,cell_bins1d_y,cell_bins1d_z,'AccumData',transpose(mass));
        count_h_m=count_h_m(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
        %     average=mean2(count);
        %     count=(count-average)/average;
        count_h_m=squeeze(count_h_m);
        
        %    cell_bins1d(end)=[];
        count_sum_h_mass=count_sum_h_mass+count_h_m;  
       
        
    end
    
    Pos=mod(Pos,nc);    
    Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
    Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
    Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
    
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    
    Pos=Rz*Pos;
    Pos=Ry*Pos;
    %Pos=Rz*Pos;
    
    liminf=-(1/(2*lenght_factor))*nc;
    limsup= (1/(2*lenght_factor))*nc;
    conditionsx=Pos(1,:)<=liminf|Pos(1,:)>=limsup;
    conditionsy=Pos(2,:)<=liminf|Pos(2,:)>=limsup;
    conditionsz=Pos(3,:)<=liminf|Pos(3,:)>=limsup;
    conditions=conditionsx|conditionsy|conditionsz;
    Pos(:,conditions)=[];
    
    Pos(1,:)=Pos(1,:)+(nc/2)+pivot(1);
    Pos(2,:)=Pos(2,:)+(nc/2)+pivot(2);
    Pos(3,:)=Pos(3,:)+(nc/2)+pivot(3);
    
    Pos=transpose(Pos);
    
    [count edges mid loc] = histcn(Pos,1,cell_bins1d_y,cell_bins1d_z);
    count=count(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
    %     average=mean2(count);
    %     count=(count-average)/average;
    count=squeeze(count);
    
    %    cell_bins1d(end)=[];
    count_sum=count_sum+count;
  
    
end

    if (isempty(Pos_h))
        
        Pos_halos=0;  
        mass_halos=0;
        radius_halos=0;
        
    end

%to particle mass unit

count_sum_h_mass=count_sum_h_mass;

average=mean2(count_sum);
count_sum_dmdc=(count_sum_h_mass-average)/average;

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));

path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/');
mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/'));

mkdir(path_out,'aux/total_dm/');
dlmwrite(strcat(path_out,'aux/total_dm/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_dm_z',num2str(z),'_data.txt'),count_sum,'delimiter','\t');
mkdir(path_out,strcat('number/'));
dlmwrite(strcat(path_out,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_n_z',num2str(z),'_data.txt'),count_sum_h_number,'delimiter','\t');
mkdir(path_out,'mass/');
dlmwrite(strcat(path_out,'mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_z',num2str(z),'_data.txt'),count_sum_h_mass,'delimiter','\t');
mkdir(path_out,'h_mass_dc_wrt_dmavr/');
dlmwrite(strcat(path_out,'h_mass_dc_wrt_dmavr/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_dmdc_z',num2str(z),'_data.txt'),count_sum_dmdc,'delimiter','\t');
mkdir(path_out,'aux/Pos_h/');
dlmwrite(strcat(path_out,'aux/Pos_h/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_Pos_halos_z',num2str(z),'_data.txt'),Pos_halos,'delimiter','\t');
mkdir(path_out,'aux/mass/');
dlmwrite(strcat(path_out,'aux/mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_mass_halos_z',num2str(z),'_data.txt'),mass_halos,'delimiter','\t');
mkdir(path_out,'aux/radius/');
dlmwrite(strcat(path_out,'aux/radius/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_radius_halos_z',num2str(z),'_data.txt'),radius_halos,'delimiter','\t');


cd('../2dproj');

end



