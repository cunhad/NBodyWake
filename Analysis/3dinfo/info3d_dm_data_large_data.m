function [  ] = info3d_dm_data_large_data( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%NO%%%(example) [  dc_bins  count_sum] = info3d_dm_data('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','','0.000xv0.dat',1,1,[0,0,0],[0,0]);
%(example) [  dc_bins  count_sum] = info3d_dm_data_large_data('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','10.000xv0.dat',1,1,[0,0,0],[0,0]);
%NO%%%%(example)  [ Pos count ] = box_statistics_dm_analysis_dwt_snan('/home/asus/Dropbox/extras/storage/guillimin/', '/home/asus/Dropbox/extras/storage/guillimin/box_stat_cubic_fast/','/home/asus/Dropbox/extras/storage/guillimin/box_stat_cubic_fast_plot/','/home/as
%us/Dropbox/extras/storage/guillimin/box_stat_cubic_fast_snan/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','','','10.000xv0.dat',2,1,[0,0,0],256,16,4,'minmax',[1],[0,1,2,3],1,[1],'sym6',[1,2,3,4],[1,2,3],[0:1]);
% 
% % cd('../preprocessing');
% % 
% % [ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
% % 
% % [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% % 
% % cell_bins1d_x=[(nc/2)-(nc/(2*lenght_factor))+pivot(1):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(1)];
% % cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
% % cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
% % cell_bins1d_x(end)=[];
% % cell_bins1d_y(end)=[];
% % cell_bins1d_z(end)=[];
% % count_sum=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));
% % 
% % theta=rot_angle(1);
% % phi=rot_angle(2);
% % 
% % %for node = 1 : 1
% % 
% % %count particle in the box
% % 
% % num_part_total=0;
% % 
% % for node = 1 : length(nodes_list)
% % 
% %     
% %     path_in=strcat(root,spec,aux_path);
% %     file_name = dir(strcat(path_in,num2str(z),'*xv',char(nodes_list(node)),'.dat'));
% %     filename=file_name.name;
% %     
% %     [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% %         
% %     Pos=mod(Pos,nc);
% %     
% %     Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
% %     Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
% %     Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
% %     
% %     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% %     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
% %     Pos=Rz*Pos;
% %     Pos=Ry*Pos;
% %     %Pos=Rz*Pos;
% %     
% %     liminf=-(1/(2*lenght_factor))*nc;
% %     limsup= (1/(2*lenght_factor))*nc;
% %     conditionsx=Pos(1,:)<=liminf|Pos(1,:)>=limsup;
% %     conditionsy=Pos(2,:)<=liminf|Pos(2,:)>=limsup;
% %     conditionsz=Pos(3,:)<=liminf|Pos(3,:)>=limsup;
% %     conditions=conditionsx|conditionsy|conditionsz;
% %     Pos(:,conditions)=[];
% %     
% %     num_part=length(Pos(1,:));
% %     
% %     num_part_total=num_part_total+num_part;
% %     
% % end
% % 
% % aver_num_pat_total_per_cell=num_part_total/((np*resol_factor/lenght_factor)^3);
% % 
% % 
% % 
% % fig0=figure('Visible', 'off');
% % set(gcf, 'Position', [0 0 800 600]);
% % ax2 = axes('Position',[0.15 0.13 0.5 0.7]);
% % 
% % 
% % for node = 1 : length(nodes_list)
% %     
% %     cell_bins1d_x=[(nc/2)-(nc/(2*lenght_factor))+pivot(1):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(1)];
% %     cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
% %     cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
% %     
% %     
% %     path_in=strcat(root,spec,aux_path);
% %     file_name = dir(strcat(path_in,num2str(z),'*xv',char(nodes_list(node)),'.dat'));
% %     filename=file_name.name;
% %     
% %     [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% %     
% %     Pos=mod(Pos,nc);
% %     
% %     Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
% %     Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
% %     Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
% %     
% %     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% %     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
% %     Pos=Rz*Pos;
% %     Pos=Ry*Pos;
% %     %Pos=Rz*Pos;
% %     
% %     liminf=-(1/(2*lenght_factor))*nc;
% %     limsup= (1/(2*lenght_factor))*nc;
% %     conditionsx=Pos(1,:)<=liminf|Pos(1,:)>=limsup;
% %     conditionsy=Pos(2,:)<=liminf|Pos(2,:)>=limsup;
% %     conditionsz=Pos(3,:)<=liminf|Pos(3,:)>=limsup;
% %     conditions=conditionsx|conditionsy|conditionsz;
% %     Pos(:,conditions)=[];
% %     
% %     Pos(1,:)=Pos(1,:)+(nc/2)+pivot(1);
% %     Pos(2,:)=Pos(2,:)+(nc/2)+pivot(2);
% %     Pos(3,:)=Pos(3,:)+(nc/2)+pivot(3);
% %     
% %     Pos=transpose(Pos);
% %     
% %     
% %     [count edges mid loc] = histcn(Pos,cell_bins1d_x,cell_bins1d_y,cell_bins1d_z);
% %     count=count(2:numel(cell_bins1d_x)-1,2:numel(cell_bins1d_y)-1,2:numel(cell_bins1d_z)-1);
% % %     %     average=mean2(count);
% % %     %     count=(count-average)/average;
% % %     count=squeeze(count);
% % %     
% % %     %    cell_bins1d(end)=[];
% % %     count_sum=count_sum+count;
% % end
% % 
% % % average=mean2(count_sum);
% % % count_sum=(count_sum-average)/average;
% % 
% % % mkdir(root_out);
% % % mkdir(root_out,strcat(spec,aux_path));
% % % 
% % % path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');
% % % mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/'));
% % % 
% % % 
% % % mkdir(path_out,'dc/');
% % % 
% % % dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_data.txt'),count_sum,'delimiter','\t');
% % 
% % 
% % % scatter3(pos(:,1),pos(:,2),pos(:,3));
% 
% cd('../3dinfo');
tic;

cd('../preprocessing');
% [~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );
% [xv_files_list,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );
[xv_files_list,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );



cd('../3dinfo');

toc;

end

