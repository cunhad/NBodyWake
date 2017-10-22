function [  ] = box_statistics_dm_from_nodes_to_total( root_per_node_out,spec,aux_path,aux_path_per_node_out,NSIDE )

%(example) box_statistics_dm_from_nodes_to_total('/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','',4);
%(example) box_statistics_dm_from_nodes_to_total('/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','',4);

pivot=[0,0,0]; %this is the position od the origin of the rotation point with respect to the center of the box
lenght_factor=2;
resol_factor=2;

path_per_node_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/','npart_per_node_1dproj/');
   


% path_out=strcat(path,spec,aux_path,strcat('Analysis/','stat/box_statistics/'));
%for guillimin
% path_out=strcat('/gs/scratch/cunhad/',spec,aux_path);


files_list = dir(strcat(path_per_node_out,'*','0_NSIDE',num2str(NSIDE),'.txt'));
%files_list = dir(strcat(path_in,'*'));
sorted_files_list={files_list.name};

cd('../../processing');

angles = dlmread(strcat('../../python/angles',num2str(NSIDE),'.txt'));
[angle_nuple,number_of_angle_nuple] = size(angles);

        theta=angles(1,:);
        phi=angles(2,:);
        
sorted_files_list=sort_nat(sorted_files_list);

[aux1 aux2] = size(num2str(NSIDE));
aux3=aux2+16;
redshift_list=cellfun(@(x) x(15:end-aux3),sorted_files_list,'UniformOutput', false);

files_list2 = dir(strcat(path_per_node_out,'1dproj_angle_z',char(redshift_list(1)),'_node','*_NSIDE',num2str(NSIDE),'.txt'));
sorted_files_list2={files_list2.name};
nodes_list=cellfun(@(x) x(4+cell2mat(strfind(sorted_files_list2(1:1), 'node')):-2+cell2mat(strfind(sorted_files_list2(1:1), 'NSIDE'))),sorted_files_list2,'UniformOutput', false);

nodes_list=sort_nat(nodes_list);




for rds = 1 : length(redshift_list)  
    

    display(char(strcat(path_per_node_out,'1dproj_angle_z',char(redshift_list(1)),'_node',char(nodes_list(1)),'_NSIDE',num2str(NSIDE),'.txt')));
    [rows columns] = size(dlmread(char(strcat(path_per_node_out,'1dproj_angle_z',char(redshift_list(1)),'_node',char(nodes_list(1)),'_NSIDE',num2str(NSIDE),'.txt'))));
    
    display(rows);
    display(columns);
    
    count=zeros(rows,columns);

    for node = 1 : length(nodes_list)

        filename=char(strcat(path_per_node_out,'1dproj_angle_z',char(redshift_list(rds)),'_node',char(nodes_list(node)),'_NSIDE',num2str(NSIDE),'.txt'));
        display(filename);
        
        [rows columns] = size(dlmread(filename));
       % display(rows)
       % display(columns);
        
     count = dlmread(filename)+count;
     
    end
   
    

    
%     

%     path_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
% 
%     mkdir(path_out,'npart_all_nodes/1d_proj/');
%     dlmwrite(strcat(path_out,'npart_all_nodes/1d_proj/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_npart_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'),count,'delimiter','\t');

    
    count=transpose(count);
% %     
    average=mean(count);
 %   average=repmat(average,columns,1);
    count(:,:)=(count(:,:)-average(1,:))./average(1,:);

    count=transpose(count);

    average=transpose(average);    
    
     path_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    
    mkdir(path_out,'dc_all_nodes_1dproj/');
   % dlmwrite(strcat(path_out,'dc_all_nodes_1dproj/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_dc_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'),count,'delimiter','\t');
     dlmwrite(strcat(path_out,'dc_all_nodes_1dproj/','_1dproj_dc_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'),count,'delimiter','\t');

     mkdir(path_out,'dc_all_nodes_1dproj/avr/');
   % dlmwrite(strcat(path_out,'dc_all_nodes_1dproj/avr/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_avr_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'),average,'delimiter','\t');
     dlmwrite(strcat(path_out,'dc_all_nodes_1dproj/avr/','_1dproj_avr_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'),average,'delimiter','\t');


    


    
end

cd('../wake_detection/box_statistics');

end

