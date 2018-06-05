function [ ] = proj1d_dm_data_out( path_in,file_in,path_out,bins)
% %Computes the 2d projections aconding to the input specifications and stores (and returns) the resulting data

% (example) proj1d_dm_data_out('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/','snapshot','/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_wakeGmu1t10m7zi31m/sample0001/gadget_out/',32);
% (example) proj1d_dm_data_out('/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/','out','/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_wakeGmu1t10m7zi31m/sample0001/picola_out/',32);

% path_in='/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/picola_out/';
% file_in='out';
% path_out='/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/plot/32Mpc_64c_64p_zi63_wakeGmu1t10m7zi31m/sample0001/picola_out/';

% path_in='/home/asus/Dropbox/extras/storage/laptop/simulations_gadget/32Mpc_64c_64p_zi63_nowakem/sample0001/gadget_out/';
% file_in='	';

mkdir(path_out); 

cd('../processing');

xv_files_list_head=dir(strcat(path_in,file_in,'*.0'));
if ~isempty(xv_files_list_head)
xv_files_list_head={xv_files_list_head.name};
    from_gadget=false;

else
    xv_files_list_head=dir(strcat(path_in,file_in,'*'));
    xv_files_list_head={xv_files_list_head.name};
    from_gadget=true;

end

for file_z = 1 : length(xv_files_list_head)
   
    cd('../processing');

    file_name_head=char(xv_files_list_head(file_z));
    if ~from_gadget
        file_name_head=file_name_head(1:end-2);
    else
        file_name_head=file_name_head;
    end
        
    xv_files_list=dir(strcat(path_in,file_name_head,'*'));
    xv_files_list={xv_files_list.name};
    xv_files_list=sort_nat(xv_files_list);
    
    count_sum=zeros(bins-1,bins-1);
    cd('../preprocessing');
    
    
    for out = 1 : length(xv_files_list)
        
        file_i=char(xv_files_list(out));
        node=file_i(strfind(file_i,'.')+1:end);
        filename_in=strcat(path_in,file_i);
        
        %read header

        fid = fopen(strcat(filename_in));
        block_size_1_head=fread(fid, 1, 'uint','l') ;
        Npart=fread(fid, [6 1], 'uint','l') ;
        Massarr=fread(fid, [6 1], 'double','l') ;
        Time=fread(fid,1, 'double','l') ;
        Redshift=fread(fid,1, 'double','l') ;
        FlagSfr=fread(fid, 1, 'int','l') ;
        FlagFeedback=fread(fid, 1, 'int','l') ;
        Nall=fread(fid, [6 1], 'int','l') ;
        FlagCooling=fread(fid, 1, 'int','l') ;
        NumFiles=fread(fid, 1, 'int','l') ;
        BoxSize=fread(fid,1, 'double','l') ;
        Omega0=fread(fid,1, 'double','l') ;
        OmegaLambda=fread(fid,1, 'double','l') ;
        HubbleParam=fread(fid,1, 'double','l') ;
        FlagAge=fread(fid, 1, 'int','l') ;
        FlagMetals=fread(fid, 1, 'int','l') ;
        NallHW=fread(fid, [6 1], 'int','l') ;
        flag_entr_ics=fread(fid, 1, 'int','l') ;
        garbage=fread(fid, [15 1], 'int','l') ;
        block_size_1_tail=fread(fid, 1, 'uint','l') ;
        
        % read positions
        
        block_size_2_head=fread(fid, 1, 'uint','l') ;
        pos_0=fread(fid, [3 Npart(1+0)], 'single','l') ;
        pos_1=fread(fid, [3 Npart(1+1)], 'single','l') ;
        pos_2=fread(fid, [3 Npart(1+2)], 'single','l') ;
        pos_3=fread(fid, [3 Npart(1+3)], 'single','l') ;
        pos_4=fread(fid, [3 Npart(1+4)], 'single','l') ;
        pos_5=fread(fid, [3 Npart(1+5)], 'single','l') ;
        block_size_2_tail=fread(fid, 1, 'uint','l') ;

          pos_1=transpose(mod(pos_1,BoxSize));
    
    
    [count edges mid loc] = histcn(pos_1,1,bins,bins);
    count=count(1:1,1:1,1:bins-1);
    count=squeeze(count);
    count=squeeze(count);
    count_sum=count_sum+count;
        
    
        
    end
    
    
    fig1=figure('Visible', 'off');
    set(gcf, 'Position', [0 0 600 600]);
    plot([1:bins-1]*BoxSize/bins,(count_sum-mean(count_sum(:)))/mean(count_sum(:)));  
    ylim([-1 2]);
    xlim([0 BoxSize]);
    xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('Density contrast', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title({strcat('Density contrast of the 1d projection'),'of dark matter mass ar z=',Redshift},'interpreter', 'latex', 'fontsize', 20);
    hold off;
    saveas(fig1,strcat(path_out,'1d_lim_dc_',file_name_head,'.png'));

    fig2=figure('Visible', 'off');
    set(gcf, 'Position', [0 0 600 600]);
    plot([1:bins-1]*BoxSize/bins,(count_sum-mean(count_sum(:)))/mean(count_sum(:)));
    hold on;    
    xlim([0 BoxSize]);
    xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('Density contrast', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title({strcat('Density contrast of the 1d projection'),'of dark matter mass at z=',Redshift},'interpreter', 'latex', 'fontsize', 20);
    hold off;
    saveas(fig2,strcat(path_out,'1d_dc',file_name_head,'.png'));
    
end




 

% %   (example) [ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','0.000xv0.dat',1,1,[0,0,0],[0,0],[1,2]);
% %   (example) [ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_dm_data_out('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','63.000xv0.dat',1,1,[0,0,0],[0,0],[1,2]);
% %   (example) [ cell_bins1d_y,cell_bins1d_z,count_sum] = proj2d_dm_data_out('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/data/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','10.000xv0.dat',1,1,[0,0,0],[0,0],[1,2]);
% 
% 
% % NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )
% 
% % if specified, data will be stored in  root_out+spec+aux_path+aux_path_out
% 
% % filename is the output file from the nbody simulation
% 
% % lenght_factor = the analysis cube will have a lateral size given by the
% % lateral size of the simulation cube divided by this number
% 
% % resol_factor= the bin will hte the particle bin size divided by this
% %number
% 
% % pivot = a 3d array containing the translation wrt to the center of the
% % cube (in grid cell units)
% 
% %rot_angle = 2d aray containing the theta and phy spherical angles pointing
% %to the direction where the new z axis will be rotated
% 
% 
% % data_stream=[0,1,2]
% % if data_stream = 0, no output
% % if data_stream = 1, binaries generated
% % if data_stream = 2, text files generated
% 
% cd('../preprocessing');
% 
% % [ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
% [~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );
% 
% 
% % [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,length(nodes_list),1);
% 
% % display(z)
% z_glob=z;
% 
% cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
% cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
% cell_bins1d_y(end)=[];
% cell_bins1d_z(end)=[];
% count_sum=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));
% 
% theta=rot_angle(1);
% phi=rot_angle(2);
% 
%     cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
%     cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
% 
% % display(cell_bins1d_y);
% % display(cell_bins1d_z);
% 
% % display(length(nodes_list));
% 
% % for node = 1 : 64
%  for node = 1 : length(nodes_list)
%     
% %     cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];
% %     cell_bins1d_z=[(nc/2)-(nc/(2*lenght_factor))+pivot(3):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(3)];
%         
%     path_in=strcat(root,spec,aux_path);
%     file_name = dir(strcat(path_in,num2str(z_glob),'*xv',char(nodes_list(node)),'.dat'));
%     filename=file_name.name;
%     
%     display(filename);
%     
%     [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename);
% %     [ ~,~,~,~,~,~,~,~,~,~,Pos ] = preprocessing_part(root,spec,aux_path,filename,length(nodes_list),node);
% 
% % display(length(Pos));   
% 
%     Pos=mod(Pos,nc);
%     
%     Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
%     Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
%     Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
%     
%     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
%     Pos=Rz*Pos;
%     Pos=Ry*Pos;
%     %Pos=Rz*Pos;
%     
%     liminf=-(1/(2*lenght_factor))*nc;
%     limsup= (1/(2*lenght_factor))*nc;
%     conditionsx=Pos(1,:)<=liminf|Pos(1,:)>=limsup;
%     conditionsy=Pos(2,:)<=liminf|Pos(2,:)>=limsup;
%     conditionsz=Pos(3,:)<=liminf|Pos(3,:)>=limsup;
%     conditions=conditionsx|conditionsy|conditionsz;
%     Pos(:,conditions)=[];
%     
%     Pos(1,:)=Pos(1,:)+(nc/2)+pivot(1);
%     Pos(2,:)=Pos(2,:)+(nc/2)+pivot(2);
%     Pos(3,:)=Pos(3,:)+(nc/2)+pivot(3);
%     
%     Pos=transpose(Pos);
%     
%     
%     [count edges mid loc] = histcn(Pos,1,cell_bins1d_y,cell_bins1d_z);
%     count=count(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
%     %     average=mean2(count);
%     %     count=(count-average)/average;
%     count=squeeze(count);
% 
% % [count,Xedges,Yedges] = histcounts2(Pos(:,2),Pos(:,3),cell_bins1d_y,cell_bins1d_z);
% 
% % h = histogram2(Pos(2,:),Pos(3,:),cell_bins1d_y,cell_bins1d_z);
% 
% 
%     %    cell_bins1d(end)=[];
%     count_sum=count_sum+count;
% %     clearvars count
% end
% 
% % average=mean2(count_sum);
% % count_sum=(count_sum-average)/average;
% 
% if ~ismember(0,data_stream)
%     
%     mkdir(root_out);
%     mkdir(root_out,strcat(spec,aux_path));
%     
%     path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');
%     mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/'));
%     
%     
%     mkdir(path_out,'dc/');
%     
%    
%     if ismember(1,data_stream)        
%         fileID = fopen(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data.bin'),'w');
%         fwrite(fileID,count_sum, 'float32','l');
%         fclose(fileID);
%     end
% 
%     if ismember(2,data_stream)
%         dlmwrite(strcat(path_out,'dc/','_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data.txt'),count_sum,'delimiter','\t');
%     end
%     
% end
% 
% 
% cd('../2dproj');
% 

cd('../1dproj');

end



