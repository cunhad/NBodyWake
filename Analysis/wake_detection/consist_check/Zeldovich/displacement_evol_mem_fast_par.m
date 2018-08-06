function [  ] = displacement_evol_mem_fast_par( root,root_data_out,root_plot_out,spec,aux_path,num_cores)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (example) []=displacement_evol_mem_fast_par( '/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/',4)


myCluster = parcluster('local');
myCluster.NumWorkers=num_cores;
saveProfile(myCluster);

p = parpool(num_cores);

% mkdir(strcat(root_out))
mkdir(strcat(root_data_out,spec,aux_path,'check/displ/'))
mkdir(strcat(root_plot_out,spec,aux_path,'check/displ/'))

cd('../../../preprocessing')

% test_particle_id=10000;

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_spec( spec);
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

tnp=np^3;

% for rds=1:length(redshift_list)
%    
%     displacement_info=tall(zeros(tnp));
%     
% end
% 
% create the empth files for the positions to be stored
% 
% for rds=1:length(redshift_list)
% 
%     filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
%     fid = fopen(filename_out,'w');
%     for node=1:length(nodes_list)
%         fwrite(fid,zeros(3,tnp/length(nodes_list)),'float32','l');        
%     end
%     fclose(fid);
% end

%stores the z positions of the particles for the no wake case

spec_nowake=strcat(string(size_box),'Mpc_',string(nc),'c_',string(np),'p_zi',string(zi),'_nowakem');

% for rds=1:length(redshift_list)
% 
%     filename_xv=cellstr(strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list),'.dat'));
%     xv_ds = fileDatastore(filename_xv,'ReadFcn',@read_bin,'FileExtensions','.dat');
%     xv=cell2mat(tall(xv_ds));
%     
%     filename_pid=cellstr(strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'PID',char(nodes_list),'.dat'));
%     pid_ds = fileDatastore(filename_pid,'ReadFcn',@read_pid,'FileExtensions','.dat');
%     pid=cell2mat(tall(pid_ds));
%     
%     pid_posZ_nowake=[pid;xv(:,3)];
%     
% end


for rds=1:length(redshift_list)
    % for rds=5:5
    
%     pos_z=[];
%     part_id=[];
    
    number_node_dim=nthroot(numel(nodes_list), 3);
    
    parfor node=1:length(nodes_list)

	filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list(node)),'_Check_Zel_wpid_posZ_displ_z',char(redshift_list(rds)),'.dat');
    fid_o = fopen(filename_out,'w');
        
        filename=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=xv(3,:)+(nc/number_node_dim)*k_node;
        
        particle_ID_cat=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=fread(fid,[1 Inf],'int64');
        fclose(fid);
        
        [sorted_ sort_inx]=sort(part_id);
        sorted_pos_z=pos_z(sort_inx);
%         pid_posZ=[part_id(sort_inx);pos_z(sort_inx)];                
        
        
        %now with wake
        
        filename=strcat(root,spec,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=xv(3,:)+(nc/number_node_dim)*k_node;
        
        particle_ID_cat=strcat(root,spec,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=fread(fid,[1 Inf],'int64');
        fclose(fid);

        [sorted_w_ sort_inx]=sort(part_id);
        sorted_pos_z_w=pos_z(sort_inx);
%         pid_posZ_w=[part_id(sort_inx);pos_z(sort_inx)];  
        
        
%         find(sorted_w_(:)==sorted_
        
%         displacement=pid_posZ_w(:,2)-pid_posZ(find(sorted_w_(:)==sorted_),2);

%         lent_out=length(find(sorted_w_(:)==sorted_));
        [comom_ID,com_id_w,com_id_nw]=intersect(sorted_w_,sorted_);
        part_z_nw=sorted_pos_z(com_id_nw);
        part_z_w=sorted_pos_z_w(com_id_w);
        displac_list=part_z_w-part_z_nw;
        displac_list(displac_list>nc/2)=displac_list(displac_list>nc/2)-nc;
        displac_list(displac_list<-nc/2)=displac_list(displac_list<-nc/2)+nc;
        
        fwrite(fid_o,[comom_ID;part_z_w;displac_list],'float32','l');
    fclose(fid_o);        
    end

end





cd('../wake_detection/consist_check/Zeldovich/')


for rds=1:length(redshift_list)
    filename_out=cellstr(strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list),'_Check_Zel_wpid_posZ_displ_z',char(redshift_list(rds)),'.dat'));
    pos_diff_ds = fileDatastore(filename_out,'ReadFcn',@read_bin,'FileExtensions','.dat');
    pos_diff=cell2mat(tall(pos_diff_ds));
    
    fig=figure('Visible', 'off');
    histogram2(mod(pos_diff(:,2),nc)*size_box/nc,pos_diff(:,3)*size_box/nc,nc/2,'DisplayStyle','tile','ShowEmptyBins','on');
    mkdir(strcat(root_plot_out,spec,aux_path,'check/displ/'));
    saveas(fig,strcat(root_plot_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_displacement_z',char(redshift_list(rds)),'_plot.png'));
    
    fig=figure('Visible', 'off');
    h=histogram(pos_diff(:,3)*size_box/nc);
    mkdir(strcat(root_plot_out,spec,aux_path,'check/displ/hist/'));    
    saveas(fig,strcat(root_plot_out,spec,aux_path,'check/displ/hist/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_displacement_z',char(redshift_list(rds)),'_plot.png'));

    posit_values=pos_diff(pos_diff(:,3)>0,3);
    mn_pos(rds)=gather(mean(posit_values));
    std_pos(rds)=gather(std(posit_values,1));
    
    negat_values=pos_diff(pos_diff(:,3)<0,3);
    mn_neg(rds)=gather(mean(negat_values));
    std_neg(rds)=gather(std(negat_values,1));
    
end

cd('../../../../parameters')
    
for rds=1:length(redshift_list)    
    [ ~, displacement, ~ ] = wake( Gmu,str2num(char(redshift_list(rds))));
    wake_displacement_zeld(rds,1)=displacement;
    
end

cd('../Analysis/wake_detection/consist_check/Zeldovich/')


%plot positive values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,mn_pos*size_box/nc,std_pos*size_box/nc)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_displacement_zeld)

xlim ([-inf inf]);
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Displacement comparizon: positive')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich")
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,'check/displ/','_Check_Zel_pos','.png'));

%plot negative values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,abs(mn_neg)*size_box/nc,std_neg*size_box/nc)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_displacement_zeld)

xlim ([-inf inf]);
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Displacement comparizon: negative')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich")
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,'check/displ/','_Check_Zel_neg','.png'));


%plot total values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,((mn_pos+abs(mn_neg))/2)*size_box/nc,((std_pos+std_neg)/2)*size_box/nc)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_displacement_zeld)

xlim ([-inf inf]);
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Displacement comparizon')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich")
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,'check/displ/','_Check_Zel','.png'));

delete(gcp('nocreate'))

end