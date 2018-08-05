function [ test_particle_displacement ] = displacement_evol_mem( root,root_data_out,root_plot_out,spec,aux_path )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (example) []=displacement_evol_mem( '/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/')


% mkdir(strcat(root_out))
mkdir(strcat(root_data_out,spec,aux_path,'check/displ/'))
mkdir(strcat(root_plot_out,spec,aux_path,'check/displ/'))

cd('../../../preprocessing')

% test_particle_id=10000;

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_spec( spec);
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

tnp=np^3;
% 
% filename_out=strcat(path_out,);
% fid = fopen(filename);

%create the empth files for the positions to be stored

for rds=1:length(redshift_list)

    filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
    fid = fopen(filename_out,'w');
    for partic=1:tnp
        fwrite(fid,[-1,-1,0],'float32','l');        
    end
    fclose(fid);
end

%stores the z positions of the particles for the no wake case

spec_nowake=strcat(string(size_box),'Mpc_',string(nc),'c_',string(np),'p_zi',string(zi),'_nowakem');


for rds=1:length(redshift_list)
    % for rds=5:5
    
%     pos_z=[];
%     part_id=[];
    
    number_node_dim=nthroot(numel(nodes_list), 3);
    
    for node=1:length(nodes_list)
        
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
        
        filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
        fid = fopen(filename_out,'r+');
        for partic=1:length(part_id)
            fseek(fid,(3*(part_id(partic)-1)*4),'bof');
            fwrite(fid,pos_z(partic),'float32','l');        
        end
        
    end
end


%stores the z positions of the particles for the wake case, and the
%displacement (with respect to the no wake case)

for rds=1:length(redshift_list)
    % for rds=5:5
    
%     pos_z=[];
%     part_id=[];
    
    number_node_dim=nthroot(numel(nodes_list), 3);
    
    for node=1:length(nodes_list)
        
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
        
        filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
        fid = fopen(filename_out,'r+');
        for partic=1:length(part_id)
%             fseek(fid,((3*(part_id(partic)-1)+1)*4),'bof');
%             fwrite(fid,pos_z(partic),'float32','l');        
            fseek(fid,((3*(part_id(partic)-1))*4),'bof');
            pos_nowake_this=fread(fid,1, 'float32','l');
            fwrite(fid,pos_z(partic),'float32','l');        
            if (pos_z(partic)-pos_nowake_this)>nc/2
            fwrite(fid,pos_z(partic)-pos_nowake_this-nc,'float32','l');        
            elseif  (pos_z(partic)-pos_nowake_this)<-nc/2
            fwrite(fid,pos_z(partic)-pos_nowake_this+nc,'float32','l');        
            else
            fwrite(fid,pos_z(partic)-pos_nowake_this,'float32','l');        
            end
        end
        
    end
end


% for rds=1:length(redshift_list)
%     filename_out=strcat(root_data_out,spec,aux_path,'check/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
%     pos_diff_ds = fileDatastore(filename_out,'ReadFcn',@read_bin,'FileExtensions','.dat');
%     pos_diff=cell2mat(tall(pos_diff_ds));
%     pos_diff=pos_diff(3,:);
%     pos_diff(pos_diff>nc/2)=pos_diff(pos_diff>nc/2)-nc;
%     pos_diff(pos_diff<-nc/2)=pos_diff(pos_diff<-nc/2)+nc;
% end

cd('../wake_detection/consist_check/Zeldovich/')


for rds=1:length(redshift_list)
    filename_out=strcat(root_data_out,spec,aux_path,'check/displ/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_pos_z',char(redshift_list(rds)),'.dat');
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

%     h_val=h.Values;
%     h_displ=h.BinEdges;
%     h_displ(end)=[];
%     hist_tot=[h_displ;h_val];
%     right_values=hist_tot(:,hist_tot(1,:)>0);
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



% wake_displacement=zeros(length(redshift_list),2);
% wake_displacement_zeld=zeros(length(redshift_list),1);
% test_particle_pos=zeros(length(redshift_list),1);
% test_particle_wpos=zeros(length(redshift_list),1);
% test_particle_displacement=zeros(length(redshift_list),1);
% 
% 
% 
% for rds=1:length(redshift_list)
%     % for rds=5:5
%     
%     pos_z=[];
%     part_id=[];
%     
%     number_node_dim=nthroot(numel(nodes_list), 3);
%     
%     for node=1:length(nodes_list)
%         
%         filename=strcat(root,spec,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
%         fid = fopen(filename);
%         fread(fid, [12 1], 'float32','l') ;
%         xv=fread(fid, [6 Inf], 'float32','l');
%         fclose(fid);
%         
%         node_ID=node-1;
%         k_node=floor(node_ID/number_node_dim^2);
%         
%         pos_z=[pos_z xv(3,:)+(nc/number_node_dim)*k_node];
%         
%         particle_ID_cat=strcat(root,spec,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
%         fid = fopen(particle_ID_cat);
%         fread(fid,6,'int64');
%         part_id=[part_id fread(fid,[1 Inf],'int64')];
%         fclose(fid);
%     end
%     
%     id_p_z=[part_id;pos_z];
%     w_sort_id_p_z=sort(transpose(id_p_z));
%     
%     test_particle_wpos(rds,1)=w_sort_id_p_z(test_particle_id,2);
%     
%     pos_z=[];
%     part_id=[];
%     
%     for node=1:length(nodes_list)
%         
%         spec_nowake=strcat(string(size_box),'Mpc_',string(nc),'c_',string(np),'p_zi',string(zi),'_nowakem');
%         
%         filename=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
%         fid = fopen(filename);
%         fread(fid, [12 1], 'float32','l') ;
%         xv=fread(fid, [6 Inf], 'float32','l');
%         fclose(fid);
%         
%         node_ID=node-1;
%         k_node=floor(node_ID/number_node_dim^2);
%         
%         pos_z=[pos_z xv(3,:)+(nc/number_node_dim)*k_node];
%         
%         particle_ID_cat=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
%         fid = fopen(particle_ID_cat);
%         fread(fid,6,'int64');
%         part_id=[part_id fread(fid,[1 Inf],'int64')];
%         fclose(fid);
%     end
%     
%     id_p_z=[part_id;pos_z];
%     nw_sort_id_p_z=sort(transpose(id_p_z));
%     
%     test_particle_pos(rds,1)=nw_sort_id_p_z(test_particle_id,2);
%     
%     displ=w_sort_id_p_z(:,2)-nw_sort_id_p_z(:,2);
%     avr=mean(displ);
%     displ(:)=displ(:)-avr;
%     
%     % fig=figure;
%     fig=figure('Visible', 'off');
%     hold on;
%     plot(nw_sort_id_p_z(:,2)*size_box/nc,displ*size_box/nc);
%     xlim ([-inf inf]);
%     xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%     ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
%     set(gca,'FontName','FixedWidth');
%     set(gca,'FontSize',16);
%     set(gca,'linewidth',2);
%     title({strcat('Displacement for z = '),char(redshift_list(rds))},'interpreter', 'latex', 'fontsize', 20);
%     hold off;
%     
%     
%     mkdir(strcat(root_plot_out,spec,aux_path,'check/'));
%     saveas(fig,strcat(root_plot_out,spec,aux_path,'check/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_z',char(redshift_list(rds)),'.png'));
%     
%     
%     positive_displ=displ(displ>=0);
%     negative_displ=displ(displ<0);
%     
%     % display(mean(positive_displ));
%     % display(std(positive_displ));
%     % display(mean(negative_displ));
%     % display(std(negative_displ));
%     
%     wake_displacement(rds,1)=(abs(mean(positive_displ))+abs(mean(negative_displ)))/2;
%     wake_displacement(rds,2)=(std(positive_displ)+std(negative_displ))/2;
%     
%     cd('../../parameters')
%     
%     [ ~, displacement, ~ ] = wake( Gmu,str2num(char(redshift_list(rds))));
%     wake_displacement_zeld(rds,1)=displacement;
%     
%     cd('../Analysis/preprocessing/')
%     
% end
% 
% % test_particle_displacement=test_particle_pos-test_particle_wpos;
% 
% cd('../wake_detection/consist_check/Zeldovich')
% 
% fig=figure('Visible', 'off');
% % fig=figure;
% errorbar((str2num(char(redshift_list))+1).^-1,wake_displacement(:,1)*size_box/nc,wake_displacement(:,2)*size_box/nc)
% hold on
% plot((str2num(char(redshift_list))+1).^-1,wake_displacement_zeld)
% 
% xlim ([-inf inf]);
% xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
%     ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
%     set(gca,'FontName','FixedWidth');
%     set(gca,'FontSize',16);
%     set(gca,'linewidth',2);
%     title({strcat('Displacement comparizon')},'interpreter', 'latex', 'fontsize', 20);
%     legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich")
%     hold off;
%     
%     saveas(fig,strcat(root_plot_out,spec,aux_path,'check/','_Check_Zel','.png'));
% figure;
% plot((str2num(char(redshift_list))+1).^-1,test_particle_displacement);
% 

% cd('../wake_detection/consist_check/Zeldovich/')
end
% 
% function data = read_bin( filename )
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% fid = fopen(filename);
% data=fread(fid, [3 Inf], 'float32','l');
% data=transpose(data);
% % data=data(3,:);
% fclose(fid);
% 
% end