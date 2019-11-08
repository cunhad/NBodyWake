function [  ] = comvel_evol_mem_fast_meanslices_par( root,root_data_out,root_plot_out,spec,aux_path,wake_type,num_cores)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (example) []=comvel_evol_mem_fast_par( '/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_wakeGmu6t10m6zi10m','/sample1001/','test/',4)
% (example) []=comvel_evol_mem_fast_meanslices_par( '/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_wakeGmu5t10m5zi63m','/sample1001/','',4)



type_folder=wake_type;

myCluster = parcluster('local');
myCluster.NumWorkers=num_cores;
saveProfile(myCluster);

p = parpool(num_cores);

% if wake_type==0
%    type_folder='' 
% end
% 
% if wake_type==1
%    type_folder='test/' 
% end
% 
% if wake_type==2
%    type_folder='no_vpert_in_wake_hard/' 
% end
% 
% if wake_type==3
%     type_folder='no_vpert_in_wake/'
% end
% 
% if wake_type==4
%     type_folder='half_lin_cutoff_half_tot_pert/'
% end
% 
% if wake_type==5
%     type_folder='quarter_lin_cutoff_half_tot_pert/'
% end
% 
% %combination of 3 and 4
% 
% if wake_type==6
%     
%     type_folder='half_lin_cutoff_half_tot_pert_nvpwh/'
%     
% end
% 
% %combination of 3 and 5
% 
% if wake_type==7
%     
%     type_folder='quarter_lin_cutoff_half_tot_pert_nvpwh/'
%     
% end
% 
% if wake_type==8
%     
%     type_folder='half_lin_cutoff_half_tot_pert_nvpw/'
%     
% end
% 
% %combination of 3 and 5
% 
% if wake_type==9
%     
%     type_folder='quarter_lin_cutoff_half_tot_pert_nvpw/'
%     
% end

% mkdir(strcat(root_out))
mkdir(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/'))
mkdir(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/'))

mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/'))
mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/'))

cd('../../../preprocessing')

% test_particle_id=10000;

display(spec)

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_spec( spec);
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,strcat(aux_path,type_folder));

cd ../../parameters/

[ h, OmegaBM ,OmegaCDM ,OmegaM ,OmegaL ,clight ,~ ,t_0 ,Hzero ,tensor_tilt ,spectral_indice ,sigma8 ,T_cmb_t0 ,Scalar_amplitude ] = cosmology(  );

cd ../Analysis/preprocessing/

LinConv=nc/size_box; %spatial linear convertion factor from comoving to simulation
% TimConv=(3*((1+z)^(2))*Hzero*(OmegaM^1/2))/2;  %Convert time in simulation units to seconds

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
    
    TimConv=(3*((1+str2double(redshift_list(rds)))^(2))*Hzero*(OmegaM^(1/2)))/2;  %Convert time in simulation units to seconds
%     TimConv=h*(3*((1+str2double(redshift_list(rds)))^(2))*Hzero*(OmegaM^(1/2)))/2;  %Convert time in simulation units to seconds (the "h" factor is to have v in Mpc/h*s units) 
%     TimConv=(3*((1+str2double(redshift_list(rds)))^(2))*Hzero*(OmegaM^(1/2)))/(h*2);  %Convert time in simulation units to seconds (the "h" factor is to have v in Mpc/h*s units) 
%     TimConv=(3*((1+str2double(redshift_list(rds)))^(3/2))*Hzero*(OmegaM^(1/2)))/2;  %Convert time in simulation units to seconds (with a tentative a dependence fixing) 
%     TimConv=(3*((1+str2double(redshift_list(rds)))^(3/2))*Hzero*(OmegaM^(1/2)))/(h*2);  %Convert time in simulation units to seconds (the "h" factor is to have v in Mpc/h*s units,with a tentative "a" dependence fixing) 

    
    
    
    parfor node=1:length(nodes_list)
        
        filename_out=strcat(root_data_out,spec,aux_path,type_folder,'check/vel/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list(node)),'_Check_Zel_wpid_posZ_vel_z',char(redshift_list(rds)),'.dat');
        fid_o = fopen(filename_out,'w');
        
        filename_out_half=strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list(node)),'_Check_Zel_wpid_half_posZ_vel_z',char(redshift_list(rds)),'.dat');
        fid_o_half = fopen(filename_out_half,'w');
        
        %         display(strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat'))
        filename=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=xv(3,:)+(nc/number_node_dim)*k_node;
        vel_z=xv(6,:)/(LinConv/TimConv);
        
        
        particle_ID_cat=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=fread(fid,[1 Inf],'int64');
        fclose(fid);
        
        [sorted_ sort_inx]=sort(part_id);
        sorted_vel_z=vel_z(sort_inx);
        %         pid_posZ=[part_id(sort_inx);pos_z(sort_inx)];
        
        
        %now with wake
        
        filename=strcat(root,spec,aux_path,type_folder,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=xv(3,:)+(nc/number_node_dim)*k_node;
        vel_z=xv(6,:)/(LinConv/TimConv);
        
        particle_ID_cat=strcat(root,spec,aux_path,type_folder,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=fread(fid,[1 Inf],'int64');
        fclose(fid);
        
        [sorted_w_ sort_inx]=sort(part_id);
        sorted_pos_z_w=pos_z(sort_inx);
        sorted_vel_z_w=vel_z(sort_inx);
        
        %         pid_posZ_w=[part_id(sort_inx);pos_z(sort_inx)];
        
        
        %         find(sorted_w_(:)==sorted_
        
        %         displacement=pid_posZ_w(:,2)-pid_posZ(find(sorted_w_(:)==sorted_),2);
        
        %         lent_out=length(find(sorted_w_(:)==sorted_));
        [comom_ID,com_id_w,com_id_nw]=intersect(sorted_w_,sorted_);
        pos_list_z_w=sorted_pos_z_w(com_id_w);
        part_z_nw=sorted_vel_z(com_id_nw);
        part_z_w=sorted_vel_z_w(com_id_w);
        vel_list=part_z_w-part_z_nw;
        %         displac_list(displac_list>nc/2)=displac_list(displac_list>nc/2)-nc;
        %         displac_list(displac_list<-nc/2)=displac_list(displac_list<-nc/2)+nc;
        
        comom_ID_half=comom_ID(pos_list_z_w>nc/4&pos_list_z_w<3*nc/4);
        pos_list_z_w_half=pos_list_z_w(pos_list_z_w>nc/4&pos_list_z_w<3*nc/4);
        vel_list_half=vel_list(pos_list_z_w>nc/4&pos_list_z_w<3*nc/4);
        
        comom_ID_half=[comom_ID_half 0];
        pos_list_z_w_half=[pos_list_z_w_half nc/2];
        vel_list_half=[vel_list_half 0];
        
        fwrite(fid_o,[comom_ID;pos_list_z_w;vel_list],'float32','l');
        fclose(fid_o);
        
        fwrite(fid_o_half,[comom_ID_half;pos_list_z_w_half;vel_list_half],'float32','l');
        fclose(fid_o_half);
        
    end
    
end





cd('../wake_detection/consist_check/Zeldovich/')


for rds=1:length(redshift_list)
    filename_out=cellstr(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list),'_Check_Zel_wpid_posZ_vel_z',char(redshift_list(rds)),'.dat'));
    vel_diff_ds = fileDatastore(filename_out,'ReadFcn',@read_bin,'FileExtensions','.dat');
    vel_diff=cell2mat(tall(vel_diff_ds));
    
    fig=figure('Visible', 'off');
    histogram2(mod(vel_diff(:,2),nc)*size_box/nc,vel_diff(:,3)*10^17,nc/2,'DisplayStyle','tile','ShowEmptyBins','on');
    mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/'));
    saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_vel_z',char(redshift_list(rds)),'_plot.png'));
    
    fig=figure('Visible', 'off');
    h=histogram(vel_diff(:,3)*10^17);
    mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/hist/'));    
    saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/hist/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_vel_z',char(redshift_list(rds)),'_plot.png'));

    posit_values=vel_diff(vel_diff(:,3)>0,3);
    mn_pos(rds)=gather(mean(posit_values));
    std_pos(rds)=gather(std(posit_values,1));
    
    negat_values=vel_diff(vel_diff(:,3)<0,3);
    mn_neg(rds)=gather(mean(negat_values));
    std_neg(rds)=gather(std(negat_values,1));
    
end

%         filename_out_half=strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list(node)),'_Check_Zel_wpid_half_posZ_vel_z',char(redshift_list(rds)),'.dat');


for rds=1:length(redshift_list)
    filename_out=cellstr(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_node',char(nodes_list),'_Check_Zel_wpid_half_posZ_vel_z',char(redshift_list(rds)),'.dat'));
    vel_diff_ds = fileDatastore(filename_out,'ReadFcn',@read_bin,'FileExtensions','.dat');
    vel_diff=cell2mat(tall(vel_diff_ds));
    
    fig=figure('Visible', 'off');
    h2=histogram2(mod(vel_diff(:,2),nc)*size_box/nc,vel_diff(:,3)*10^17,nc/2,'DisplayStyle','tile','ShowEmptyBins','on');
    mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half'));
    saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_vel_z',char(redshift_list(rds)),'_plot.png'));
    
    fig=figure('Visible', 'off');
    h=histogram(vel_diff(:,3)*10^17);
    mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/hist/'));    
    saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/hist/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_vel_z',char(redshift_list(rds)),'_plot.png'));

    half_posit_values=vel_diff(vel_diff(:,3)>0,3);
    half_mn_pos(rds)=gather(mean(half_posit_values));
    half_std_pos(rds)=gather(std(half_posit_values,1));
    
    half_negat_values=vel_diff(vel_diff(:,3)<0,3);
    half_mn_neg(rds)=gather(mean(half_negat_values));
    half_std_neg(rds)=gather(std(half_negat_values,1));
    

    Pos_Bin_Centers=mean([h2.XBinEdges(1:end-1);h2.XBinEdges(2:end)]);
    Vel_Bin_Centers=mean([h2.YBinEdges(1:end-1);h2.YBinEdges(2:end)]);
    Hist2_val=h2.Values;

    mean_vel_per_slice=(Hist2_val*Vel_Bin_Centers')./sum(Hist2_val')';     
    mean_vel_per_slice(isinf(mean_vel_per_slice)|isnan(mean_vel_per_slice)) = 0;
    fig=figure('Visible', 'off');
    plot(Pos_Bin_Centers,mean_vel_per_slice)
    mkdir(strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/mean_pos/'));    
    saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/mean_pos/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_vel_z',char(redshift_list(rds)),'_plot.png'));
    
    half_mp_mn(rds)=(mean(abs(mean_vel_per_slice)));
    half_mp_std(rds)=(std(abs(mean_vel_per_slice)));    
    
    
end


cd('../../../../parameters')
    
for rds=1:length(redshift_list)    
    [ ~, ~, vel_pert ] = wake( Gmu,str2num(char(redshift_list(rds))));
    wake_vel_pert_zeld(rds,1)=vel_pert;
    
end

cd('../Analysis/wake_detection/consist_check/Zeldovich/')


%plot positive values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,mn_pos*10^17,std_pos*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon: positive')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel_pos','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel_pos.txt'),[(str2num(char(redshift_list))+1).^-1,mn_pos'*10^17,std_pos'*10^17],'delimiter','\t')

%plot negative values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,abs(mn_neg)*10^17,std_neg*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon: negative')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel_neg','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel_neg.txt'),[(str2num(char(redshift_list))+1).^-1,abs(mn_neg')*10^17,std_neg'*10^17],'delimiter','\t')


%plot total values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,((mn_pos+abs(mn_neg))/2)*10^17,((std_pos+std_neg)/2)*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/','_Check_vel_Zel.txt'),[(str2num(char(redshift_list))+1).^-1,((mn_pos'+abs(mn_neg'))/2)*10^17,((std_pos'+std_neg')/2)*10^17],'delimiter','\t')


%same for half plot

%plot total values using mean_position way

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,(half_mp_mn)*10^17,(half_mp_std)*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/mean_pos/','_Check_vel_Zel','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/mean_pos/','_Check_vel_Zel.txt'),[(str2num(char(redshift_list))+1).^-1,((half_mn_pos'+abs(half_mn_neg'))/2)*10^17,((half_std_pos'+half_std_neg')/2)*10^17],'delimiter','\t')


%plot positive values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,half_mn_pos*10^17,half_std_pos*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon: positive')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel_pos','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel_pos.txt'),[(str2num(char(redshift_list))+1).^-1,half_mn_pos'*10^17,half_std_pos'*10^17],'delimiter','\t')

%plot negative values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,abs(half_mn_neg)*10^17,half_std_neg*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon: negative')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel_neg','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel_neg.txt'),[(str2num(char(redshift_list))+1).^-1,abs(half_mn_neg')*10^17,half_std_neg'*10^17],'delimiter','\t')


%plot total values

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,((half_mn_pos+abs(half_mn_neg))/2)*10^17,((half_std_pos+half_std_neg)/2)*10^17)
hold on
plot((str2num(char(redshift_list))+1).^-1,wake_vel_pert_zeld*10^17)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity ((Mpc/h)/s)*10^-17', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon')},'interpreter', 'latex', 'fontsize', 20);
legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northwest')
hold off;

saveas(fig,strcat(root_plot_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel','.png'));
dlmwrite(strcat(root_data_out,spec,aux_path,type_folder,'check/vel/half/','_Check_vel_Zel.txt'),[(str2num(char(redshift_list))+1).^-1,((half_mn_pos'+abs(half_mn_neg'))/2)*10^17,((half_std_pos'+half_std_neg')/2)*10^17],'delimiter','\t')



delete(gcp('nocreate'))

end