function [ test_particle_displacement ] = displacement_evol( root,root_out,spec,aux_path )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% (example) []=displacement_evol( '/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/check/','64Mpc_96c_48p_zi255_wakeGmu5t10m6zi63m','/sample1001/')


cd('../../../preprocessing')

test_particle_id=10000;

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_spec( spec);

% [ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path);
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

wake_displacement=zeros(length(redshift_list),2);
wake_displacement_zeld=zeros(length(redshift_list),1);
test_particle_pos=zeros(length(redshift_list),1);
test_particle_wpos=zeros(length(redshift_list),1);
test_particle_displacement=zeros(length(redshift_list),1);



for rds=1:length(redshift_list)
    % for rds=5:5
    
    pos_z=[];
    part_id=[];
    
    number_node_dim=nthroot(numel(nodes_list), 3);
    
    for node=1:length(nodes_list)
        
        filename=strcat(root,spec,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=[pos_z xv(3,:)+(nc/number_node_dim)*k_node];
        
        particle_ID_cat=strcat(root,spec,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=[part_id fread(fid,[1 Inf],'int64')];
        fclose(fid);
    end
    
    id_p_z=[part_id;pos_z];
    w_sort_id_p_z=sort(transpose(id_p_z));
    
    test_particle_wpos(rds,1)=w_sort_id_p_z(test_particle_id,2);
    
    pos_z=[];
    part_id=[];
    
    for node=1:length(nodes_list)
        
        spec_nowake=strcat(string(size_box),'Mpc_',string(nc),'c_',string(np),'p_zi',string(zi),'_nowakem');
        
        filename=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat');
        fid = fopen(filename);
        fread(fid, [12 1], 'float32','l') ;
        xv=fread(fid, [6 Inf], 'float32','l');
        fclose(fid);
        
        node_ID=node-1;
        k_node=floor(node_ID/number_node_dim^2);
        
        pos_z=[pos_z xv(3,:)+(nc/number_node_dim)*k_node];
        
        particle_ID_cat=strcat(root,spec_nowake,aux_path,char(redshift_list(rds)),'PID',char(nodes_list(node)),'.dat');
        fid = fopen(particle_ID_cat);
        fread(fid,6,'int64');
        part_id=[part_id fread(fid,[1 Inf],'int64')];
        fclose(fid);
    end
    
    id_p_z=[part_id;pos_z];
    nw_sort_id_p_z=sort(transpose(id_p_z));
    
    test_particle_pos(rds,1)=nw_sort_id_p_z(test_particle_id,2);
    
    displ=w_sort_id_p_z(:,2)-nw_sort_id_p_z(:,2);
    avr=mean(displ);
    displ(:)=displ(:)-avr;
    
    % fig=figure;
    fig=figure('Visible', 'off');
    hold on;
    plot(nw_sort_id_p_z(:,2)*size_box/nc,displ*size_box/nc);
    xlim ([-inf inf]);
    xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('Displacement (Mpc/h)', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title({strcat('Displacement for z = '),char(redshift_list(rds))},'interpreter', 'latex', 'fontsize', 20);
    hold off;
    
    
    mkdir(strcat(root_out,spec,aux_path,'check/'));
    saveas(fig,strcat(root_out,spec,aux_path,'check/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_Check_Zel_z',char(redshift_list(rds)),'.png'));
    
    
    positive_displ=displ(displ>=0);
    negative_displ=displ(displ<0);
    
    % display(mean(positive_displ));
    % display(std(positive_displ));
    % display(mean(negative_displ));
    % display(std(negative_displ));
    
    wake_displacement(rds,1)=(abs(mean(positive_displ))+abs(mean(negative_displ)))/2;
    wake_displacement(rds,2)=(std(positive_displ)+std(negative_displ))/2;
    
    cd('../../parameters')
    
    [ ~, displacement, ~ ] = wake( Gmu,str2num(char(redshift_list(rds))));
    wake_displacement_zeld(rds,1)=displacement;
    
    cd('../Analysis/preprocessing/')
    
end

% test_particle_displacement=test_particle_pos-test_particle_wpos;

cd('../wake_detection/consist_check/Zeldovich')

fig=figure('Visible', 'off');
% fig=figure;
errorbar((str2num(char(redshift_list))+1).^-1,wake_displacement(:,1)*size_box/nc,wake_displacement(:,2)*size_box/nc)
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
    
    saveas(fig,strcat(root_out,spec,aux_path,'check/','_Check_Zel','.png'));
figure;
plot((str2num(char(redshift_list))+1).^-1,test_particle_displacement);

end

