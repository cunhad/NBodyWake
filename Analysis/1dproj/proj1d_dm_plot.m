function [  ] = proj2d_dm_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim) 
% reads data of the 2d projections aconding to the input specifications and plot the result

%(example) proj1d_dm_plot('/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','','63.000xv0.dat',1,1,[0,0,0],[0,0],'minmax');
%(example) proj1d_dm_plot('/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','','65.000xv0.dat',1,1,[0,0,0],[0,0],[-1 2]);



path_in=strcat(root,spec,aux_path);

cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,filename);

cd('../1dproj');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));
mkdir(strcat(root_data_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/'));
path_out=strcat(root_data_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');


path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');

proj1d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);
proj1d=dlmread(char(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.txt')));

%computes the density contrast

average=mean2(proj1d);
proj1d=(proj1d-average)/average;

%cell_bins1d=[0:nc/(np*resol_factor):nc/lenght_factor];
%cell_bins1d(end)=[];
 
fig=figure('Visible', 'off');

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];

if (~ischar(lim))
    plot(cell_bins1d_z,proj1d,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
   % xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

else 
    plot(cell_bins1d_z,proj1d,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    %xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

end

xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

title({strcat('Density contrast of the 1d projection'),strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
hold off;

if (~ischar(lim))
    mkdir(path_out,strcat(num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('minmax/'));
    saveas(fig,strcat(path_out,'minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
end



end

