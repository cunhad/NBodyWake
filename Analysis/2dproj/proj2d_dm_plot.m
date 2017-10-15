function [  ] = proj2d_dm_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim) 
% reads data of the 2d projections aconding to the input specifications and plot the result

%(example) proj2d_dm_plot('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','','0.000xv0.dat',1,1,[0,0,0],[0,0],'minmax');



path_in=strcat(root,spec,aux_path);

cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,filename);

cd('../2dproj');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));
mkdir(strcat(root_data_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/'));
path_out=strcat(root_data_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');


path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');

proj2d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);
proj2d=dlmread(char(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_data.txt')));

%computes the density contrast

average=mean2(proj2d);
proj2d=(proj2d-average)/average;

%cell_bins1d=[0:nc/(np*resol_factor):nc/lenght_factor];
%cell_bins1d(end)=[];
 
fig=figure('Visible', 'off');

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
if (~ischar(lim))
    clims = lim;
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d,clims);
else 
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d);
end
set(gca,'dataAspectRatio',[1 1 1]);
colorbar;
xlabel('$z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Density contrast of the 2d projection'),strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
hold off;

if (~ischar(lim))
    mkdir(path_out,strcat(num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('minmax/'));
    saveas(fig,strcat(path_out,'minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
end



end

