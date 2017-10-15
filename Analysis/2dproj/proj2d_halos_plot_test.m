function [ Pos_halos Radius ] = proj2d_halos_plot_test( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim) 
% reads data of the halos 2d projections aconding to the input specifications and plot the result

%   (example) proj2d_halos_plot( '/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','/home/asus/Dropbox/extras/storage/guillimin/old/','32Mpc_96c_48p_zi63_nowakes','/','','31.000halo0.dat',1,1,[0,0,0],[0,0],'minmax')
%   (example) proj2d_halos_plot( '/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/', '/home/asus/Dropbox/extras/storage/','40Mpc_192c_96p_zi65_nowakes','/','','0.000halo0.dat',1,1,[0,0,0],[0,0],[-1 2])
%   (example) proj2d_halos_plot( '/gs/project/smj-701-aa/disrael/cubep3m/simulations/', '/gs/project/smj-701-aa/disrael/cubep3m/simulations/', '/gs/project/smj-701-aa/disrael/cubep3m/simulations/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','4.000halo0.dat',1,1,[0,0,0],[0,0],[-1 2])


% /sample0001/
% /gs/project/smj-701-aa/disrael/cubep3m/simulations/
% 64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m

% set(gcf,'Renderer','zbuffer');

path_in=strcat(root,spec,aux_path);
path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));
mkdir(strcat(root_data_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/'));
path_out=strcat(root_data_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/');

cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
filename_xv=filename;
[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,strcat(filename_xv(1:strfind(filename,'halo')-1),'xv',filename_xv(strfind(filename,'halo')+4:end)));

cd('../2dproj');

proj2d_halos_data_out( root,root_data_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);

proj2d=dlmread(char(strcat(path_data,'aux/total_dm/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_dm_z',num2str(z),'_data.txt')));

Pos_halos=dlmread(char(strcat(path_data,'aux/Pos_h/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_Pos_halos_z',num2str(z),'_data.txt')));
if (length(Pos_halos)~=1) 
    Pos_halos=Pos_halos*size_box/nc;
else
    Pos_halos=[];
end

Radius=dlmread(char(strcat(path_data,'aux/radius/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_radius_halos_z',num2str(z),'_data.txt')));
if (length(Pos_halos)~=1) 
    Radius=Radius*size_box/nc;
else
    Radius=[];
end

proj2d_h_mass_dc_wrt_dmavr=dlmread(char(strcat(path_data,'h_mass_dc_wrt_dmavr/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_dmdc_z',num2str(z),'_data.txt')));

proj2d_h_mass_dc=dlmread(char(strcat(path_data,'mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_z',num2str(z),'_data.txt')));

proj2d_h_number=dlmread(char(strcat(path_data,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_n_z',num2str(z),'_data.txt')));
%plot the 2d projection of the halo on top of DM density contrast


%computes the density contrast

average=mean2(proj2d);
proj2d=(proj2d-average)/average;

%cell_bins1d=[0:nc/(np*resol_factor):nc/lenght_factor];
%cell_bins1d(end)=[];
 
%fig=figure('Visible', 'off');
fig=figure('Visible', 'on');

set(gcf,'Renderer','painters');


set(fig, 'DefaultFigureRenderer', 'painters');
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
title({strcat('Density contrast of the 2d projection'),'plus halos',strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);

%plot halos
if(~isempty(Pos_halos))
    scatter(Pos_halos(3,:),Pos_halos(2,1:100),10*Radius(1,1:100),'r','filled');
end



hold off;

if (~ischar(lim))
    mkdir(path_out,strcat('dmdc_plus_haloscatter/',strcat(num2str(lim(1)),'_',num2str(lim(2))),'lim','/'));
    saveas(fig,strcat(path_out,'dmdc_plus_haloscatter/',strcat(num2str(lim(1)),'_',num2str(lim(2))),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_dmdc_plus_haloscatter_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('dmdc_plus_haloscatter/minmax/'));
    saveas(fig,strcat(path_out,'dmdc_plus_haloscatter/minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_dmdc_plus_haloscatter_z',num2str(z),'_plot.png'));
end


%image(fig);


%plot the 2d projection of the halo density contrast with respect to the
%total DM content.

fig=figure('Visible', 'off');

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
if (~ischar(lim))
    clims = lim;
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_h_mass_dc_wrt_dmavr,clims);
else 
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_h_mass_dc_wrt_dmavr);
end
set(gca,'dataAspectRatio',[1 1 1]);
colorbar;
xlabel('$z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Density contrast of the 2d halo projection'),strcat('with respect to total DM'),strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
hold off;

if (~ischar(lim))
    mkdir(path_out,strcat('h_mass_dc_wrt_dmavr/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,'h_mass_dc_wrt_dmavr/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_h_mass_dc_wrt_dmavr_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('h_mass_dc_wrt_dmavr/minmax/'));
    saveas(fig,strcat(path_out,'h_mass_dc_wrt_dmavr/minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_h_mass_dc_wrt_dmavr_z',num2str(z),'_plot.png'));
end

% image(fig);


%plot the density contrast 2d projection of the halo mass

fig=figure('Visible', 'off');

average=mean2(proj2d_h_mass_dc);
proj2d_h_mass_dc=(proj2d_h_mass_dc-average)/average;


hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
if (~ischar(lim))
    clims = lim;
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_h_mass_dc,clims);
else 
    imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_h_mass_dc);
end
set(gca,'dataAspectRatio',[1 1 1]);
colorbar;
xlabel('$z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Density contrast of the'),strcat('2d halo mass projection'),strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
hold off;

if (~ischar(lim))
    mkdir(path_out,strcat('mass_dc/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,'mass_dc/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_dc_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('mass_dc/minmax/'));
    saveas(fig,strcat(path_out,'mass_dc/minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_dc_z',num2str(z),'_plot.png'));
end


%plot the 2d projection of the halo number

fig=figure('Visible', 'off');

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])

imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_h_number);
set(gca,'dataAspectRatio',[1 1 1]);
colorbar;
%colormap(hot);
xlabel('$z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('$y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Halo number'),strcat('of 2d count projection'),strcat('at z =',num2str(z),' for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
hold off;


mkdir(path_out,strcat('number/'));
saveas(fig,strcat(path_out,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_number_z',num2str(z),'_plot.png'));

colormap default;


end

