function [  ] = proj1d_halos_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim) 
% reads data of the halos 1d projections aconding to the input specifications and plot the result

%   (example) proj1d_halos_plot( '/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','','0.000halo0.dat',1,1,[0,0,0],[0,0],'minmax')


path_in=strcat(root,spec,aux_path);
path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/halos/');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));
mkdir(strcat(root_data_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/halos/'));
path_out=strcat(root_data_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/halos/');

cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
filename_xv=filename;
[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,strcat(filename_xv(1:strfind(filename,'halo')-1),'xv',filename_xv(strfind(filename,'halo')+4:end)));

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/1dproj');

proj1d_halos_data_out( root,root_data_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);

proj1d=dlmread(char(strcat(path_data,'aux/total_dm/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_dm_z',num2str(z),'_data.txt')));

Pos_halos=dlmread(char(strcat(path_data,'aux/Pos_h/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_Pos_halos_z',num2str(z),'_data.txt')));
if (length(Pos_halos)~=1) 
    Pos_halos=Pos_halos*size_box/nc;
else
    Pos_halos=[];
end

Radius=dlmread(char(strcat(path_data,'aux/radius/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_radius_halos_z',num2str(z),'_data.txt')));
if (length(Pos_halos)~=1) 
    Radius=Radius*size_box/nc;
else
    Radius=[];
end

proj1d_h_mass_dc_wrt_dmavr=dlmread(char(strcat(path_data,'h_mass_dc_wrt_dmavr/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_dmdc_z',num2str(z),'_data.txt')));

proj1d_h_mass_dc=dlmread(char(strcat(path_data,'mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_mass_z',num2str(z),'_data.txt')));

proj1d_h_number=dlmread(char(strcat(path_data,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_n_z',num2str(z),'_data.txt')));


%plot the 1d projection of the halo density contrast with respect to the
%total DM content.

fig=figure('Visible', 'off');
set(gcf, 'Position', [0 0 800 400])
ax2 = axes('Position',[0.2 0.2 0.6 0.6]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

xlim ([-inf inf]);

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
if (~ischar(lim))
    plot(ax2,cell_bins1d_z,proj1d_h_mass_dc_wrt_dmavr,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
else 
    plot(ax2,cell_bins1d_z,proj1d_h_mass_dc_wrt_dmavr,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
end

xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax2,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

title(ax2,{strcat('Density contrast of the halo 1d projection'),'with respect to total DM'},'interpreter', 'latex', 'fontsize', 20);

descr = {strcat('z = ',num2str(z));
    strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
    strcat('lenghtFactor = ',num2str(lenght_factor));
    strcat('resolFactor = ',num2str(resol_factor));
    strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
    strcat('box displ wrt centre  = ');
    strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
    strcat('boxDensContr = ');
    num2str((sum(proj1d)-(np)^3)/((np)^3));
    strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc'); 
    strcat('cell/dim = ',num2str(np/lenght_factor));
    strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
    strcat('avrPartInHalo/slice = ');
    strcat(num2str( mean(proj1d)*(mean(proj1d_h_mass_dc_wrt_dmavr)+1)/(resol_factor)));
    strcat('expectedWakeThick = ');
    strcat( num2str(displacement),'\ Mpc');
    strcat('wakeThickResol = ');
    strcat( num2str(displacement/(size_box/(np/(resol_factor)))));
    strcat('peak =',num2str(max(proj1d_h_mass_dc_wrt_dmavr)))};
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');

hold off;

if (~ischar(lim))
    mkdir(path_out,strcat('h_mass_dc_wrt_dmavr/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,'h_mass_dc_wrt_dmavr/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_h_mass_dc_wrt_dmavr_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('h_mass_dc_wrt_dmavr/minmax/'));
    saveas(fig,strcat(path_out,'h_mass_dc_wrt_dmavr/minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_h_mass_dc_wrt_dmavr_z',num2str(z),'_plot.png'));
end




%plot the density contrast 1d projection of the halo mass

fig=figure('Visible', 'off');
set(gcf, 'Position', [0 0 800 400])
ax2 = axes('Position',[0.2 0.2 0.6 0.6]);

average=mean2(proj1d_h_mass_dc);
proj1d_h_mass_dc=(proj1d_h_mass_dc-average)/average;


hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;


xlim ([-inf inf]);


if (~ischar(lim))
    plot(ax2,cell_bins1d_z,proj1d_h_mass_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
else 
    plot(ax2,cell_bins1d_z,proj1d_h_mass_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
end

xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax2,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

title(ax2,{strcat('Density contrast of the'),strcat('1d projection of the halo mass')},'interpreter', 'latex', 'fontsize', 20);



descr = {strcat('z = ',num2str(z));
    strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
    strcat('lenghtFactor = ',num2str(lenght_factor));
    strcat('resolFactor = ',num2str(resol_factor));
    strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
    strcat('box displ wrt centre  = ');
    strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
    strcat('boxDensContr = ');
    num2str((sum(proj1d)-(np)^3)/((np)^3));
    strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc'); 
    strcat('cell/dim = ',num2str(np/lenght_factor));
    strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
    strcat('avrPartInHalo/slice = ');
    strcat(num2str( mean(proj1d)*(mean(proj1d_h_mass_dc_wrt_dmavr)+1)/(resol_factor)));
    strcat('expectedWakeThick = ');
    strcat( num2str(displacement),'\ Mpc');
    strcat('wakeThickResol = ');
    strcat( num2str(displacement/(size_box/(np/(resol_factor)))));
    strcat('peak =',num2str(max(proj1d_h_mass_dc)));
    strcat('$\sigma$ = ',num2str(std(proj1d_h_mass_dc)));
    strcat('$peak/ \sigma$ = ',num2str(max(proj1d_h_mass_dc)/std(proj1d_h_mass_dc)))};
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');

hold off;

if (~ischar(lim))
    mkdir(path_out,strcat('mass_dc/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
    saveas(fig,strcat(path_out,'mass_dc/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_mass_dc_z',num2str(z),'_plot.png'));
else
    mkdir(path_out,strcat('mass_dc/minmax/'));
    saveas(fig,strcat(path_out,'mass_dc/minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_mass_dc_z',num2str(z),'_plot.png'));
end


%plot the 1d projection of the halo number

fig=figure('Visible', 'off');
set(gcf, 'Position', [0 0 800 400])
ax2 = axes('Position',[0.2 0.2 0.6 0.6]);


hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

xlim ([-inf inf]);

proj1d_h_number_dc=(proj1d_h_number-mean(proj1d_h_number))/mean(proj1d_h_number);

plot(ax2,cell_bins1d_z,proj1d_h_number_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);

xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax2,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

title(ax2,{strcat('Density contrast of the'),strcat('1d projection of the halo number')},'interpreter', 'latex', 'fontsize', 20);
descr = {strcat('z = ',num2str(z));
    strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
    strcat('lenghtFactor = ',num2str(lenght_factor));
    strcat('resolFactor = ',num2str(resol_factor));
    strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
    strcat('box displ wrt centre  = ');
    strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
    strcat('boxDensContr = ');
    num2str((sum(proj1d)-(np)^3)/((np)^3));
    strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc'); 
    strcat('cell/dim = ',num2str(np/lenght_factor));
    strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
    strcat('avrHalo/slice = ');
    strcat(num2str( mean(proj1d_h_number)/(resol_factor)));
    strcat('expectedWakeThick = ');
    strcat( num2str(displacement),'\ Mpc');
    strcat('wakeThickResol = ');
    strcat( num2str(displacement/(size_box/(np))));
    strcat('peak =',num2str(max(proj1d_h_number_dc)));
    strcat('$\sigma$ = ',num2str(std(proj1d_h_number_dc)));
    strcat('$peak/ \sigma$ = ',num2str(max(proj1d_h_number_dc)/std(proj1d_h_number_dc)))};
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');
hold off;


mkdir(path_out,strcat('number/'));
saveas(fig,strcat(path_out,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_halos_number_z',num2str(z),'_plot.png'));

colormap default;


end

