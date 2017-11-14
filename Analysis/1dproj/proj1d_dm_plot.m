function [  ] = proj1d_dm_plot( root,root_data_out,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim,info) 
% reads data of the 2d projections aconding to the input specifications and plot the result

%(example) proj1d_dm_plot('/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','','0.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[0,1,2]);
%(example) proj1d_dm_plot('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','15.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[0,1,2]);

path_in=strcat(root,spec,aux_path);

cd('../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header i_node j_node k_node number_node_dim ] = preprocessing_nodes_all_but_phasespace( root,spec,aux_path,filename);

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/1dproj');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));
mkdir(strcat(root_data_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/'));
path_out=strcat(root_data_out,spec,aux_path,'plot/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/');


path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');

proj1d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle);
proj1d=dlmread(char(strcat(path_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.txt')));

%computes the density contrast

average=mean2(proj1d);
proj1d_dc=(proj1d-average)/average;

%cell_bins1d=[0:nc/(np*resol_factor):nc/lenght_factor];
%cell_bins1d(end)=[];



if ismember(0,info)


fig0=figure('Visible', 'off');
set(gcf, 'Position', [0 0 800 400]);

ax2 = axes('Position',[0.2 0.2 0.6 0.6]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
xlim ([-inf inf]);

    if (~ischar(lim))
    his=histogram( ax2,proj1d_dc,'BinLimits',lim );
    else
    his=histogram( ax2,proj1d_dc);        
    end


xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax2,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

st = std(proj1d_dc);

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title(ax2,{'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

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
    strcat('expectedWakeThick = ');
    strcat( num2str(displacement),'\ Mpc');
    strcat('wakeThickResol = ');
    strcat( num2str(displacement/(size_box/(np))));
    strcat('expecPartic/slice = ');
    strcat(num2str(((np/lenght_factor)^2)/resol_factor));
    strcat('peak =',num2str(max(proj1d_dc)));
    strcat('$\sigma$ = ',num2str(st));
    strcat('$peak/ \sigma$ = ',num2str(max(proj1d_dc)/st));
        strcat('skewness = ',num2str(skewness(proj1d_dc(:))));
        strcat('kurtosis = ',num2str(kurtosis(proj1d_dc(:))));
        strcat('num of bins = ',num2str(his.NumBins))};


ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');


hold off;

end

if ismember(1,info)


fig1=figure('Visible', 'off');
%figure('position', [0 0 1 2]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 600 400]);
set(gcf, 'Position', [0 0 600 400]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
xlim ([-inf inf]);

if (~ischar(lim))
    plot(cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
   % xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

else 
    plot(cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    %xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

end


xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

st = std(proj1d_dc);

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title({'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

hold off;

end


if ismember(2,info)


fig2=figure('Visible', 'off');
%figure('position', [0 0 1 2]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 600 400]);
set(gcf, 'Position', [0 0 800 400]);

ax2 = axes('Position',[0.2 0.2 0.6 0.6]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
xlim ([-inf inf]);

if (~ischar(lim))
    plot(ax2,cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
   % xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

else 
    plot(ax2,cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    %xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

end


xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax2,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

st = std(proj1d_dc);

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title(ax2,{'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

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
    strcat('expectedWakeThick = ');
    strcat( num2str(displacement),'\ Mpc');
    strcat('wakeThickResol = ');
    strcat( num2str(displacement/(size_box/(np))));
    strcat('expecPartic/slice = ');
    strcat(num2str(((np/lenght_factor)^2)/resol_factor));
    strcat('peak =',num2str(max(proj1d_dc)));
    strcat('$\sigma$ = ',num2str(st));
    strcat('$peak/ \sigma$ = ',num2str(max(proj1d_dc)/st))};
%axes(ax1); % sets ax1 to current axes
%fig.CurrentAxes = ax1;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');

hold off;

end

if ismember(0,info)
    
    if (~ischar(lim))
        mkdir(path_out,strcat('dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig0,strcat(path_out,'dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(path_out,strcat('dm_0/','minmax/'));
        saveas(fig0,strcat(path_out,'dm_0/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(1,info)
    
    if (~ischar(lim))
        mkdir(path_out,strcat('dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig1,strcat(path_out,'dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(path_out,strcat('dm_1/','minmax/'));
        saveas(fig1,strcat(path_out,'dm_1/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(2,info)
    
    if (~ischar(lim))
        mkdir(path_out,strcat('dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig2,strcat(path_out,'dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(path_out,strcat('dm_2/','minmax/'));
        saveas(fig2,strcat(path_out,'dm_2/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end



end

