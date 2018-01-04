function [  ] = proj1d_dm_analysis( root,root_data_out,root_plot_out,root_snan_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,aux_path_snan_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim,data_stream,info,analysis) 
% reads (and/or generate) data of the 1d projections aconding to the input specifications and plot the result

%(example) proj1d_dm_analysis('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','/home/asus/Dropbox/extras/storage/graham/small_res/snan/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','','','0.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[1,3],[0,1,2,3],1);
%(example) proj1d_dm_analysis('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/data/','/home/asus/Dropbox/extras/storage/guillimin/plot/','/home/asus/Dropbox/extras/storage/guillimin/snan/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','','','15.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[1,3],[0,1,2,3],1);

% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% plot will be stored in  root_plot_out+spec+aux_path+aux_path_plot_out

% if specified, data will be stored in  root_data_out+spec+aux_path+aux_path_data_out

% if specified, signal to noise analysis will be stored in  root_snan_out+spec+aux_path+aux_path_snan_out

% filename is the output file from the nbody simulation

% lenght_factor = the analysis cube will have a lateral size given by the
% lateral size of the simulation cube divided by this number

% resol_factor= the bin will hte the particle bin size divided by this
%number

% pivot = a 3d array containing the translation wrt to the center of the
% cube (in grid cell units)

%rot_angle = 2d aray containing the theta and phy spherical angles pointing
%to the direction where the new z axis will be rotated

%lim= limits on the y axis of the plot, in array format. If set to 'minmax'
%will display between the min and max values

% data_stream=[1,2,3]
% if data_stream = 0, no data output generated and readed, the data is
% passed directily to this program
% if data_stream = 1, reads data binaries 
% if data_stream = 2, reads data text 
% if data_stream = 3, generates the data output in binary or text if 1 or 2 options are given, respectively


% info=[0,1,2,3]
% if info=0, histogram of each plot is generated
% if info=1, minimal plots are generated
% if info=2 complete plots are generated
% if info=3 complete plots plus info text are generated


% analysis=1 -> create a textfile with signal to noise data (peak, std, peak/std)


path_in=strcat(root,spec,aux_path);

cd('../preprocessing');

[~,redshift_list,~,size_box,nc,np,zi,~,~,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/1dproj');

mkdir(root_plot_out);
mkdir(root_plot_out,strcat(spec,aux_path));
mkdir(strcat(root_plot_out,spec,aux_path),strcat('plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/'));
tot_plot_path_out=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/');

if ismember(0,data_stream)    
   [proj1d] = proj1d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,0);
else    
    path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/dm/');   
    if ismember(2,data_stream)
        if ismember(3,data_stream)
            proj1d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,2);
        end
        proj1d=dlmread(char(strcat(path_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.txt')));    
    end
    if ismember(1,data_stream)
        if ismember(3,data_stream)
            proj1d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,1);
        end
        fileID = fopen(strcat(path_data,'nc/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_data.bin'));
        proj1d=fread(fileID,'float32','l');
        fclose(fileID);
    end
end
%computes the density contrast

average=mean2(proj1d);
proj1d_dc=(proj1d-average)/average;

st_proj1d_dc = std(proj1d_dc);
max_proj1d_dc=max(proj1d_dc);

%cell_bins1d=[0:nc/(np*resol_factor):nc/lenght_factor];
%cell_bins1d(end)=[];



if ismember(0,info)


fig0=figure('Visible', 'off');
set(gcf, 'Position', [0 0 800 450]);

ax3 = axes('Position',[0.15 0.2 0.6 0.6]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
xlim ([-inf inf]);

    if (~ischar(lim))
    his=histogram( ax3,proj1d_dc,'BinLimits',lim );
    else
    his=histogram( ax3,proj1d_dc);        
    end


xlabel(ax3,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax3,'frequency', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title(ax3,{'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

descr = {strcat('z = ',num2str(z));
    strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
    strcat('z of wake insertion = ',num2str(ziw));
    strcat('z of simulation init = ',num2str(zi));
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
    strcat('peak =',num2str(max_proj1d_dc));
    strcat('$\sigma$ = ',num2str(st_proj1d_dc));
    strcat('$peak/ \sigma$ = ',num2str(max_proj1d_dc/st_proj1d_dc));
        strcat('skewness = ',num2str(skewness(proj1d_dc(:))));
        strcat('kurtosis = ',num2str(kurtosis(proj1d_dc(:))));
        strcat('num of bins = ',num2str(his.NumBins))};


ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.8,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');


hold off;

end

if ismember(1,info)


fig1=figure('Visible', 'off');
%figure('position', [0 0 1 2]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 600 400]);
set(gcf, 'Position', [0 0 800 400]);

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


% xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

% title({'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

hold off;

end

if ismember(2,info)


fig2=figure('Visible', 'off');
%figure('position', [0 0 1 2]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 600 400]);
set(gcf, 'Position', [0 0 800 400]);

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

%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title({'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

hold off;

end


if ismember(3,info)


fig3=figure('Visible', 'off');
%figure('position', [0 0 1 2]);
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperPosition', [0 0 600 400]);
set(gcf, 'Position', [0 0 800 400]);

ax3 = axes('Position',[0.2 0.2 0.6 0.6]);

hold on;

gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

cell_bins1d_z=[(size_box/2)-(size_box/(2*lenght_factor))+pvz:size_box/(np*resol_factor):(size_box/2)+(size_box/(2*lenght_factor))+pvz];
cell_bins1d_z(end)=[];
xlim ([-inf inf]);

if (~ischar(lim))
    plot(ax3,cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    ylim(lim);
   % xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

else 
    plot(ax3,cell_bins1d_z,proj1d_dc,'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
    %xlim([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy]);

end


xlabel(ax3,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
ylabel(ax3,'Density contrast', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);


%legend('Location','eastoutside');
%legend([strcat('peak =',num2str(maxi)),'//',strcat('sigma =',num2str(st)),strcat('peak/sigma =',num2str(maxi/st))],'interpreter', 'latex');

title(ax3,{'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);

descr = {strcat('z = ',num2str(z));
    strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
    strcat('z of wake insertion = ',num2str(ziw));
    strcat('z of simulation init = ',num2str(zi));
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
    strcat('peak =',num2str(max_proj1d_dc));
    strcat('$\sigma$ = ',num2str(st_proj1d_dc));
    strcat('$peak/ \sigma$ = ',num2str(max_proj1d_dc/st_proj1d_dc))};
%axes(ax1); % sets ax1 to current axes
%fig.CurrentAxes = ax1;
ax1 = axes('Position',[0 0 1 1],'Visible','off');
txt=text(0.82,0.5,descr);
set(txt,'Parent',ax1,'interpreter', 'latex');

hold off;

end

if ismember(0,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig0,strcat(tot_plot_path_out,'dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_0/','minmax/'));
        saveas(fig0,strcat(tot_plot_path_out,'dm_0/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(1,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig1,strcat(tot_plot_path_out,'dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_1/','minmax/'));
        saveas(fig1,strcat(tot_plot_path_out,'dm_1/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(2,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig2,strcat(tot_plot_path_out,'dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_2/','minmax/'));
        saveas(fig2,strcat(tot_plot_path_out,'dm_2/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(3,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_3/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig3,strcat(tot_plot_path_out,'dm_3/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_3/','minmax/'));
        saveas(fig3,strcat(tot_plot_path_out,'dm_3/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(1,analysis)
    mkdir(root_snan_out);
    mkdir(root_snan_out,strcat(spec,aux_path));
    mkdir(strcat(root_snan_out,spec,aux_path),strcat('snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/'));
    tot_snan_path_out=strcat(root_snan_out,spec,aux_path,'snan/',aux_path_snan_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/');
    snan=[max_proj1d_dc st_proj1d_dc (max_proj1d_dc)/(st_proj1d_dc)];  
    mkdir(tot_snan_path_out,strcat('dm/'));
    dlmwrite(strcat(tot_snan_path_out,'dm/','_',num2str(find(str2num(char(redshift_list))==z)),'_snan_1dproj_z',num2str(z),'_data.txt'),snan,'delimiter','\t');
end



end
