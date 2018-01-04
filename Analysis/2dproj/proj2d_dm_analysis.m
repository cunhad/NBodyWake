function [  ] = proj2d_dm_analysis( root,root_data_out,root_plot_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim,data_stream,info)
% reads data of the 2d projections aconding to the input specifications and plot the result

%(example) proj2d_dm_analysis('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data/','/home/asus/Dropbox/extras/storage/graham/small_res/plot/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','','0.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[1,3],[0,1,2,3]);
%(example) proj2d_dm_analysis('/home/asus/Dropbox/extras/storage/guillimin/','/home/asus/Dropbox/extras/storage/guillimin/data/','/home/asus/Dropbox/extras/storage/guillimin/plot/','64Mpc_1024c_512p_zi63_wakeGmu1t10m7zi31m','/sample0001/','','','15.000xv0.dat',1,1,[0,0,0],[0,0],'minmax',[1,3],[0,1,2,3]);

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


path_in=strcat(root,spec,aux_path);

cd('../preprocessing');

[~,redshift_list,~,size_box,nc,np,zi,~,~,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/2dproj');

mkdir(root_plot_out);
mkdir(root_plot_out,strcat(spec,aux_path));
mkdir(strcat(root_plot_out,spec,aux_path),strcat('plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/'));
tot_plot_path_out=strcat(root_plot_out,spec,aux_path,'plot/',aux_path_plot_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/');



if ismember(0,data_stream)    
   [ ~,~,proj2d] = proj2d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,0);
else    
    path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/dm/');
    if ismember(2,data_stream)
        if ismember(3,data_stream)
            proj2d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,2);
        end
        proj2d=dlmread(char(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_data.txt')));
    end
    if ismember(1,data_stream)
        if ismember(3,data_stream)
            proj2d_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,rot_angle,1);
        end
        cell_bins1d_y=[(nc/2)-(nc/(2*lenght_factor))+pivot(2):nc/(np*resol_factor):(nc/2)+(nc/(2*lenght_factor))+pivot(2)];

        fileID = fopen(strcat(path_data,'dc/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_data.bin'));
        proj2d=fread(fileID,[length(cell_bins1d_y)-1,length(cell_bins1d_y)-1],'float32','l');
        fclose(fileID);
    end
end

%computes the density contrast

average=mean2(proj2d);
proj2d_dc=(proj2d-average)/average;




gcc_to_mpc=size_box/nc;
pvy=pivot(2)*gcc_to_mpc;
pvz=pivot(3)*gcc_to_mpc;

if ismember(0,info)
    
    fig0=figure('Visible', 'off');
    set(gcf, 'Position', [0 0 800 600]);
    ax2 = axes('Position',[0.15 0.13 0.5 0.7]);
    
    if (~ischar(lim))
    his=histogram( ax2,proj2d_dc,'BinLimits',lim );
    else
    his=histogram( ax2,proj2d_dc);        
    end
    
    hold on;
    
    %axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
    
    %set(gca,'dataAspectRatio',[1 1 1]);
    %colorbar;
    xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel(ax2,'$Y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title(ax2,{strcat('Histogram of the density contrast'),'of the 2d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);
    descr = {strcat('z = ',num2str(z));
        strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
        strcat('lenghtFactor = ',num2str(lenght_factor));
        strcat('resolFactor = ',num2str(resol_factor));
        strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
        strcat('box displ wrt centre  = ');
        strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
        strcat('boxDensContr = ');
        num2str((sum(sum(proj2d))-(np)^3)/((np)^3));
        strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc');
        strcat('cell/dim = ',num2str(np/lenght_factor));
        strcat('resolution = ',num2str(size_box/(np)),'\ Mpc');
        strcat('expectedWakeThick = ');
        strcat( num2str(displacement),'\ Mpc');
        strcat('wakeThickResol = ');
        strcat( num2str(displacement/(size_box/(np/(resol_factor)))));
        strcat('$\sigma$ = ',num2str(std(proj2d_dc(:))));
        strcat('skewness = ',num2str(skewness(proj2d_dc(:))));
        strcat('kurtosis = ',num2str(kurtosis(proj2d_dc(:))));
        strcat('num of bins = ',num2str(his.NumBins))};
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    txt=text(0.75,0.5,descr);
    set(txt,'Parent',ax1,'interpreter', 'latex');
    
    
end

if ismember(1,info)
    
    fig1=figure('Visible', 'off');
    set(gcf, 'Position', [0 0 600 600]);
    
    
    hold on;
    axes('Position',[0 0 1 1],'Visible','off');
    
    axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
    if (~ischar(lim))
        clims = lim;
        imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc,clims);
    else
        imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc);
    end
    
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    
    if (~ischar(lim))
        fig_cb1=figure('Visible', 'off');
        caxis([lim(1) lim(2)]);
        colorbar;
        axis off
        left=800; bottom=100 ; width=60 ; height=600;
        pos=[left bottom width height];
        set(fig_cb1,'OuterPosition',pos);
        
        fig_cb2=figure('Visible', 'off');
        caxis([lim(1) lim(2)]);
        colorbar('location','Southoutside');
        axis off
        left=800; bottom=110 ; width=600 ; height=130;
        pos=[left bottom width height];
        set(fig_cb2,'OuterPosition',pos);
    else
        fig_cb1=figure('Visible', 'off');
        caxis([min(min(proj2d_dc)) max(max(proj2d_dc))]);
        colorbar;
        axis off
        left=800; bottom=100 ; width=60 ; height=600;
        pos=[left bottom width height];
        set(fig_cb1,'OuterPosition',pos);
        
        fig_cb2=figure('Visible', 'off');
        caxis([min(min(proj2d_dc)) max(max(proj2d_dc))]);
        colorbar('location','Southoutside');
        axis off
        left=800; bottom=110 ; width=600 ; height=130;
        pos=[left bottom width height];
        set(fig_cb2,'OuterPosition',pos);
        
    end
    
end

if ismember(2,info)
    
    fig2=figure('Visible', 'off');
    
    hold on;
    
    axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
    if (~ischar(lim))
        clims = lim;
        imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc,clims);
    else
        imagesc([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc);
    end
    set(gca,'dataAspectRatio',[1 1 1]);
    colorbar;
    xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel('$Y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title({strcat('Density contrast of the 2d projection'),'of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);
    
end

if ismember(3,info)
    
    fig3=figure('Visible', 'off');
    set(gcf, 'Position', [0 0 800 600]);
    ax2 = axes('Position',[0.05 0.13 0.7 0.7]);
    
    peak2d=max(max(proj2d_dc));
    [row_peak2d,col_peak2d] = find(proj2d_dc==peak2d);
    row_peak2d= -size_box/(2*lenght_factor)+size_box/(2)+pvy+(size_box/lenght_factor)*(row_peak2d/(np*resol_factor/lenght_factor));
    col_peak2d= -size_box/(2*lenght_factor)+size_box/(2)+pvz+(size_box/lenght_factor)*(col_peak2d/(np*resol_factor/lenght_factor));
    
    hold on;
    
    axis([-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz])
    if (~ischar(lim))
        clims = lim;
        imagesc(ax2,[-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc,clims);
    else
        imagesc(ax2,[-size_box/(2*lenght_factor)+size_box/(2)+pvy size_box/(2*lenght_factor)+size_box/(2)+pvy],[ -size_box/(2*lenght_factor)+size_box/(2)+pvz size_box/(2*lenght_factor)+size_box/(2)+pvz],proj2d_dc);
    end
    set(gca,'dataAspectRatio',[1 1 1]);
    colorbar;
    xlabel(ax2,'$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    ylabel(ax2,'$Y(Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
    set(gca,'FontName','FixedWidth');
    set(gca,'FontSize',16);
    set(gca,'linewidth',2);
    title(ax2,{strcat('Density contrast of the 2d projection'),'of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);
    descr = {strcat('z = ',num2str(z));
        strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
        strcat('lenghtFactor = ',num2str(lenght_factor));
        strcat('resolFactor = ',num2str(resol_factor));
        strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
        strcat('box displ wrt centre  = ');
        strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
        strcat('boxDensContr = ');
        num2str((sum(sum(proj2d))-(np)^3)/((np)^3));
        strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc');
        strcat('cell/dim = ',num2str(np/lenght_factor));
        strcat('resolution = ',num2str(size_box/(np)),'\ Mpc');
        strcat('expectedWakeThick = ');
        strcat( num2str(displacement),'\ Mpc');
        strcat('wakeThickResol = ');
        strcat( num2str(displacement/(size_box/(np/(resol_factor)))));
        strcat('peak = ',num2str(peak2d));
        strcat('peak location = ');
        strcat('(',num2str(col_peak2d),',',num2str(row_peak2d),')')};
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    txt=text(0.75,0.5,descr);
    set(txt,'Parent',ax1,'interpreter', 'latex');
    
    
end

if ismember(0,info)
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig0,strcat(tot_plot_path_out,'dm_0/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_hist_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_0/','minmax/'));
        saveas(fig0,strcat(tot_plot_path_out,'dm_0/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_hist_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(1,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig1,strcat(tot_plot_path_out,'dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
        saveas(fig_cb1,strcat(tot_plot_path_out,'dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','colormap_v_.png'));
        saveas(fig_cb2,strcat(tot_plot_path_out,'dm_1/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','colormap_h.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_1/','minmax/'));
        saveas(fig1,strcat(tot_plot_path_out,'dm_1/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
        saveas(fig_cb1,strcat(tot_plot_path_out,'dm_1/','minmax/','colormap_v_.png'));
        saveas(fig_cb2,strcat(tot_plot_path_out,'dm_1/','minmax/','colormap_h.png'));
    end
    
end

if ismember(2,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig2,strcat(tot_plot_path_out,'dm_2/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_2/','minmax/'));
        saveas(fig2,strcat(tot_plot_path_out,'dm_2/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
    end
    
end

if ismember(3,info)
    
    if (~ischar(lim))
        mkdir(tot_plot_path_out,strcat('dm_3/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/'));
        saveas(fig3,strcat(tot_plot_path_out,'dm_3/',num2str(lim(1)),'_',num2str(lim(2)),'lim','/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
    else
        mkdir(tot_plot_path_out,strcat('dm_3/','minmax/'));
        saveas(fig3,strcat(tot_plot_path_out,'dm_3/','minmax/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_z',num2str(z),'_plot.png'));
    end
    
end



end

