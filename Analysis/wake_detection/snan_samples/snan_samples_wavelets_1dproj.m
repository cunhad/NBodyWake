function [   ] = snan_samples_wavelets_1dproj(root,root_snan_in,root_snan_out,spec,aux_path,aux_path_snan_in,aux_path_snan_out,lenght_factor,resol_factor,pivot,rot_angle,cutoff,z_id_range,sample_id_range,info,analysis)

%reads the data from wavelets_proj1d_dm_analysis for the signal to noise analysis
%and creates the corresponding figures

%(example) snan_samples_wavelets_1dproj('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/snan/','/home/asus/Dropbox/extras/storage/graham/small_res/snan_out/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','',1,1,[0,0,0],[0,0],10,'all','all',[1,2,3],[1,2,3]);

% NBody output should be stored as root+spec+aux_path (root directory, specification in the form size_numberofcellsperdimension_number_particlesperdimension_initialredshift_wakespecification&multiplicity, aux_path is the sample number )

% data will be readed in  root_snan_in+spec+aux_path+aux_path_snan_in

% table will be saved in  aux_path_snan_out+spec+aux_path+aux_path_snan_out

% lenght_factor = the analysis cube will have a lateral size given by the
% lateral size of the simulation cube divided by this number

% resol_factor= the bin will hte the particle bin size divided by this
%number

% pivot = a 3d array containing the translation wrt to the center of the
% cube (in grid cell units)

%rot_angle = 2d aray containing the theta and phy spherical angles pointing
%to the direction where the new z axis will be rotated


% z_id_range = array with the redshift id of the requested plots and
% analysis, which starts
% with the highest one equals to 1 and decreasing by unit as the redshift
% is decreased for the id convention. If set to "all" will do for every
% redshift

%sample_id_range: an array containing the id of the samples to be analyzed.
%If set to 'all' every sample will be accounted

% info=[0,1,2,3]
% if info=0, histogram of each plot is generated
% if info=1, minimal plots are generated
% if info=2 complete plots are generated

% analysis=1 -> create a textfile with signal to noise data (peak, std, peak/std)



cd('../../preprocessing');

path_in=strcat(root,spec,aux_path);
[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path);

redshift_list=flip(redshift_list);

 cd('../processing');


if ischar(z_id_range)
    z_id_range=[1 : length(redshift_list)];
end

path_samples_out=strcat(root_snan_out,spec,aux_path_snan_out);
mkdir(strcat(root_snan_out));
mkdir(strcat(root_snan_out),strcat(spec,aux_path_snan_out));

path_samples_in=strcat(root_snan_in,spec,aux_path_snan_in);
sample_list=dir(strcat(path_samples_in,'/sample*'));
sample_list={sample_list.name};
sample_list=sort_nat(sample_list);

if ischar(sample_id_range)
    sample_id_range=[1 : length(sample_list)];
end


 cd('../preprocessing');


snan_data=zeros(length(z_id_range),length(sample_id_range),3);

if ismember(1,analysis)    
    if ismember(1,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/1d_dm_1')));
    end    
    if ismember(2,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/','1d_dm_2')));
    end
    if ismember(3,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/','1d_dm_3')));
    end    
end

if ismember(2,analysis)    
    if ismember(1,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_1')));
    end  
    if ismember(2,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_2')));
    end      
    if ismember(3,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_3')));
    end  
end

if ismember(3,analysis)    
    if ismember(1,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_1')));
    end  
    if ismember(2,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_2')));
    end      
    if ismember(3,info)
        mkdir(strcat(path_samples_out,strcat('/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_3')));
    end  
end

    sample_recount=1;
for sample=sample_id_range
    rds_recount=1;
    for rds = z_id_range
        tot_snan_path_in=strcat(root_snan_in,spec,aux_path_snan_in,'/',char(sample_list(sample)),'/snan/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','1dproj/');
        %the bove is snan_data(rds,sample,:) from 1 to 3 -> peak, average,
        %peak/aver
        snan_data(rds_recount,sample_recount,:)=dlmread(char(strcat(tot_snan_path_in,'dm/','wavelet_filtered_abs_',num2str(cutoff),'MpcCut/','_',num2str(find(str2num(char(redshift_list))==str2num(char(redshift_list(rds))))),'_snan_1dproj_cwt_z',num2str(str2num(char(redshift_list(rds)))),'_data.txt')));
        data1(rds_recount,1)=redshift_list(1,rds);
        data1(rds_recount,2:4)=num2cell(squeeze(snan_data(rds_recount,sample_recount,:)));        
        [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_path( root,spec,aux_path);
        rds_recount=rds_recount+1;
    end
    if ismember(1,analysis)
        if ismember(1,info)
            fig1=figure('Visible', 'off');
%             data1(:,1)=transpose(redshift_list);
%             data1(:,2:4)=num2cell(squeeze(snan_data(:,sample,:)));
            t = uitable(fig1,'Data',data1,'FontSize',8);
            t.ColumnName = {'Redshift','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+30 t.Position(4)+30]);
            saveas(fig1,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/1d_dm_1'),'/snan_eachSample_manyRed_',char(sample_list(sample))),'.png'));
        end
        
        if ismember(2,info)
            fig2=figure('Visible', 'off');
            ax2 = axes('Position',[0 0 1 1],'Visible','off');
%             data1(:,1)=transpose(redshift_list);
%             data1(:,2:4)=num2cell(squeeze(snan_data(:,sample,:)));
            t = uitable(fig2,'Data',data1,'FontSize',8);
            t.ColumnName = {'Redshift','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+50 t.Position(4)+130]);
            tit = {strcat('Signal to noise analysis'),strcat('for $\ $',string(sample_list(sample)))};
            txt=text(0.1,0.8,tit);
            set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
            saveas(fig2,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/1d_dm_2'),'/snan_eachSample_manyRed_',char(sample_list(sample))),'.png'));
            
        end
        if ismember(3,info)
            fig3=figure('Visible', 'off');
            ax2 = axes('Position',[0 0 1 1],'Visible','off');
%             data1(:,1)=transpose(redshift_list);
%             data1(:,2:4)=num2cell(squeeze(snan_data(:,sample,:)));
            t = uitable(fig3,'Data',data1,'FontSize',8);
            t.ColumnName = {'Redshift','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+250 t.Position(4)+130]);
            tit = {strcat('Signal to noise analysis'),strcat('for $\ $',string(sample_list(sample)))};
            txt=text(0.1,0.8,tit);
            set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
            descr = {strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
                strcat('z of wake insertion = ',num2str(ziw));
                strcat('z of simulation init = ',num2str(zi));
                strcat('lenghtFactor = ',num2str(lenght_factor));
                strcat('resolFactor = ',num2str(resol_factor));
                strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
                strcat('box displ wrt centre  = ');
                strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
                strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc');
                strcat('cell/dim = ',num2str(np/lenght_factor));
                strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
                strcat('expecPartic/slice = ');
                strcat(num2str(((np/lenght_factor)^2)/resol_factor))};
            txt2=text(0.65,0.5,descr);
            set(txt2,'Parent',ax2,'interpreter', 'latex','fontsize', 10);            
            saveas(fig3,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_eachSample_manyRed/1d_dm_3'),'/snan_eachSample_manyRed_',char(sample_list(sample))),'.png'));
            
        end        
    end
    
    sample_recount=sample_recount+1;
    
end

if ismember(2,analysis)
    sample_recount=1;
    for sample=sample_id_range
    data2(sample_recount,1)=sample_list(sample);
    sample_recount=sample_recount+1;
    end
    
%     data2(:,1)=transpose(sample_list);
    data2(length(data2)+1,1)={'average'};
    data2(length(data2)+1,1)={'std'};
    rds_recount=1;
    
    for rds = z_id_range        
        data2(1:length(sample_id_range),2)=num2cell(squeeze(snan_data(rds_recount,1:length(sample_id_range),1)));
        data2(1:length(sample_id_range),3)=num2cell(squeeze(snan_data(rds_recount,1:length(sample_id_range),2)));
        data2(1:length(sample_id_range),4)=num2cell(squeeze(snan_data(rds_recount,1:length(sample_id_range),3)));        
        data2(length(sample_id_range)+1,2:4)=num2cell(squeeze(mean(snan_data(rds_recount,1:length(sample_id_range),1:3))));
        data2(length(sample_id_range)+2,2:4)=num2cell(squeeze(std(snan_data(rds_recount,1:length(sample_id_range),1:3))));
        
        if ismember(1,info)                        
            fig1=figure('Visible', 'off');
%             part(:,2:4)=num2cell(squeeze(snan_data(:,sample,:)));
            t = uitable(fig1,'Data',data2,'FontSize',8);
            t.ColumnName = {'','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+30 t.Position(4)+30]);
            saveas(fig1,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_1'),'/snan_manySample_eachRed',num2str(str2num(char(redshift_list(rds))))),'.png'));
        end 
        if ismember(2,info)            
            fig2=figure('Visible', 'off');
            ax2 = axes('Position',[0 0 1 1],'Visible','off');
            t = uitable(fig2,'Data',data2,'FontSize',8);
            t.ColumnName = {'','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+50 t.Position(4)+130]);
            tit = {strcat('Signal to noise analysis'),strcat('for z =$\ $',num2str(str2num(char(redshift_list(rds)))))};
            txt=text(0.1,0.8,tit);
            set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
            saveas(fig2,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_2'),'/snan_manySample_eachRed',num2str(str2num(char(redshift_list(rds))))),'.png'));
        end
        if ismember(3,info)                       
            fig2=figure('Visible', 'off');
            ax2 = axes('Position',[0 0 1 1],'Visible','off');
            t = uitable(fig2,'Data',data2,'FontSize',8);
            t.ColumnName = {'','Peak','Deviation','Peak/Deviation'};
            t.RowName = {};
            t.Position(3:4) = t.Extent(3:4);
            set(gcf, 'Position', [0 0 t.Position(3)+250 t.Position(4)+130]);
            tit = {strcat('Signal to noise analysis'),strcat('for z =$\ $',num2str(str2num(char(redshift_list(rds)))))};
            txt=text(0.1,0.8,tit);
            set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
                descr = {strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
                strcat('z of wake insertion = ',num2str(ziw));
                strcat('z of simulation init = ',num2str(zi));
                strcat('lenghtFactor = ',num2str(lenght_factor));
                strcat('resolFactor = ',num2str(resol_factor));
                strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
                strcat('box displ wrt centre  = ');
                strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
                strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc');
                strcat('cell/dim = ',num2str(np/lenght_factor));
                strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
                strcat('expecPartic/slice = ');
                strcat(num2str(((np/lenght_factor)^2)/resol_factor))};
            txt2=text(0.65,0.5,descr);
            set(txt2,'Parent',ax2,'interpreter', 'latex','fontsize', 10);              
            saveas(fig2,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_eachRed/1d_dm_3'),'/snan_manySample_eachRed',num2str(str2num(char(redshift_list(rds))))),'.png'));
        end        
        rds_recount=rds_recount+1;
    end    
end

if ismember(3,analysis)
    rds_recount=1;
    for rds = z_id_range
        data3(rds_recount,1)=redshift_list(1,rds);
        rds_recount=rds_recount+1;
    end
    data3(1:length(z_id_range),2)=num2cell(squeeze(mean(snan_data(1:length(z_id_range),1:length(sample_id_range),1),2)));
    data3(1:length(z_id_range),3)=num2cell(squeeze(std(snan_data(1:length(z_id_range),1:length(sample_id_range),1),0,2)));
    data3(1:length(z_id_range),4)=num2cell(squeeze(mean(snan_data(1:length(z_id_range),1:length(sample_id_range),2),2)));
    data3(1:length(z_id_range),5)=num2cell(squeeze(std(snan_data(1:length(z_id_range),1:length(sample_id_range),2),0,2)));
    data3(1:length(z_id_range),6)=num2cell(squeeze(mean(snan_data(1:length(z_id_range),1:length(sample_id_range),3),2)));
    data3(1:length(z_id_range),7)=num2cell(squeeze(std(snan_data(1:length(z_id_range),1:length(sample_id_range),3),0,2)));
    
    if ismember(1,info)
        
        fig1=figure('Visible', 'off');
        %             part(:,2:4)=num2cell(squeeze(snan_data(:,sample,:)));
        t = uitable(fig1,'Data',data3,'FontSize',8);
        t.ColumnName = {'redshift','E(peak)','std(peak)','E(Deviation)','std(Deviation)','E(Peak/Deviation)','std(Peak/Deviation)'};
        t.RowName = {};
        t.Position(3:4) = t.Extent(3:4);
        set(gcf, 'Position', [0 0 t.Position(3)+30 t.Position(4)+30]);
        saveas(fig1,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_1'),'/snan_manySample_manyRed'),'.png'));
    end
    if ismember(2,info)
        
        fig1=figure('Visible', 'off');
        ax2 = axes('Position',[0 0 1 1],'Visible','off');
        t = uitable(fig1,'Data',data3,'FontSize',8);
        t.ColumnName = {'redshift','E(peak)','std(peak)','E(Deviation)','std(Deviation)','E(Peak/Deviation)','std(Peak/Deviation)'};
        t.RowName = {};
        t.Position(3:4) = t.Extent(3:4);
        set(gcf, 'Position', [0 0 t.Position(3)+50 t.Position(4)+130]);
        tit = {strcat('Signal to noise analysis'),strcat('summary')};
        txt=text(0.1,0.8,tit);
        set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
        saveas(fig1,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_2'),'/snan_manySample_manyRed'),'.png'));
    end
    if ismember(3,info)
        fig1=figure('Visible', 'off');
        ax2 = axes('Position',[0 0 1 1],'Visible','off');
        t = uitable(fig1,'Data',data3,'FontSize',8);
        t.ColumnName = {'redshift','E(peak)','std(peak)','E(Deviation)','std(Deviation)','E(Peak/Deviation)','std(Peak/Deviation)'};
        t.RowName = {};
        t.Position(3:4) = t.Extent(3:4);
        set(gcf, 'Position', [0 0 t.Position(3)+250 t.Position(4)+130]);
        tit = {strcat('Signal to noise analysis'),strcat('summary')};
        txt=text(0.1,0.8,tit);
        set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
        descr = {strcat('$G\mu = $ ',num2str(Gmu,'%.1E'));
                strcat('z of wake insertion = ',num2str(ziw));
                strcat('z of simulation init = ',num2str(zi));
                strcat('lenghtFactor = ',num2str(lenght_factor));
                strcat('resolFactor = ',num2str(resol_factor));
                strcat('$(\theta,\phi)$ = (',num2str(rot_angle(1)),',',num2str(rot_angle(2)),')' );
                strcat('box displ wrt centre  = ');
                strcat('(',num2str(pivot(1)),',',num2str(pivot(2)),',',num2str(pivot(3)),')',' (cell unit)');
                strcat('boxSize/dim = ',num2str(size_box/lenght_factor),'\ Mpc');
                strcat('cell/dim = ',num2str(np/lenght_factor));
                strcat('sliceSize = ',num2str(size_box/(np/(resol_factor))),'\ Mpc');
                strcat('expecPartic/slice = ');
                strcat(num2str(((np/lenght_factor)^2)/resol_factor));
                strcat('number of samples =',num2str(length(sample_id_range)))};
            txt2=text(0.8,0.5,descr);
            set(txt2,'Parent',ax2,'interpreter', 'latex','fontsize', 10);  
        saveas(fig1,strcat(strcat(path_samples_out,'/',num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/',strcat('wavelet_filtered_abs_',num2str(cutoff),'MpcCut','/snan_manySample_manyRed/1d_dm_3'),'/snan_manySample_manyRed'),'.png'));
    end
end


cd('../wake_detection/snan_samples');


% f = figure;
% part(:,1)=transpose(redshift_list);
% part(:,2:4)=num2cell(squeeze(snan_data(:,2,:)));
% t = uitable(f,'Data',part);
% t.ColumnName = {'Redshift','Peak','Deviation','Peak/Deviation'};
% t.RowName = {};
% t.Position(3:4) = t.Extent(3:4);

% % % fig2=figure('Visible', 'off');
% fig2=figure;
% % % set(gcf, 'Position', [0 0 800 400]);
% ax2 = axes('Position',[0 0 1 1],'Visible','off');
% % 
% % hold on;
% part(:,1)=transpose(redshift_list);
% part(:,2:4)=num2cell(squeeze(snan_data(:,2,:)));
% t = uitable(fig2,'Data',part);
% t.ColumnName = {'Redshift','Peak','Deviation','Peak/Deviation'};
% t.RowName = {};
% t.Position(3:4) = t.Extent(3:4);
% % % title(fig2,{'Density contrast of the 1d projection','of dark matter mass'},'interpreter', 'latex', 'fontsize', 20);
% tit = {strcat('Signal to noise analysis'),strcat('for $\ $ ',string(sample_list(1)))};
% txt=text(0.1,0.9,tit);
% set(txt,'Parent',ax2,'interpreter', 'latex','fontsize', 20);
% % hold off;

    

end

