function [ maxi ] = box_statistics_dm_analysis( root,root_data_out,root_plot_out,spec,aux_path,aux_path_data_out,aux_path_plot_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,lim,data_stream,info,analysis ,filter,cutoff)

%(example)  [ maxi ] = box_statistics_dm_analysis('/home/asus/Dropbox/extras/storage/graham/small_res/', '/home/asus/Dropbox/extras/storage/graham/small_res/data_test/','/home/asus/Dropbox/extras/storage/graham/small_res/test_plot_box/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','','','0.000xv0.dat',2,1,[0,0,0],4,1,4,'minmax',[1,3],[0,1,2,3],1,[0,1],10);

% 
% path_total_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj/');
% path_analysis_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj_analysis/');

cd('../../preprocessing');

[~,redshift_list,~,size_box,nc,np,zi,~,~,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

[  ~,~,~,~,z ] = preprocessing_filename_info( root,spec,aux_path,filename);

cd('../../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu,z);

cd('../Analysis/wake_detection/box_statistics');

if ismember(0,data_stream)    
   [stn_proj1d_angles,stn_filtered_proj1d_angles] =box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,0,filter,cutoff);
else    
    path_data=strcat(strcat(root_data_out,spec,aux_path),'data/',aux_path_data_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/');
    if ismember(2,data_stream)
        if ismember(3,data_stream)
            box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,2,filter,cutoff);
        end
%         proj1d_angles=dlmread(strcat(path_data,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
          stn_proj1d_angles=dlmread(strcat(path_data,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.txt'));
    end
    if ismember(1,data_stream)
        if ismember(3,data_stream)
            box_statistics_dm_data_out( root,root_data_out,spec,aux_path,aux_path_data_out,filename,lenght_factor,resol_factor,pivot,NSIDE,part,num_cores,1,filter,cutoff);
        end
%         bins=[-(nc/(2*lenght_factor)):nc/(np*resol_factor):(nc/(2*lenght_factor))];

        fileID = fopen(strcat(path_data,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
%         proj1d_angles=fread(fileID,[12*NSIDE^2,length(bins)-1],'float32','l');
        %display(strcat(path_data,'box_1dproj/','_',num2str(find(str2num(char(redshift_list))==z)),'_stn_1dproj_angle_z',num2str(z),'_parts',num2str(part),'_NSIDE',num2str(NSIDE),'.bin'));
        stn_proj1d_angles=fread(fileID,[1,12*NSIDE^2],'float32','l');
        fclose(fileID);
    end
end


    
    
    
    %analys=count;
    
    %proj1d_angles=transpose(proj1d_angles);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %computes the peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f=transpose(stn_proj1d_angles);
    
    %Now we can use the s2let package
    
    addpath('/home/asus/Programs/s2let/src/main/matlab','/home/asus/Programs/ssht/src/matlab','/home/asus/Programs/so3/src/matlab','/home/asus');
    
%     s2let_hpx_plot_mollweide(maxi);
    
    sz = size(f);
nsideguessed = sqrt(max(sz)/12);
    L = 2*nsideguessed;
    B=2;
    J_min=0;
    
    [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f,'B',B,'L',L,'J_min',J_min);
    
    % Plot
J = s2let_jmax(L, B);
zoomfactor = 1.2;
ns = ceil(sqrt(2+J-J_min+1)) ;
ny = ns - 1 ;
nx = ns ;
figure('Position',[100 100 1300 1000])

subplot(nx, ny, 1);
s2let_hpx_plot_mollweide(f);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Initial band-limited data')

subplot(nx, ny, 2);
s2let_hpx_plot_mollweide(f_scal);
campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
title('Scaling fct')

for j = J_min:J
   subplot(nx, ny, j-J_min+3);
   s2let_hpx_plot_mollweide(f_wav{j-J_min+1});
   campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
   title(['Wavelet scale : ',int2str(j)-J_min+1])
end 
    
    
%     
%     mkdir(path_analysis_out,strcat('max_1dproj/'));
%     dlmwrite(strcat(path_analysis_out,'max_1dproj/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_max_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi,'delimiter','\t');
% 
% % %     % histogram fo the peaks 
% %     
% % %     maxi=transpose(maxi);
% %     
% %     bins=[0:0.1:10];
% % 
% %     [hist_max edges mid loc] = histcn(maxi,bins);
% %     hist_max=transpose(hist_max);
% % %     
% % %     bins=[0:0.001:10];
% % %     non_zero_ind_of_maxi=find(~hist_max(1,:)==0);
% % %     
% % %      fig=figure('Visible', 'off');
% % %     hold on;
% % %     
% % %     %xlim([0 4]);
% % %     bins(end)=[];
% % %     plot(bins(non_zero_ind_of_maxi(1:end)),hist_max(non_zero_ind_of_maxi(1:end)),'DisplayName',strcat('z = ',char(redshift_list(rds))),'LineWidth',2);
% % %     legend('show');
% % %     xlabel('maximum density', 'interpreter', 'latex', 'fontsize', 20);
% % %     ylabel('number count', 'interpreter', 'latex', 'fontsize', 20);
% % %     set(gca,'FontName','FixedWidth');
% % %     set(gca,'FontSize',16);
% % %     set(gca,'linewidth',2);
% % %     
% % %    title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E')),strcat('skewness',num2str(skewness(hist_max(1,:)))),strcat('kurtosis',num2str(kurtosis(hist_max(1,:))))},'interpreter', 'latex', 'fontsize', 20);
% % % %     title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
% % % 
% % %     mkdir(path_analysis_out,strcat('max_1dproj/histogram_zoom/'));
% % %     path_file_out=strcat(path_analysis_out,'max_1dproj/histogram_zoom/','_',num2str(rds),'_histogram_z',char(redshift_list(rds)),'.png');
% % %     saveas(fig,path_file_out);
% % %     
% % %     hold off;
% % %     
% %          fig=figure('Visible', 'off');
% %     hold on;
% %     
% %     %fix range of the above
% %     %xlim([0 4]);
% % %     bins=[0:0.1:10];
% %     bins(end)=[];
% % %     axis([0 10 0 50000]);
% %     plot(bins(1,:),hist_max(1,:),'LineWidth',2);
% %     legend(strcat('z = ',char(redshift_list(rds))),strcat('skewness',num2str(skewness(hist_max(1,:)))));
% %     xlabel('maximum density', 'interpreter', 'latex', 'fontsize', 20);
% %     ylabel('number count', 'interpreter', 'latex', 'fontsize', 20);
% %     set(gca,'FontName','FixedWidth');
% %     set(gca,'FontSize',16);
% %     set(gca,'linewidth',2);
% %     
% %    title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E')),strcat('skewness',num2str(skewness(hist_max(1,:)))),strcat('kurtosis',num2str(kurtosis(hist_max(1,:))))},'interpreter', 'latex', 'fontsize', 20);
% % %     title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
% % 
% %     mkdir(path_analysis_out,strcat('max_1dproj/histogram/'));
% %     path_file_out=strcat(path_analysis_out,'max_1dproj/histogram/','_',num2str(rds),'_histogram_z',char(redshift_list(rds)),'.png');
% %     saveas(fig,path_file_out);
% %     
% %     hold off;
% 
%     
%       % % signal to noise stat
%   
% % %     %computes the histogram
% % %     
% %   %  cd('../../preprocessing');
% % %     
% %    %  bins=[-1:0.01:10];
% %   %   for i = 1:12*NSIDE^2
% % %     
% %  %    [result(i) edges mid loc] = histcn(count(:,i),bins);
% % %     
% % %     sk(:,i) = skewness(count_1d(:,i));
% % %     kur(:,i)= kurtosis(count_1d(:,i));
% %             
% %         fig=figure('Visible', 'off');
% %         %fig=figure;
% %         
% %      hold on;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %computes the noise and sig to noise
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       st = std(proj1d_angles);
%       st = transpose(st);
% 
%       mkdir(path_analysis_out,strcat('sigma/'));  
%       dlmwrite(strcat(path_analysis_out,'sigma/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_std_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),st,'delimiter','\t');
%       
%       mkdir(path_analysis_out,strcat('sig_to_noise/'));  
%       dlmwrite(strcat(path_analysis_out,'sig_to_noise/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_sig_to_noise_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi(:,1)./st(:,1),'delimiter','\t')
%       
%    
% 
%     
% 
% 
% cd('../wake_detection/box_statistics');

end

