function [  ] = box_statistics_dm_from_tot_to_analysis_wwavel( root_per_node_out,spec,aux_path,aux_path_per_node_out,NSIDE,cutoff )

%(example) box_statistics_dm_from_tot_to_analysis_wwavel('/home/asus/Dropbox/extras/storage/guillimin/test/','64Mpc_96c_48p_zi63_nowakes','/','',4,10);


pivot=[0,0,0]; %this is the position od the origin of the rotation point with respect to the center of the box
lenght_factor=2;
resol_factor=2;
    
%reads the specifications and extract the information on variables
spec_arr = strsplit(spec,'_');

%extract the box size

size_box = spec_arr(1);
size_box = char(size_box);
size_box = size_box(1:end-3);
size_box = str2num(size_box);

%extract the number of cells per dimension

nc = spec_arr(2);
nc = char(nc);
nc = nc(1:end-1);
nc = str2num(nc);
 
%extract the number of particle per dimension

np = spec_arr(3);
np = char(np);
np = np(1:end-1);
np = str2num(np);

%extract the initial redshift of the simulation
  
 zi = spec_arr(4);
 zi = char(zi);
 zi = zi(3:end);
 zi = str2num(zi);
 
  % extracts the informations of the wake if there is one
 
 wake_spec = spec_arr(5);
 wake_spec = char(wake_spec);
 if wake_spec(1)=='n'
    wake_or_no_wake='no wake';
    multiplicity_of_files=wake_spec(end); 
    Gmu=0;
    ziw=0;
 end
 if wake_spec(1)=='w'
     wake_or_no_wake='wake';
     wake_spec2=strsplit(wake_spec,{'u','t10m','zi'},'CollapseDelimiters',true);
     Gmu=str2num(char(wake_spec2(2)))*10^(-str2num(char(wake_spec2(3))));
     ziw=char(wake_spec2(4));
     ziw=str2num(ziw(1:end-1));
     multiplicity_of_files=char(wake_spec(end));
 end
 
    

path_total_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj/');
path_analysis_out=strcat(strcat(root_per_node_out,spec,aux_path),'data/',aux_path_per_node_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv','/','stat/box_statistics/dm/dc_all_nodes_1dproj_analysis/');



files_list = dir(strcat(path_total_out,'*','_NSIDE',num2str(NSIDE),'.txt'));
sorted_files_list={files_list.name};

cd('../../processing');

angles = dlmread(strcat('../../python/angles',num2str(NSIDE),'.txt'));
[angle_nuple,number_of_angle_nuple] = size(angles);

        theta=angles(1,:);
        phi=angles(2,:);

sorted_files_list=sort_nat(sorted_files_list);

[aux1 aux2] = size(num2str(NSIDE));
aux3=aux2+22;
redshift_list=cellfun(@(x) x(19:end-aux3),sorted_files_list,'UniformOutput', false);

cd('../preprocessing');

for rds = 1 : length(redshift_list)  
    
    dmdc=dlmread(strcat(path_total_out,'_1dproj_dc_angle_z',char(redshift_list(rds)),'_total_nodes','_NSIDE',num2str(NSIDE),'.txt'));
    
    %analys=count;
    
    dmdc=transpose(dmdc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %computes the peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxi=max(dmdc);
    maxi=transpose(maxi);
    
    mkdir(path_analysis_out,strcat('max_1dproj/'));
    dlmwrite(strcat(path_analysis_out,'max_1dproj/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_max_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi,'delimiter','\t');

% %     % histogram fo the peaks 
%     
% %     maxi=transpose(maxi);
%     
%     bins=[0:0.1:10];
% 
%     [hist_max edges mid loc] = histcn(maxi,bins);
%     hist_max=transpose(hist_max);
% %     
% %     bins=[0:0.001:10];
% %     non_zero_ind_of_maxi=find(~hist_max(1,:)==0);
% %     
% %      fig=figure('Visible', 'off');
% %     hold on;
% %     
% %     %xlim([0 4]);
% %     bins(end)=[];
% %     plot(bins(non_zero_ind_of_maxi(1:end)),hist_max(non_zero_ind_of_maxi(1:end)),'DisplayName',strcat('z = ',char(redshift_list(rds))),'LineWidth',2);
% %     legend('show');
% %     xlabel('maximum density', 'interpreter', 'latex', 'fontsize', 20);
% %     ylabel('number count', 'interpreter', 'latex', 'fontsize', 20);
% %     set(gca,'FontName','FixedWidth');
% %     set(gca,'FontSize',16);
% %     set(gca,'linewidth',2);
% %     
% %    title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E')),strcat('skewness',num2str(skewness(hist_max(1,:)))),strcat('kurtosis',num2str(kurtosis(hist_max(1,:))))},'interpreter', 'latex', 'fontsize', 20);
% % %     title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
% % 
% %     mkdir(path_analysis_out,strcat('max_1dproj/histogram_zoom/'));
% %     path_file_out=strcat(path_analysis_out,'max_1dproj/histogram_zoom/','_',num2str(rds),'_histogram_z',char(redshift_list(rds)),'.png');
% %     saveas(fig,path_file_out);
% %     
% %     hold off;
% %     
%          fig=figure('Visible', 'off');
%     hold on;
%     
%     %fix range of the above
%     %xlim([0 4]);
% %     bins=[0:0.1:10];
%     bins(end)=[];
% %     axis([0 10 0 50000]);
%     plot(bins(1,:),hist_max(1,:),'LineWidth',2);
%     legend(strcat('z = ',char(redshift_list(rds))),strcat('skewness',num2str(skewness(hist_max(1,:)))));
%     xlabel('maximum density', 'interpreter', 'latex', 'fontsize', 20);
%     ylabel('number count', 'interpreter', 'latex', 'fontsize', 20);
%     set(gca,'FontName','FixedWidth');
%     set(gca,'FontSize',16);
%     set(gca,'linewidth',2);
%     
%    title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E')),strcat('skewness',num2str(skewness(hist_max(1,:)))),strcat('kurtosis',num2str(kurtosis(hist_max(1,:))))},'interpreter', 'latex', 'fontsize', 20);
% %     title({strcat('Maximum density of the '),strcat('1d projection for $G\mu=$ ',num2str(Gmu,'%.1E'))},'interpreter', 'latex', 'fontsize', 20);
% 
%     mkdir(path_analysis_out,strcat('max_1dproj/histogram/'));
%     path_file_out=strcat(path_analysis_out,'max_1dproj/histogram/','_',num2str(rds),'_histogram_z',char(redshift_list(rds)),'.png');
%     saveas(fig,path_file_out);
%     
%     hold off;

    
      % % signal to noise stat
  
% %     %computes the histogram
% %     
%   %  cd('../../preprocessing');
% %     
%    %  bins=[-1:0.01:10];
%   %   for i = 1:12*NSIDE^2
% %     
%  %    [result(i) edges mid loc] = histcn(count(:,i),bins);
% %     
% %     sk(:,i) = skewness(count_1d(:,i));
% %     kur(:,i)= kurtosis(count_1d(:,i));
%             
%         fig=figure('Visible', 'off');
%         %fig=figure;
%         
%      hold on;

      st = std(dmdc);
      st = transpose(st);

      mkdir(path_analysis_out,strcat('sigma/'));  
      dlmwrite(strcat(path_analysis_out,'sigma/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_std_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),st,'delimiter','\t');
      
      mkdir(path_analysis_out,strcat('sig_to_noise/'));  
      dlmwrite(strcat(path_analysis_out,'sig_to_noise/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_sig_to_noise_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi(:,1)./st(:,1),'delimiter','\t')
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %computes the pleak, noise and sig to noise with wavelet filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [dc_cwt1,periods] = cwt(dmdc(:,1),seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);
    wavelet_coef= size(dc_cwt1,1);
    angles_n=size(dmdc,2);
    bin_n=size(dmdc,1);
    
    dc_cwt=zeros(wavelet_coef,bin_n,angles_n);
    i_dc_cwt=zeros(bin_n,angles_n);
    
    for i=1:angles_n
    [dc_cwt(:,:,i),periods] = cwt(dmdc(:,i),seconds(size_box/(np*resol_factor)),'waveletparameters',[3 3.01]);   
    i_dc_cwt(:,i) = icwt(dc_cwt(:,:,i),periods,[periods(1) seconds(cutoff)],'waveletparameters',[3 3.01]);    
    end
    
    
    
    maxi_wf=max(i_dc_cwt);
    maxi_wf=transpose(maxi_wf);
    
    mkdir(path_analysis_out,strcat('max_1dproj/','wav_filter_',num2str(cutoff),'Mpc_cut/'));
    dlmwrite(strcat(path_analysis_out,'max_1dproj/','wav_filter_',num2str(cutoff),'Mpc_cut/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_max_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi_wf,'delimiter','\t');

    st_wf = std(i_dc_cwt);
      st_wf = transpose(st_wf);

      mkdir(path_analysis_out,strcat('sigma/','wav_filter_',num2str(cutoff),'Mpc_cut/'));  
      dlmwrite(strcat(path_analysis_out,'sigma/','wav_filter_',num2str(cutoff),'Mpc_cut/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_std_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),st_wf,'delimiter','\t');
      
      mkdir(path_analysis_out,strcat('sig_to_noise/','wav_filter_',num2str(cutoff),'Mpc_cut/'));  
      dlmwrite(strcat(path_analysis_out,'sig_to_noise/','wav_filter_',num2str(cutoff),'Mpc_cut/','_',num2str(rds,strcat('%0',num2str(1+floor(length(redshift_list)/10)),'d')),'_1dproj_sig_to_noise_z',char(redshift_list(rds)),'_NSIDE',num2str(NSIDE),'.txt'),maxi(:,1)./st_wf(:,1),'delimiter','\t')

    
end

cd('../wake_detection/box_statistics');

end

