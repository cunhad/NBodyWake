function [  ] = wavelets_1d_dm_plot( root,root_data_out,root_out,spec,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle,lim)
    
%(example)  qu
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%p = parpool(num_cores);
%tic;


path_analysis_out=strcat(root,spec,aux_path,'Analysis/stat/box_statistics/');

cd('../../preprocessing');

[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );

path_in=strcat(root,spec,aux_path);
filename = dir(strcat(path_in,char(redshift_list(1)),'xv',char(nodes_list(1)),'.dat'));
filename=filename.name;

[ size_box nc zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,filename,1.0);


%for k = 3:-1:1
%for k = 1  :   length(sorted_files_list)
for rds = 1 : length(redshift_list)
    
    bins=[-nc/2:2:nc/2];
    bins(end)=[];
    count_sum=zeros(1,numel(bins));
    
    for node = 1 : length(nodes_list)
        
    
        path_in=strcat(root,spec,aux_path);
        file_name = dir(strcat(path_in,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat'));
        file_name=file_name.name;
        
        
        [ size_box nc zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos] = preprocessing_nodes( root,spec,aux_path,file_name,1.0);
        
        Pos=transpose(Pos);
        
        bins=[-nc/2:2:nc/2];
        
        pivot=[0,0,0]; %this is the position od the origin of the rotation point with respect to the center of the box
        
        
        theta=0;
        phi=0;
        
        rx=[];
        
        %         rx(1,:)=Pos(1,counter(j)+1:counter(j+1))-(nc/2)-pivot(1);
        %         rx(2,:)=Pos(2,counter(j)+1:counter(j+1))-(nc/2)-pivot(2);
        %         rx(3,:)=Pos(3,counter(j)+1:counter(j+1))-(nc/2)-pivot(3);
        
        rx(1,:)=Pos(1,:)-(nc/2)-pivot(1);
        rx(2,:)=Pos(2,:)-(nc/2)-pivot(2);
        rx(3,:)=Pos(3,:)-(nc/2)-pivot(3);
        
        
        Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
        %
        rx=Rz*rx;
        rx=Ry*rx;
        
        liminf=-(1/2)*nc;
        limsup= (1/2)*nc;
        conditionsx=rx(1,:)<=liminf|rx(1,:)>=limsup;
        conditionsy=rx(2,:)<=liminf|rx(2,:)>=limsup;
        conditionsz=rx(3,:)<=liminf|rx(3,:)>=limsup;
        conditions=conditionsx|conditionsy|conditionsz;
        rx(:,conditions)=[];
        
        rx=transpose(rx);
        
        %display(rx);
        
        if(~isempty(rx))
            
            [count edges mid loc] = histcn(rx,1,1,bins);
            % display(count);
            % display(length(bins));
            count=count(1:1,1:1,1:length(bins)-1);
            %     average=mean2(count);
            %     count=(count-average)/average;
            count=squeeze(count);
            count=squeeze(count);
            
            count=transpose(count);
            
            % display(count);
            
            
            count_sum=count_sum+count;
        
        
        
         end
        

        end
                
        average=mean(count_sum);
        count_sum=(count_sum-average)/average;
    
     
    %wavelet analysis
    

    %hold on;
    fig=figure('Visible', 'off');
    
    [mcwt,periods] = cwt(count_sum,seconds(2*size_box/nc),'waveletparameters',[3 3.01]);
%     hp = pcolor( 0:size_box/length(count_sum):size_box-size_box/length(count_sum),seconds(periods),abs(mcwt)); 
%     hp.EdgeColor = 'none';

imagesc(0:size_box/length(count_sum):size_box-size_box/length(count_sum),seconds(periods),abs(mcwt));
colorbar;

set(gca,'YScale','log');
xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20); ylabel('Scale parameter (Mpc)', 'interpreter', 'latex', 'fontsize', 20);
title(strcat('Continuous wavelet transformation (Morse)'),'interpreter', 'latex', 'fontsize', 20);
mkdir(path_analysis_out,strcat('wavelets/mtwt/'));
path_file_out=strcat(path_analysis_out,'wavelets/mtwt/','_',num2str(rds),'_mcwt_z',num2str(z),'.png');
saveas(fig,path_file_out);
hold off;

fig=figure('Visible', 'off');
clims = [0 1];
imagesc(0:size_box/length(count_sum):size_box-size_box/length(count_sum),seconds(periods),abs(mcwt),clims);
colorbar;
set(gca,'YScale','log');
xlabel('$Z(Mpc)$', 'interpreter', 'latex', 'fontsize', 20); ylabel('Scale parameter (Mpc)', 'interpreter', 'latex', 'fontsize', 20);
title(strcat('Continuous wavelet transformation (Morse)'),'interpreter', 'latex', 'fontsize', 20);
mkdir(path_analysis_out,strcat('wavelets/mtwt_lims/'));
path_file_out=strcat(path_analysis_out,'wavelets/mtwt_lims/','_',num2str(rds),'_mcwt_lims_z',num2str(z),'.png');
saveas(fig,path_file_out);
hold off;
    
    
end

cd('../wake_detection/lets');

%toc;

%delete(gcp('nocreate'))

end

