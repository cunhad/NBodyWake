function [  ] = b_0_001t_box_statistics_halos_withpartmass( path,spec,aux_path,NSIDE,node ,num_cores)
    
%(example) b_0_001t_box_statistics_halos_withpartmas('/home/asus/Dropbox/extras/storage/','40Mpc_192c_zi65_nowakes','/',4,0,4);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



p = parpool(num_cores);
tic;

path_in=strcat(path,spec,aux_path);
files_list = dir(strcat(path_in,'*halo',num2str(node),'.dat'));
sorted_files_list={files_list.name};

cd('../../processing');

sorted_files_list=sort_nat(sorted_files_list);

%display(sorted_files_list)

cd('../preprocessing');

angles = dlmread(strcat('../../python/angles',num2str(NSIDE),'.txt'));
[angle_nuple,number_of_angle_nuple] = size(angles);

number_of_redshifts=length(sorted_files_list);



%for k = 3:-1:1
for k = 1  :   length(sorted_files_list)
%for k = 1  : 1
    
    cd('../preprocessing');
    
    %proj=[];
    
    filename=char(sorted_files_list(k));
    
    [ size_box nc zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos mass Radiusd halos] = preprocessing_halo_nodes(path,spec,aux_path,filename,1.0 );
   % display(strcat(path,spec,aux_path,filename));
    %     [ size_box nc zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( path,spec,aux_path,filename,1.0);
    if (halos~=0)
    Pos=mod(Pos,nc);
    
    [m n] = size(Pos);
    
    int_div=fix(m/num_cores);
     
    for i=1:num_cores+1
        counter(i)=((i-1)*int_div);
    end
        counter(num_cores+2)=m;
        
        Pos=transpose(Pos);
    
    bins=[-nc/4:2:nc/4];
    
    pivot=[0,0,0]; %this is the position od the origin of the rotation point with respect to the center of the box
    
   proj1d_angles=zeros(number_of_angle_nuple,length(bins)-1);
    
    parfor i=1:number_of_angle_nuple
%     for i=1:number_of_angle_nuple
 %   for i=1:1
         

         
        theta=angles(1,i);
        phi=angles(2,i);
        
        hist1d_cor=zeros(1,length(bins)-1);
        
       for j=1:num_cores+1
       % for j=1:1000
       
       rx=[];
        
        rx(1,:)=Pos(1,counter(j)+1:counter(j+1))-(nc/2)-pivot(1);
        rx(2,:)=Pos(2,counter(j)+1:counter(j+1))-(nc/2)-pivot(2);
        rx(3,:)=Pos(3,counter(j)+1:counter(j+1))-(nc/2)-pivot(3);

        Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]; 
%         
        rx=Rz*rx;
        rx=Ry*rx;      
        
        liminf=-(1/4)*nc;
        limsup= (1/4)*nc;
        conditionsx=rx(1,:)<=liminf|rx(1,:)>=limsup;
        conditionsy=rx(2,:)<=liminf|rx(2,:)>=limsup;
        conditionsz=rx(3,:)<=liminf|rx(3,:)>=limsup;
        conditions=conditionsx|conditionsy|conditionsz;
        rx(:,conditions)=[];
        
         rx=transpose(rx);
         
         %display(rx);
         
         if(~isempty(rx))
        
        [count edges mid loc] = histcn(rx,1,1,bins,'AccumData',transpose(mass));   %computes the halo mas in each bin, in grid units
       % display(count);
       % display(length(bins));
        count=count(1:1,1:1,1:length(bins)-1);
   %     average=mean2(count);
   %     count=(count-average)/average;
        count=squeeze(count);
        count=squeeze(count);
        
        count=transpose(count);
        
       % display(count);
        
        
        
        hist1d_cor=hist1d_cor+count;
        
        %display(hist1d_cor);
      
        
        
         end
         
%          display(hist1d_cor);
%          display(mass);
        

       end


       
       % proj=transpose(proj);
        
        %display(proj);
        
        % [proj edges mid loc] = histcn(proj,bins);
         
        
        %proj=transpose(proj);
        
       % display(size(proj1d_angles));
       % display(size(hist1d_cor));
       
        % average=mean2(hist1d_cor);
        % hist1d_cor=(hist1d_cor-average)/average;
         
         proj1d_angles(i,:)=hist1d_cor(:);

        fprintf('done for z= %f and  i= %d\n',z, i);
        %display(proj);
 
    end
    
        % display(proj1d_angles);
        % display(mass); 
     
    path_out=path_in;
    %for guillimin
%     path_out=strcat('/gs/scratch/cunhad/',spec,aux_path);
    
     mkdir(path_out,strcat('data/','stat/box_statistics/halos/'));
      dlmwrite(strcat(path_out,'data/','stat/box_statistics/halos/','1dproj_angle_halos_z',num2str(z),'_node',num2str(node),'_NSIDE',num2str(NSIDE),'.txt'),proj1d_angles,'delimiter','\t');
%     save(strcat(path_in,'Analysis/','stat/box_statistics/','1dproj_angle_z',num2str(z),'_node',num2str(node),'_NSIDE',num2str(NSIDE),'.txt'),'proj1d_angles', '-ascii');
    end
end

cd('../wake_detection/box_statistics');

toc;

delete(gcp('nocreate'))

end

