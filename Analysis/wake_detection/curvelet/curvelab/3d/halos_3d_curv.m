function [ cell_bins1d_x,cell_bins1d_y, cell_bins1d_z , Pos_halos ,mass_halos ,radius_halos, count_sum_h_number,list] = halos_3d_curv( root,root_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,pivot,rot_angle)
%Computes the 2d projections aconding to the input specifications and stores (and returns) the resulting data


%   (example) halos_3d_curv( '/home/asus/Dropbox/extras/storage/graham/ht/', '/home/asus/Dropbox/extras/storage/graham/ht/data/','4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m','/sample3001/half_lin_cutoff_half_tot_pert_nvpw/','','3.000halo0.dat',1,1/4,[0,0,0],[pi/2,0,0]);

% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct3d'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab-2.1.3/fdct_wrapping_matlab'));
 
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));

cd('../../../../preprocessing');


%[ nodes_list redshift_list ] = preprocessing_many_nodes(root,spec,aux_path );
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path ); 

% [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos_h mass Radiusd halos] = preprocessing_halo_nodes( root,spec,aux_path,filename);
[ size_box ,nc ,np ,zi ,wake_or_no_wake ,multiplicity_of_files ,Gmu ,ziw ,z ,path_file_in, ~ ,~, ~ ] = preprocessing_halo_part( root,spec,aux_path,filename,length(nodes_list),1);


%to particle mass

% mass=mass*(np/nc)^3;

cell_bins1d_x=[(np/2)-(np/(2*lenght_factor))+pivot(1)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(1)*(np*resol_factor)/(nc)];
cell_bins1d_y=[(np/2)-(np/(2*lenght_factor))+pivot(2)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(2)*(np*resol_factor)/(nc)];
cell_bins1d_z=[(np/2)-(np/(2*lenght_factor))+pivot(3)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(3)*(np*resol_factor)/(nc)];
cell_bins1d_x(end)=[];
cell_bins1d_y(end)=[];
cell_bins1d_z(end)=[];

count_sum_h_number=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z),numel(cell_bins1d_x));
count_sum_h_mass=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z),numel(cell_bins1d_x));
% count_sum=zeros(numel(cell_bins1d_y),numel(cell_bins1d_z));

psi=rot_angle(3);
theta=rot_angle(2);
phi=rot_angle(1);

Pos_halos=[];
mass_halos=[];
radius_halos=[];

path_in=strcat(root,spec,aux_path);


% for part = 1 : 1
for part = 1 : length(nodes_list)
    
    cell_bins1d_x=[(np/2)-(np/(2*lenght_factor))+pivot(1)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(1)*(np*resol_factor)/(nc)];
    cell_bins1d_y=[(np/2)-(np/(2*lenght_factor))+pivot(2)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(2)*(np*resol_factor)/(nc)];
    cell_bins1d_z=[(np/2)-(np/(2*lenght_factor))+pivot(3)*(np*resol_factor)/(nc):1/resol_factor:(np/2)+(np/(2*lenght_factor))+pivot(3)*(np*resol_factor)/(nc)];

    
%     file_name = dir(strcat(path_in,num2str(z),'*halo',char(nodes_list(node)),'.dat'));
%     filename=file_name.name;
    display(part);
    
    [ ~, ~ ,~, ~, ~, ~, ~, ~, ~, ~, Pos_h ,mass ,Radiusd] = preprocessing_halo_part( root,spec,aux_path,filename,length(nodes_list),part);
%     filename_pos=filename;
%     [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( root,spec,aux_path,strcat(filename_pos(1:strfind(filename,'halo')-1),'xv',filename_pos(strfind(filename,'halo')+4:end)));
    
    %to particle mass
    
    %mass=mass*(np/nc)^3;
    %mass=mass/2;
    
    Pos_h=mod(Pos_h,nc); 
    Pos_h=Pos_h*(np*resol_factor)/(nc);
%     Pos_h(1,:)=Pos_h(1,:)-(nc/2)-pivot(1);
%     Pos_h(2,:)=Pos_h(2,:)-(nc/2)-pivot(2);
%     Pos_h(3,:)=Pos_h(3,:)-(nc/2)-pivot(3);
%     
    Pos_h(1,:)=Pos_h(1,:)-(np*resol_factor/2)-pivot(1)*(np*resol_factor)/(nc);
    Pos_h(2,:)=Pos_h(2,:)-(np*resol_factor/2)-pivot(2)*(np*resol_factor)/(nc);
    Pos_h(3,:)=Pos_h(3,:)-(np*resol_factor/2)-pivot(3)*(np*resol_factor)/(nc);
    
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rx = [ 1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    Rz = [ cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0; 0 0 1];

%     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
%     Rx = [ 1 0 0; 0 cos(psi) -sin(psi); 0 sin(psi) cos(psi)];
%     
%     Pos_h=Rz*Pos_h;
%     Pos_h=Ry*Pos_h;
%     Pos_h=Rx*Pos_h;
    
%     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%     Rx = [ 1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
%     Rz = [ cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0; 0 0 1];
    

    
    R1=Ry*Rx;
    
    R=Rz*R1;
    
%     Pos=Rz*Pos;
%     Pos=Rx*Pos;
%     Pos=Ry*Pos;
    Pos_h=R*Pos_h;
    
%     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
%     Pos_h=Rz*Pos_h;
%     Pos_h=Ry*Pos_h;
%     %Pos=Rz*Pos;
    
    liminf=-(1/(2*lenght_factor))*np;
    limsup= (1/(2*lenght_factor))*np;
    conditionsx=Pos_h(1,:)<=liminf|Pos_h(1,:)>=limsup;
    conditionsy=Pos_h(2,:)<=liminf|Pos_h(2,:)>=limsup;
    conditionsz=Pos_h(3,:)<=liminf|Pos_h(3,:)>=limsup;
    conditions=conditionsx|conditionsy|conditionsz;
    Pos_h(:,conditions)=[];
    mass(conditions)=[];
    Radiusd(conditions)=[];
    
    Pos_h(1,:)=Pos_h(1,:)+(np*resol_factor/2)+pivot(1)*(np*resol_factor)/(nc);
    Pos_h(2,:)=Pos_h(2,:)+(np*resol_factor/2)+pivot(2)*(np*resol_factor)/(nc);
    Pos_h(3,:)=Pos_h(3,:)+(np*resol_factor/2)+pivot(3)*(np*resol_factor)/(nc);
        
    Pos_halos=[ Pos_halos , Pos_h];  
    mass_halos=[mass_halos , mass];
    radius_halos=[radius_halos , Radiusd];
        
    if (~isempty(Pos_h))
        
        
        [count_h_n edges mid loc] = histcn(transpose(Pos_h),cell_bins1d_y,cell_bins1d_z,cell_bins1d_x);
%         count_h_n=count_h_n(1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1,1:1);
        %     average=mean2(count);
        %     count=(count-average)/average;
%         count_h_n=squeeze(count_h_n);
        
        %    cell_bins1d(end)=[];
        count_sum_h_number=count_sum_h_number+count_h_n;
        
        
        [count_h_m edges mid loc] = histcn(transpose(Pos_h),cell_bins1d_y,cell_bins1d_z,cell_bins1d_x,'AccumData',transpose(mass));
%         count_h_m=count_h_m(1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1,1:1);
        %     average=mean2(count);
        %     count=(count-average)/average;
%         count_h_m=squeeze(count_h_m);
        
        %    cell_bins1d(end)=[];
        count_sum_h_mass=count_sum_h_mass+count_h_m;  
       
        
    end
    
 
    
%     Pos=mod(Pos,nc);    
%     Pos(1,:)=Pos(1,:)-(nc/2)-pivot(1);
%     Pos(2,:)=Pos(2,:)-(nc/2)-pivot(2);
%     Pos(3,:)=Pos(3,:)-(nc/2)-pivot(3);
%     
%     Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
%     Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
%     
%     Pos=Rz*Pos;
%     Pos=Ry*Pos;
%     %Pos=Rz*Pos;
%     
%     liminf=-(1/(2*lenght_factor))*nc;
%     limsup= (1/(2*lenght_factor))*nc;
%     conditionsx=Pos(1,:)<=liminf|Pos(1,:)>=limsup;
%     conditionsy=Pos(2,:)<=liminf|Pos(2,:)>=limsup;
%     conditionsz=Pos(3,:)<=liminf|Pos(3,:)>=limsup;
%     conditions=conditionsx|conditionsy|conditionsz;
%     Pos(:,conditions)=[];
%     
%     Pos(1,:)=Pos(1,:)+(nc/2)+pivot(1);
%     Pos(2,:)=Pos(2,:)+(nc/2)+pivot(2);
%     Pos(3,:)=Pos(3,:)+(nc/2)+pivot(3);
%     
%     Pos=transpose(Pos);
%     
%     [count edges mid loc] = histcn(Pos,1,cell_bins1d_y,cell_bins1d_z);
%     count=count(1:1,1:numel(cell_bins1d_y)-1,1:numel(cell_bins1d_z)-1);
%     %     average=mean2(count);
%     %     count=(count-average)/average;
%     count=squeeze(count);
%     
%     %    cell_bins1d(end)=[];
%     count_sum=count_sum+count;
  
    
end

    if (isempty(Pos_h))
        
        Pos_halos=0;  
        mass_halos=0;
        radius_halos=0;
        
    end
    
       
    
% F2 = ones(numel(cell_bins1d_y)-1,numel(cell_bins1d_z)-1,numel(cell_bins1d_x)-1);
% X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
% C2 = fdct3d_forward(X2);
% E2 = cell(size(C2));
% for s=1:length(C2)
%     E{s} = cell(size(C2{s}));
%     for w=1:length(C2{s})
%         A2 = C2{s}{w};
%         E2{s}{w} = sqrt(sum(sum(A2.*conj(A2))) / prod(size(A2)));
%     end
% end
% 
% 
% 
% F2=zeros(numel(cell_bins1d_y)-1,numel(cell_bins1d_z)-1,numel(cell_bins1d_x)-1);
% C_zero2 = fdct3d_forward(F2);

count_sum_h_number=count_sum_h_number+1;
count_sum_h_number=log(count_sum_h_number);

C = fdct3d_forward(count_sum_h_number);
% coef=real(C{6}{1});
% coef=coef(:);
% 
% fig=figure;
% kurtosis(coef(:))
% histogram(coef(:))
% title({strcat('max=',num2str(max(coef(:)))),strcat('kurtosis = ',num2str(kurtosis(coef(:))))})
% 
% 
% coef(coef==0)=[];
% fig=figure;
% kurtosis(coef(:))
% histogram(coef(:))
% title({strcat('max=',num2str(max(coef(:)))),strcat('kurtosis = ',num2str(kurtosis(coef(:))))})
% 
% fig=figure;
% coef(coef<=1e-2)=[];
% histogram(coef(:))
% title({strcat('max=',num2str(max(coef(:)))),strcat('kurtosis = ',num2str(kurtosis(coef(:))))})


list=cell(length(C),1);

% for s = 1:length(C)
%     for w = 1:length(C{s})
% 
%         list{s} = [list{s};abs(C{s}{w}(:))];
%     end
%     fig=figure;
%     histogram(coef(:))
%     title({strcat('level=',num2str(s),strcat('max=',num2str(max(list{s}))),strcat('kurtosis = ',num2str(kurtosis(list{s})))})
% end

    spec
    aux_path

for s = length(C):-1:1
    for w = 1:length(C{s})
%         w/length(C{s})
        list{s} = [list{s};(C{s}{w}(:))];
    end
%     fig=figure;
%     histogram(abs(list{s}))
    title({strcat('level=',num2str(s)),strcat('max=',num2str(max(abs(list{s})))),strcat('kurtosis = ',num2str(kurtosis(abs(list{s}))))})
    strcat('level=',num2str(s))
    strcat('abs max=',num2str(max(abs(list{s}))))
    strcat('abs kurtosis = ',num2str(kurtosis(abs(list{s}))))
    strcat('max=',num2str(max([real(list{s}),imag(list{s})])))
    strcat('kurtosis = ',num2str(kurtosis([real(list{s}),imag(list{s})])))
end


%to particle mass unit


% average=mean2(count_sum_h_number);
% count_sum_dmdc=(count_sum_h_mass-average)/average;

% mkdir(root_out);
% mkdir(root_out,strcat(spec,aux_path));
% 
% path_out=strcat(strcat(root_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/');
% mkdir(strcat(root_out,spec,aux_path),strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf_',strcat(num2str(pivot(1)),'-',num2str(pivot(2)),'-',num2str(pivot(3))),'pv_',strcat(num2str(rot_angle(1)),'-',num2str(rot_angle(2))),'ra','/','2dproj/halos/'));
% 
% %mkdir(path_out,'aux/total_dm/');
% %dlmwrite(strcat(path_out,'aux/total_dm/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_dm_z',num2str(z),'_data.txt'),count_sum,'delimiter','\t');
% mkdir(path_out,strcat('number/'));
% % dlmwrite(strcat(path_out,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_n_z',num2str(z),'_data.txt'),count_sum_h_number,'delimiter','\t');
% mkdir(path_out,'mass/');
% % dlmwrite(strcat(path_out,'mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_z',num2str(z),'_data.txt'),count_sum_h_mass,'delimiter','\t');
% %mkdir(path_out,'h_mass_dc_wrt_dmavr/');
% %dlmwrite(strcat(path_out,'h_mass_dc_wrt_dmavr/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_dmdc_z',num2str(z),'_data.txt'),count_sum_dmdc,'delimiter','\t');
% 
%     for count_slice=1:slice
%         
%             fileID = fopen(strcat(path_out,'number/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_n_z',num2str(z),'_data_sl',num2str(count_slice),'.bin'),'w');
%             fwrite(fileID,count_sum_h_number(:,:,count_slice), 'float32','l');
%             fclose(fileID);
%             
%             fileID = fopen(strcat(path_out,'mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_halos_mass_z',num2str(z),'_data_sl',num2str(count_slice),'.bin'),'w');
%             fwrite(fileID,count_sum_h_mass(:,:,count_slice), 'float32','l');
%             fclose(fileID);
%     end
% 
% mkdir(path_out,'aux/Pos_h/');
% dlmwrite(strcat(path_out,'aux/Pos_h/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_Pos_halos_z',num2str(z),'_data.txt'),Pos_halos,'delimiter','\t');
% mkdir(path_out,'aux/mass/');
% dlmwrite(strcat(path_out,'aux/mass/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_mass_halos_z',num2str(z),'_data.txt'),mass_halos,'delimiter','\t');
% mkdir(path_out,'aux/radius/');
% dlmwrite(strcat(path_out,'aux/radius/','_',num2str(find(str2num(char(redshift_list))==z)),'_2dproj_radius_halos_z',num2str(z),'_data.txt'),radius_halos,'delimiter','\t');


cd('../wake_detection/curvelet/curvelab/3d/');

end



