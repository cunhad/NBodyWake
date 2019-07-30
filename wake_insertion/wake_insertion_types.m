function [  ] = wake_insertion_types( path,spec,aux_path ,z_insert,Gmu_insert,type)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%(example ) wake_insertion_types('/home/asus/Dropbox/extras/storage/graham/small_res/','64Mpc_96c_48p_zi255_nowakem','/sample1001/',10,6E-6,1);




if type==0
   type_folder='' 
end

if type==1
   type_folder='/test/' 
end

if type==2
   type_folder='/no_vpert_in_wake_hard/' 
end

if type==3
    type_folder='/no_vpert_in_wake/'
end

if type==4
    
    type_folder='/half_lin_cutoff_half_tot_pert/'
    
end

if type==5
    
    type_folder='/quarter_lin_cutoff_half_tot_pert/'
    
end

%combination of 3 and 4

if type==6
    
    type_folder='/half_lin_cutoff_half_tot_pert_nvpwh/'
    
end

%combination of 3 and 5

if type==7
    
    type_folder='/quarter_lin_cutoff_half_tot_pert_nvpwh/'
    
end

if type==8
    
    type_folder='/half_lin_cutoff_half_tot_pert_nvpw/'
    
end

%combination of 3 and 5

if type==9
    
    type_folder='/quarter_lin_cutoff_half_tot_pert_nvpw/'
    
end



cd('../parameters')

[ vSgammaS displacement vel_pert] = wake( Gmu_insert ,z_insert);

[ h OmegaBM OmegaCDM OmegaM OmegaL clight zi t_0 Hzero tensor_tilt spectral_indice sigma8 T_cmb_t0 Scalar_amplitude ] = cosmology(  );


cd('../Analysis/preprocessing');

% [ nodes_list redshift_list  ] = preprocessing_many_nodes(path,spec,aux_path );
[~,redshift_list,nodes_list,~,~,~,~,~,~,~,~] = preprocessing_info(path,spec,aux_path );


path_in=strcat(path,spec,aux_path);
files_list = dir(strcat(path_in,num2str(z_insert,'%.3f'),'*','xv*.dat'));
files_list={files_list.name};
filename=cell2mat(files_list(1));

files_list_pid = dir(strcat(path_in,num2str(z_insert,'%.3f'),'*','PID*.dat'));
files_list_pid={files_list_pid.name};

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( path,spec,aux_path,filename);

LinConv=nc/size_box; %spatial linear convertion factor from comoving to simulation
TimConv=(3*((1+z)^(2))*Hzero*(OmegaM^(1/2)))/2;  %Convert time in simulation units to seconds


cd('../processing');

files_list=sort_nat(files_list);
files_list_pid=sort_nat(files_list_pid);


cd('../preprocessing');

Gmu_exp=floor(log10(Gmu_insert));
Gmu_coef=Gmu_insert/(1*10^Gmu_exp);

spec_out=strcat(num2str(size_box),'Mpc_',num2str(nc),'c_',num2str(np),'p_','zi',num2str(zi),'_wakeGmu',num2str(Gmu_coef),'t10m',num2str(-Gmu_exp),'zi',num2str(z_insert),multiplicity_of_files,aux_path,type_folder);
mkdir(path,spec_out);


    %for node = 1 : 1
    for node = 1 : length(nodes_list)
        
        %read data
%         fid = fopen(strcat(path_in,cell2mat(files_list(node))));
%         directory = dir(strcat(path_in,cell2mat(files_list(node))));
%         particles=(directory.bytes-48)/24;
%         headear = fread(fid, [12 1], 'float32','l') ;
%         data=fread(fid, [6 particles], 'float32','l');

        filename=cell2mat(files_list(node));
        [ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in header Pos Vel i_node j_node k_node number_node_dim] = preprocessing_nodes_all( path,spec,aux_path,filename);
        data=Pos;
         data(4:6,:)=Vel;
%         
        %displace towards the wake
        
        dist_to_wake3=data(3,:)-nc/2; %is the vector that points to the wake at Z=nc/2 plane
        
        %if type=2(no_vel_pert_in_wake_hard) do
        
%         if type==2|type==6|type==7
%             
%             dist_to_wake3((dist_to_wake3<1)&(dist_to_wake3>-1))=0;
%             
%         end
%         
%         %if type=3(no_vel_pert_in_wake) do
%         
%         if type==3|type==8|type==9
%             
%             dist_to_wake3((dist_to_wake3<displacement*LinConv)&(dist_to_wake3>-displacement*LinConv))=0;
%             
%         end
        
       
        
        displacement_to_wake3=-sign(dist_to_wake3(1,:))*displacement*LinConv; %the particles will be displaced towards the wake
        
        %if type=4(half_lin_cutoff_half_tot_pert) do
        
        if type==4|type==6|type==8
            
            cuttoff=ones(1,length(dist_to_wake3));
            
            cuttoff(dist_to_wake3<-nc/4)=(dist_to_wake3(dist_to_wake3<-nc/4)+nc/2)/(nc/4);
            cuttoff(dist_to_wake3>nc/4)=(-dist_to_wake3(dist_to_wake3>nc/4)+nc/2)/(nc/4);            
            displacement_to_wake3=displacement_to_wake3.*cuttoff;
        end
        
        %if type=5(quarter_lin_cutoff_half_tot_pert) do
        
        if type==5|type==7|type==9
            
            cuttoff=ones(1,length(dist_to_wake3));            
            cuttoff(dist_to_wake3<=-3*nc/8)=0;
            cuttoff(dist_to_wake3>-3*nc/8&dist_to_wake3<-nc/4)=(dist_to_wake3(dist_to_wake3>-3*nc/8&dist_to_wake3<-nc/4)+3*nc/8)/(nc/8);
            cuttoff(dist_to_wake3<3*nc/8&dist_to_wake3>nc/4)=(-dist_to_wake3(dist_to_wake3<3*nc/8&dist_to_wake3>nc/4)+3*nc/8)/(nc/8);
            cuttoff(dist_to_wake3>=3*nc/8)=0;
            displacement_to_wake3=displacement_to_wake3.*cuttoff;
        end
        
        data(3,:)=data(3,:)+displacement_to_wake3(1,:);
        
        %transform back to local coordinates
        
        data(1,:)=data(1,:)-(nc/number_node_dim)*i_node;
        data(2,:)=data(2,:)-(nc/number_node_dim)*j_node;
        data(3,:)=data(3,:)-(nc/number_node_dim)*k_node;
        

        %give hte velocity kick
        
        kick_to_wake3=-sign(dist_to_wake3(1,:))*vel_pert*LinConv/TimConv;
        
        %if type=4(half_lin_cutoff_half_tot_pert) do
        
        if type==2|type==6|type==7
            
            kick_to_wake3((dist_to_wake3<1)&(dist_to_wake3>-1))=0;
            
        end
        
        %if type=3(no_vel_pert_in_wake) do
        
        if type==3|type==8|type==9
            
            kick_to_wake3((dist_to_wake3<2*displacement*LinConv)&(dist_to_wake3>-2*displacement*LinConv))=0;
            
        end

        
        if type==4|type==6|type==8
            
            kick_to_wake3=kick_to_wake3.*cuttoff;
            
        end
        
        %if type=5(quarter_lin_cutoff_half_tot_pert) do
        
        if type==5|type==7|type==9
            
           kick_to_wake3=kick_to_wake3.*cuttoff;
           
        end
        
        
        data(6,:)=data(6,:)+kick_to_wake3(1,:);
        
        %write data
        file_out = fopen(strcat(path,spec_out,cell2mat(files_list(node))),'w');
        fwrite(file_out,header,'float32','l');
        fwrite(file_out,data,'float32','l');
        fclose(file_out);
        
        %copy pids
        
        
        
        copyfile(strcat(path,spec,aux_path,cell2mat(files_list_pid(node))),strcat(path,spec_out,cell2mat(files_list_pid(node)))) ;
        
        
%         file_name = dir(strcat(path_in,char(redshift_list(rds)),'xv',char(nodes_list(node)),'.dat'));
%         filename=file_name.name;
%         [ size_box nc zi wake_or_no_wake multiplicity_of_files Gmu ziw z path_file_in Pos ] = preprocessing_nodes( path,spec,aux_path,filename,percentage_analysed);
    end
    





cd('../../wake_insertion');

end

