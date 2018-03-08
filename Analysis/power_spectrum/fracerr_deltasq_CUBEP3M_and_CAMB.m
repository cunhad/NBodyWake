function [  ] = fracerr_deltasq_CUBEP3M_and_CAMB( root,root_out,spec,aux_path,aux_path_out )

%(example) fracerr_deltasq_CUBEP3M_and_CAMB('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/check/','64Mpc_96c_48p_zi255_nowakem','/sample1001/','');

cd('../preprocessing');

[~,redshift_list,~,~,~,~,~,~,~,~,~] = preprocessing_info(root,spec,aux_path );

path_output_delta=strcat(root,spec,aux_path);
files_list = dir(strcat(path_output_delta,'*new.dat'));
sorted_files_list={files_list.name};

path_CAMB='../../CAMB/transfer_functions/';
% files_list_CAMB = dir(strcat(path_CAMB,'*.dat'));
% sorted_files_list_CAMB={files_list_CAMB.name};

% cd('../processing');
% sorted_files_list=sort_nat(sorted_files_list);
% sorted_files_list_CAMB=sort_nat(sorted_files_list_CAMB);

% cd('../preprocessing');

for rds = 1 :  length(redshift_list)
    %for k = 1 : 1
    
    fig=figure('Visible', 'off');

    filename_CUBEP3M = strcat(char(redshift_list(rds)),'ngpps_new.dat');
    
%     filename=char(sorted_files_list(rds));
    dat_output = dlmread(strcat(path_output_delta,filename_CUBEP3M));
    %[ dat_output ] = import_( strcat(path_output_delta,filename), '%f %f %f $f $f',5 );
%     z = filename(1:end-15);
    
%     pattern=strcat('z',z);
%     logical_op=contains(sorted_files_list_CAMB,pattern);
    
    
%     filename=char(sorted_files_list_CAMB(logical_op));
      filename_CAMB=strcat('camb_matterpower_z',char(redshift_list(rds)),'.dat');



    %filename=char(sorted_files_list_CAMB(k));
    [ dat_CAMB ] = import_( strcat(path_CAMB,filename_CAMB), '%f %f ',2 );
    
    

%form the dimensionless power spectrum for CAMB

dat_CAMB(:,2)=(1/(2*pi^2))*(dat_CAMB(:,1).^3).*dat_CAMB(:,2);

%take into account h=0.7
%h=0.7;
%dat_CAMB(:,1)=dat_CAMB(:,1)/h;
%dat_input2(:,1)=dat_input(:,1);
%dat_input(:,2)=dat_input(:,2)/(h^3);
%dat_input(:,1)=dat_input(:,1)/h;

%let make then to have the same range (take the nbody simulation as
%standard

liminf=min(dat_output(:,1));
limsup=max(dat_output(:,1));
cond=dat_CAMB(:,1)<=liminf|dat_CAMB(:,1)>=limsup;
dat_CAMB(cond,:)=[];

%now compare the relative separation

comp(:,1)=dat_CAMB(:,1);
dat_output2(:,2)=interp1(dat_output(:,1),dat_output(:,2),dat_CAMB(:,1));
dat_output2(:,1)=comp(:,1);
comp(:,2)=(dat_output2(:,2)-dat_CAMB(:,2))./dat_CAMB(:,2);



semilogx(comp(:,1),comp(:,2));

ylim([-1 1]);

%semilogx(dat_input(:,1),dat_input(:,2),'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
%loglog(dat_output(:,1),dat_output(:,2),'DisplayName','CUBEP3M','LineWidth',2);
    %hold on;
%loglog(dat_CAMB(:,1),dat_CAMB(:,2),'DisplayName','CAMB','LineWidth',2);

%title({strcat('Fractional error of the'),strcat('dimensionless power spectrum at ',strcat(' z = ',char(redshift_list(rds))))},'interpreter', 'latex', 'fontsize', 15);
ylabel('$\frac{\Delta^{2}_{Nbody}-\Delta^{2}_{CAMB}}{\Delta^{2}_{CAMB}}\ (k)$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$k (h/Mpc)$', 'interpreter', 'latex', 'fontsize', 20);
%legend('show');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));

path_out=strcat(strcat(root_out,spec,aux_path),'plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/fractional_error/');
mkdir(strcat(root_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/fractional_error/'));

path_file_out=strcat(path_out,'_',num2str(rds),'_fraceer_deltasq_CAMBcomp_z',char(redshift_list(rds)),'.jpg');
saveas(fig,path_file_out);

%hold off;

end

cd('../../Analysis/power_spectrum');


end