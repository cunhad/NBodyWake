function [  ] = deltasq_CUBEP3M( root,root_out,spec,aux_path,aux_path_out )

%(example) deltasq_CUBEP3M('/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','40Mpc_240c_120p_zi65_nowakes','/','' );

path_output_delta=strcat(root,spec,aux_path);
files_list = dir(strcat(path_output_delta,'*new.dat'));
sorted_files_list={files_list.name};
cd('../processing');
sorted_files_list=sort_nat(sorted_files_list);

cd('../preprocessing');

for k = 1 :  length(sorted_files_list)
    %for k = 1 : 1
    
    fig=figure('Visible', 'off');
    filename=char(sorted_files_list(k));
    dat_output = dlmread(strcat(path_output_delta,filename));
    %[ dat_output ] = import_( strcat(path_output_delta,filename), '%f %f %f $f $f',5 );
    z = filename(1:end-15);
    

%form the dimensionless power spectrum

%dat_input(:,2)=(1/(2*pi^2))*(dat_input(:,1).^3).*dat_input(:,2);

%semilogx(dat_input(:,1),dat_input(:,2),'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
loglog(dat_output(:,1),dat_output(:,2),'DisplayName',strcat('z = ',z),'LineWidth',2);


title({strcat('Dimensionless power spectrum'),strcat(' from the Nbody Simulation')},'interpreter', 'latex', 'fontsize', 20);
ylabel('$\Delta^2\ (k)$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$k (Mpc^{-1})$', 'interpreter', 'latex', 'fontsize', 20);
legend('show');

mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));

path_out=strcat(strcat(root_out,spec,aux_path),'plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/delta_sq/');
mkdir(strcat(root_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/delta_sq/'));


path_file_out=strcat(path_out,'_',num2str(k),'_deltasq_z',z,'.jpg');
saveas(fig,path_file_out);



end

cd('../power_spectrum');


end
