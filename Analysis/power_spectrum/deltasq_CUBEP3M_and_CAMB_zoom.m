function [  ] = deltasq_CUBEP3M_and_CAMB_zoom( root,root_out,spec,aux_path,aux_path_out  )

%(example) deltasq_CUBEP3M_and_CAMB_zoom('/home/asus/Dropbox/extras/storage/','/home/asus/Dropbox/extras/storage/','40Mpc_240c_120p_zi65_nowakes','/','' );

path_output_delta=strcat(root,spec,aux_path);
files_list = dir(strcat(path_output_delta,'*new.dat'));
sorted_files_list={files_list.name};

path_CAMB='../../CAMB/transfer_functions/';
files_list_CAMB = dir(strcat(path_CAMB,'*.dat'));
sorted_files_list_CAMB={files_list_CAMB.name};

cd('../processing');
sorted_files_list=sort_nat(sorted_files_list);
sorted_files_list_CAMB=sort_nat(sorted_files_list_CAMB);

cd('../preprocessing');

for k = 1 :  length(sorted_files_list)
    %for k = 1 : 1
    
    fig=figure('Visible', 'off');

    
    filename=char(sorted_files_list(k));
    dat_output = dlmread(strcat(path_output_delta,filename));
    %[ dat_output ] = import_( strcat(path_output_delta,filename), '%f %f %f $f $f',5 );
    z = filename(1:end-15);
    
    pattern=strcat('z',z);
    logical_op=contains(sorted_files_list_CAMB,pattern);
    filename=char(sorted_files_list_CAMB(logical_op));
    
    %filename=char(sorted_files_list_CAMB(k));
    [ dat_CAMB ] = import_( strcat(path_CAMB,filename), '%f %f ',2 );
    
    

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

%semilogx(dat_input(:,1),dat_input(:,2),'DisplayName',strcat('z = ',num2str(z)),'LineWidth',2);
loglog(dat_output(:,1),dat_output(:,2),'DisplayName','CUBEP3M','LineWidth',2);
    hold on;
loglog(dat_CAMB(:,1),dat_CAMB(:,2),'DisplayName','CAMB','LineWidth',2);

title(strcat('Dimensionless power spectrum at ',strcat(' z = ',z)),'interpreter', 'latex', 'fontsize', 15);
ylabel('$\Delta^2\ (k)$', 'interpreter', 'latex', 'fontsize', 20);
xlabel('$k (Mpc^{-1})$', 'interpreter', 'latex', 'fontsize', 20);
legend('show');
hold off;


mkdir(root_out);
mkdir(root_out,strcat(spec,aux_path));

path_out=strcat(strcat(root_out,spec,aux_path),'plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/deltasq_CAMBcomp_zoom/');
mkdir(strcat(root_out,spec,aux_path),strcat('plot/',aux_path_out,num2str(1),'lf_',num2str(1),'rf_',strcat(num2str(0),'-',num2str(0),'-',num2str(0)),'pv_',strcat(num2str(0),'-',num2str(0)),'ra','/','power_spectrum/deltasq_CAMBcomp_zoom/'));


path_file_out=strcat(path_out,'_',num2str(k),'_deltasq_CAMBcomp_zoom_z',z,'.jpg');
saveas(fig,path_file_out);

end

cd('../../Analysis/power_spectrum');


end