function [ ] = comvel_evol_mem_fast_meanslices_par_all( root_data_in,spec,wake_type)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% (example) []=comvel_evol_mem_fast_meanslices_par_all( '/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_96c_48p_zi255_wakeGmu5t10m7zi63m','')

% root='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/'
root_data_in='/home/asus/Dropbox/extras/storage/graham/ht/data_paper3newRange2/check/'
% root_plot_out='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/plot_paper2/'
spec='16Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
% aux_path='/sample3001/'
wake_type='/half_lin_cutoff_half_tot_pert_nvpw_v0p6/'




type_folder=wake_type;

%aux_path,type_folder,'check/vel/half/'

cd('../../../preprocessing')

% test_particle_id=10000;

[ size_box nc np zi wake_or_no_wake multiplicity_of_files Gmu ziw ] = preprocessing_from_spec( spec);

cd('../wake_detection/consist_check/Zeldovich/')

specs_path_list_nowake=strcat(root_data_in,spec)
sample_list_wake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_wake={sample_list_wake.name};


for sample = 1:length(sample_list_wake)
    
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'/check/vel/half/','_Check_mn_Vel_Zel.txt')
%                 info_med = dlmread(filename_wake);
    
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_mod_Vel_Zel.txt')
%                 info_med = dlmread(filename_wake);
                
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'/check/vel/half/','_Check_med_Vel_Zel.txt')
%                 info_med = dlmread(filename_wake);

%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_mn_Vel_Zel_quart.txt')
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_mn_Vel_Zel.txt')
%                 info_med = dlmread(filename_wake);

%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_mod_Vel_Zel_quart.txt')
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_mod_Vel_Zel.txt')
%                 info_med = dlmread(filename_wake);

                filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_med_Vel_Zel_quart.txt')
%                 filename_wake=strcat('',specs_path_list_nowake,'/',string(sample_list_wake(sample)),type_folder,'check/vel/half/','_Check_med_Vel_Zel.txt')
                info_med = dlmread(filename_wake);
                
%                 mn(sample,:)=info(2,:);
                med(sample,:)=info_med(:,2);
%                 mod(sample,:)=info_mod(2,:);



end




a_factor=info_med(:,1);

% mn_mn=mean(mn);
% mn_std=std(mn);

med_mn=mean(med);
med_std=std(med);

% mod_mn=mean(mod);
% mod_std=std(mod);

cd('../../../../parameters')
    
for rds=1:length(a_factor)    
    [ ~, ~, vel_pert ] = wake( Gmu,-1+1/a_factor(rds));
    wake_vel_pert_zeld(rds,1)=vel_pert;
    
end

cd('../Analysis/wake_detection/consist_check/Zeldovich/')

% fig=figure('Visible', 'off');
fig=figure;
errorbar(a_factor,med_mn*10^-17,med_std*10^-17)
hold on
plot(a_factor,wake_vel_pert_zeld)

%xlim ([-inf inf]);
xlim ([0.08 0.26]);    %for paper
ylim ([2 4.5]*10^-19);    %for paper
xlabel('Scale factor', 'interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity (Mpc/h)/s', 'interpreter', 'latex', 'fontsize', 20);
set(gca,'FontName','FixedWidth');
set(gca,'FontSize',16);
set(gca,'linewidth',2);
title({strcat('Velocity comparizon')},'interpreter', 'latex', 'fontsize', 20);
% legend(strcat('G\mu = ',num2str(Gmu,'%.1E')),"Zel'dovich",'Location','northeast')
legend("Simulation","Zel'dovich",'Location','northeast')
hold off;



end