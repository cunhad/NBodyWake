function [ ] = comvel_evol_mem_fast_meanslices_par_all( root_data_in,spec,wake_type)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% (example) []=comvel_evol_mem_fast_meanslices_par_all( '/home/asus/Dropbox/extras/storage/graham/small_res/data/','64Mpc_96c_48p_zi255_wakeGmu5t10m7zi63m','')

type_folder=wake_type;

%aux_path,type_folder,'check/vel/half/'

specs_path_list_nowake=strcat(root_data_out,spec)
sample_list_wake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};


end

