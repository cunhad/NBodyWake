z_insert=15.0
Gmu_insert="8E-7"

path_=/gs/project/smj-701-aa/disrael/cubep3m/simulations/
spec_=64Mpc_1024c_512p_zi63_nowakem
aux_path_=/sample0001/

path="'$path_'"
spec="'$spec_'"
aux_path="'$aux_path_'"


matlab -nosplash -nodesktop -r "wake_insertion($path,$spec,$aux_path,$z_insert,$Gmu_insert);  exit;"
