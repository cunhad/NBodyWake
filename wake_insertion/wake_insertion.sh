z_insert=31.0
Gmu_insert="1E-5"

path_=/gs/project/smj-701-aa/disrael/cubep3m/simulations/test/
spec_=64Mpc_96c_48p_zi63_nowakes
aux_path_=/

path="'$path_'"
spec="'$spec_'"
aux_path="'$aux_path_'"


matlab -nosplash -nodesktop -r "wake_insertion($path,$spec,$aux_path,$z_insert,$Gmu_insert);  exit;"
