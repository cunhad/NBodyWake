z_insert=31.0
Gmu_insert="6E-6"

path_=/home/acer/Documents/storage/guillimin/test2/
spec_=64Mpc_256c_zi63_nowakes
aux_path_=/
percentage_analysed_=1.0

path="'$path_'"
spec="'$spec_'"
aux_path="'$aux_path_'"


matlab -nosplash -nodesktop -r "wake_insertion($path,$spec,$aux_path,$percentage_analysed_,$z_insert,$Gmu_insert);  exit;"
