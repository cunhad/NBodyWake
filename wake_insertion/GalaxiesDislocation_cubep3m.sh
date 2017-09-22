initial_redshift=10.0
log10_Gmu=-6
path_in=/home/acer/Documents/storage/16Mpc_128c_zi63_nowakes/
path_out=/home/acer/Documents/storage/16Mpc_128c_zi63_wakeGmu1t10m6zi10s/
path_nowake=/home/acer/Documents/storage/16Mpc_128c_zi63_nowakes/
my_path=/home/acer/Dropbox/Disrael/Doutorado/Research/NBodyWake/Analysis/from_terminal/
file_pos_name_in=00xv0.dat
file_pos_name_nowake=00xv_nowake.dat
file_pos_name_out=00xv0_wake.dat
box_size=16.0
n_cells=128.0

cp $path_in$initial_redshift$file_pos_name_in $path_nowake$initial_redshift$file_pos_name_nowake

g++ GalaxiesDislocation_cubep3m.cpp -o InsertWake

./InsertWake $initial_redshift $log10_Gmu $path_in $path_out $file_pos_name_in $file_pos_name_out $box_size $n_cells

mv $path_nowake$initial_redshift$file_pos_name_nowake  $path_nowake$initial_redshift$file_pos_name_in
mv $path_out$initial_redshift$file_pos_name_out $path_out$initial_redshift$file_pos_name_in

cd $path_in
#cp pk* seed* delta* $path_out 
