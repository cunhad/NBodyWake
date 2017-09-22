initial_redshift=20.0
log10_Gmu=-6
path_in=/scratch/irulan/camargod/25Mpc_128c_zi20_nowakes/
path_out=/scratch/irulan/camargod/25Mpc_128c_zi20_wakeGmu3t10m6zi20s/
path_nowake=/scratch/irulan/camargod/25Mpc_128c_zi20_nowakes/
my_path=/homes/jivaro/camargod/Programs/WakeVizualization/
file_pos_name_in=00xv0.dat
file_pos_name_nowake=00xv_nowake.dat
file_pos_name_out=00xv0_wake.dat
box_size=25.0
n_cells=128.0

cp $path_in$initial_redshift$file_pos_name_in $path_nowake$initial_redshift$file_pos_name_nowake

g++ GalaxiesDislocation_cubep3m_irulan.cpp -o InsertWake

#echo $initial_redshift $log10_Gmu $path_in $path_out $file_pos_name_in $file_pos_name_out $box_size $n_cells

./InsertWake $initial_redshift $log10_Gmu $path_in $path_out $file_pos_name_in $file_pos_name_out $box_size $n_cells

mv $path_nowake$initial_redshift$file_pos_name_nowake  $path_nowake$initial_redshift$file_pos_name_in
mv $path_out$initial_redshift$file_pos_name_out $path_out$initial_redshift$file_pos_name_in

cd $path_in
#cp pk* seed* delta* $path_out 
