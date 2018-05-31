export LD_LIBRARY_PATH=/usr/local/lib
mpif90 -cpp -g  -O3 -ffree-line-length-none  -DDIAG -DNGP -DPPINT -DOPENMP  -DMPI_TIME  -DLRCKCORR -DNGPH -DDISP_MESH -DPID_FLAG -fopenmp hist3d.f90 -I/usr/local/include -o hist3d -L/usr/local/lib -lfftw3f_mpi -lfftw3f -lfftw3f_threads  -lm -ldl
