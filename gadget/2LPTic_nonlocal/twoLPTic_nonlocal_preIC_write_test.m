function [  ] = twoLPTic_nonlocal_preIC_read(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filename='/home/asus/Programs/N-GenIC/dummy_glass_original.dat';
path_out='/home/asus/Programs/2LPTic_nonlocal/inputs/';
% mkdir(path_out,'wake');
filename_out=strcat(path_out,'test_out_32');

file_out = fopen(filename_out,'w');

bytesleft=120;

N  = 64;
Ntot = N * N * N;
BoxSize = 1200;


npart=zeros(6,1);
fwrite(file_out,npart,'long','l');
massarr=zeros(6,1);
fwrite(file_out, massarr, 'double','l') ;
time=0;
fwrite(file_out, time, 'double','l') ;
redshift=0;
fwrite(file_out, redshift, 'double','l') ;
flag_sfr=0;
fwrite(file_out,flag_sfr, 'long','l') ;
flag_feedback=0;
fwrite(file_out,flag_feedback, 'long','l') ;
npartall=zeros(1,6);
fwrite(file_out, npartall, 'long','l') ;
flag_cooling= 0;
fwrite(file_out,flag_cooling, 'long','l') ;
num_files= 1;
fwrite(file_out,num_files, 'long','l') ;
% BoxSize = 320;
fwrite(file_out,BoxSize, 'double','l') ;
la =zeros(1,bytesleft/2);
fwrite(file_out,la, 'short','l') ;
% fread(fid, 2, 'int32','l') ;
%pos= fread(fid, [3 Ntot], 'single','l') ;

for i= 0: N-1 
  for j=0: N-1
    for k=0:N-1
      pos(1, (i*N+j)*N+k+1) = (i+0.0)/N * BoxSize;
      pos(2, (i*N+j)*N+k+1) = (j+0.0)/N * BoxSize;
      pos(3, (i*N+j)*N+k+1) = (k+0.0)/N * BoxSize;
    end
  end
end

fwrite(file_out,pos, 'single','l') ;

fclose(fid);

end

