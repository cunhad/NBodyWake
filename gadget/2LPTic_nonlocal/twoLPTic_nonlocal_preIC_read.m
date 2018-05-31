function [  ] = twoLPTic_nonlocal_preIC_read(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% filename='/home/asus/Programs/N-GenIC/dummy_glass_original.dat';
filename='/home/asus/Programs/2LPTic_nonlocal/inputs/glass1_le';
% path_out='/home/asus/Programs/2LPTic_nonlocal/';
% filename=strcat('/home/asus/Programs/2LPTic_nonlocal/test_out_32');

fid = fopen(strcat(filename));

bytesleft=120

N  = 280;
Ntot = N * N * N;

block_size_1_head=fread(fid, 1, 'uint','l') ;
npart=fread(fid, [6 1], 'long','l') ;
massarr=fread(fid, [6 1], 'double','l') ;
time=fread(fid, 1, 'double','l') ;
redshift=fread(fid, 1, 'double','l') ;
flag_sfr=fread(fid, 1, 'long','l') ;
flag_feedback=fread(fid, 1, 'long','l') ;
npartall=fread(fid, [6 1], 'long','l') ;
flag_cooling= fread(fid, 1, 'long','l') ;
num_files= fread(fid, 1, 'long','l') ;
BoxSize = fread(fid, 1, 'double','l') ;
la = fread(fid, [1 bytesleft/2], 'short','l') ;
block_size_1_tail=fread(fid, 1, 'uint','l') ;

% fread(fid, 2, 'int32','l') ;
%pos= fread(fid, [3 Ntot], 'single','l') ;

block_size_2_head=fread(fid, 1, 'uint','l') ;

pos_0=fread(fid, [3 npart(1+0)], 'single','l') ;
pos_1=fread(fid, [3 npart(1+1)], 'single','l') ;
pos_2=fread(fid, [3 npart(1+2)], 'single','l') ;
pos_3=fread(fid, [3 npart(1+3)], 'single','l') ;
pos_4=fread(fid, [3 npart(1+4)], 'single','l') ;
pos_5=fread(fid, [3 npart(1+5)], 'single','l') ;
block_size_2_tail=fread(fid, 1, 'uint','l') ;

% pos= fread(fid, [3 Inf], 'single','l') ;

fclose(fid);

end

