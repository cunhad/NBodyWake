function [  ] = read_output(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


 filename='/home/asus/Programs/Gadget-2.0.7/cosmology/old/snapshot_000';
% filename='/home/asus/Programs/MG-PICOLA-PUBLIC-master/files/snapshot_000';
% filename='/home/asus/Programs/MG-PICOLA-PUBLIC-master/output/out_z0p000.0';



%header

fid = fopen(strcat(filename));
block_size_1_head=fread(fid, 1, 'uint','l') ;
Npart=fread(fid, [6 1], 'uint','l') ;
Massarr=fread(fid, [6 1], 'double','l') ;
Time=fread(fid,1, 'double','l') ;
Redshift=fread(fid,1, 'double','l') ;
FlagSfr=fread(fid, 1, 'int','l') ;
FlagFeedback=fread(fid, 1, 'int','l') ;
Nall=fread(fid, [6 1], 'int','l') ;
FlagCooling=fread(fid, 1, 'int','l') ;
NumFiles=fread(fid, 1, 'int','l') ;
BoxSize=fread(fid,1, 'double','l') ;
Omega0=fread(fid,1, 'double','l') ;
OmegaLambda=fread(fid,1, 'double','l') ;
HubbleParam=fread(fid,1, 'double','l') ;
FlagAge=fread(fid, 1, 'int','l') ;
FlagMetals=fread(fid, 1, 'int','l') ;
NallHW=fread(fid, [6 1], 'int','l') ;
flag_entr_ics=fread(fid, 1, 'int','l') ;
garbage=fread(fid, [15 1], 'int','l') ;
block_size_1_tail=fread(fid, 1, 'uint','l') ;

%positions

block_size_2_head=fread(fid, 1, 'uint','l') ;
pos_0=fread(fid, [3 Npart(1+0)], 'single','l') ;
pos_1=fread(fid, [3 Npart(1+1)], 'single','l') ;
pos_2=fread(fid, [3 Npart(1+2)], 'single','l') ;
pos_3=fread(fid, [3 Npart(1+3)], 'single','l') ;
pos_4=fread(fid, [3 Npart(1+4)], 'single','l') ;
pos_5=fread(fid, [3 Npart(1+5)], 'single','l') ;
block_size_2_tail=fread(fid, 1, 'uint','l') ;

%velocities

block_size_3_head=fread(fid, 1, 'uint','l') ;
vel_0=fread(fid, [3 Npart(1+0)], 'single','l') ;
vel_1=fread(fid, [3 Npart(1+1)], 'single','l') ;
vel_2=fread(fid, [3 Npart(1+2)], 'single','l') ;
vel_3=fread(fid, [3 Npart(1+3)], 'single','l') ;
vel_4=fread(fid, [3 Npart(1+4)], 'single','l') ;
vel_5=fread(fid, [3 Npart(1+5)], 'single','l') ;
block_size_3_tail=fread(fid, 1, 'uint','l') ;

%particle ids

block_size_4_head=fread(fid, 1, 'uint','l') ;
pid_0=fread(fid, [1 Npart(1+0)], 'uint','l') ;
pid_1=fread(fid, [1 Npart(1+1)], 'uint','l') ;
pid_2=fread(fid, [1 Npart(1+2)], 'uint','l') ;
pid_3=fread(fid, [1 Npart(1+3)], 'uint','l') ;
pid_4=fread(fid, [1 Npart(1+4)], 'uint','l') ;
pid_5=fread(fid, [1 Npart(1+5)], 'uint','l') ;
block_size_4_tail=fread(fid, 1, 'uint','l') ;

%from now on a more carefull test must be made, since only the pure dm case
%was tested

%Variable particle masses

Nm=int8(dot(single(Massarr==0),single(Npart~=0)));

block_size_5_head=fread(fid, int8(Nm~=0), 'uint','l') ;
var_masses=fread(fid, [1 Nm], 'single','l') ;
block_size_5_tail=fread(fid, Nm, 'uint','l') ;

%internal energy gas

block_size_6_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
u=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_6_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%density gas

block_size_7_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
rho=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_7_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%smoothing length gas

block_size_8_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
hsml=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_8_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%%%%%
%from now one there will be problems if we whant to recover those following
%quantities (will not affect a sucesfull copy-paste) since its presence
%(maximum four blocks) is determined in compilation and not know here
%%%%%

%gravitational potential

block_size_9_head=fread(fid, 1, 'uint','l') ;
pot=fread(fid, [1 block_size_9_head/4], 'single','l') ;
block_size_9_tail=fread(fid, 1, 'uint','l') ;

%Accelerations:

block_size_10_head=fread(fid, 1, 'uint','l') ;
acc=fread(fid, [1 block_size_10_head/4], 'single','l') ;
block_size_10_tail=fread(fid, 1, 'uint','l') ;

%Rate of entropy production

block_size_11_head=fread(fid, 1, 'uint','l') ;
dAdt=fread(fid, [1 block_size_11_head/4], 'single','l') ;
block_size_11_tail=fread(fid, 1, 'uint','l') ;

%Timesteps of particles

block_size_12_head=fread(fid, 1, 'uint','l') ;
dt=fread(fid, [1 block_size_12_head/4], 'single','l') ;
block_size_12_tail=fread(fid, 1, 'uint','l') ;



fclose(fid);


end

