function [  ] = wake_instertion(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

filename_in='/home/asus/Programs/Gadget-2.0.7/cosmology/old/snapshot_000';
path_out='/home/asus/Programs/Gadget-2.0.7/ICs/';
mkdir(path_out,'wake');
filename_out=strcat(path_out,'wake/','ics_cosmo_64p.dat');


Gmu=1E-5;
z_insert=31;

Mpc_to_km=3.086e+19;

cd('../../parameters')

[ h OmegaBM OmegaCDM OmegaM OmegaL clight zi t_0 Hzero tensor_tilt spectral_indice sigma8 T_cmb_t0 Scalar_amplitude ] = cosmology(  );

vSgammaS=clight*(sqrt(3))/3; %speed of cosmic string times Lorentz factor in Mpc per second units*/

displacement=((12*3.14)/5)*Gmu*t_0*vSgammaS*(sqrt(1+zi))/(1+z_insert);    %displacement in comoving coordinates

vel_pert=((8.*3.14)/5.)*(Gmu)*vSgammaS*(sqrt(1+zi))*(sqrt(1+(z_insert)));  %velocity perturbation in comoving coordinates



%%now we will read, modify and write the snapshot


%read header

fid = fopen(strcat(filename_in));
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



%write header

file_out = fopen(filename_out,'w');
fwrite(file_out,block_size_1_head,'uint','l');
fwrite(file_out,Npart,'uint','l');
fwrite(file_out,Massarr,'double','l');
fwrite(file_out,Time,'double','l');
fwrite(file_out,Redshift,'double','l');
fwrite(file_out,FlagSfr,'int','l');
fwrite(file_out,FlagFeedback,'int','l');
fwrite(file_out,Nall,'int','l');
fwrite(file_out,FlagCooling,'int','l');
fwrite(file_out,NumFiles,'int','l');
fwrite(file_out,BoxSize,'double','l');
fwrite(file_out,Omega0,'double','l');
fwrite(file_out,OmegaLambda,'double','l');
fwrite(file_out,HubbleParam,'double','l');
fwrite(file_out,FlagAge,'int','l');
fwrite(file_out,FlagMetals,'int','l');
fwrite(file_out,NallHW,'int','l');
fwrite(file_out,flag_entr_ics,'int','l');
fwrite(file_out,garbage,'int','l');
fwrite(file_out,block_size_1_tail,'int','l');

% read positions

block_size_2_head=fread(fid, 1, 'uint','l') ;
pos_0=fread(fid, [3 Npart(1+0)], 'single','l') ;
pos_1=fread(fid, [3 Npart(1+1)], 'single','l') ;
pos_2=fread(fid, [3 Npart(1+2)], 'single','l') ;
pos_3=fread(fid, [3 Npart(1+3)], 'single','l') ;
pos_4=fread(fid, [3 Npart(1+4)], 'single','l') ;
pos_5=fread(fid, [3 Npart(1+5)], 'single','l') ;
block_size_2_tail=fread(fid, 1, 'uint','l') ;


%displace towards the wake
        
        dist_to_wake3=pos_1(3,:)-BoxSize/2; %is the vector that points to the wake at Z=nc/2 plane
%         displacement_to_wake3=-sign(pos_1(3,:)-BoxSize/2)*displacement; %the particles will be displaced towards the wake
        pos_1(3,:)=pos_1(3,:)-sign(dist_to_wake3)*displacement;

% write positions

fwrite(file_out,block_size_2_head,'uint','l');
fwrite(file_out,pos_0,'single','l');
fwrite(file_out,pos_1,'single','l');
fwrite(file_out,pos_2,'single','l');
fwrite(file_out,pos_3,'single','l');
fwrite(file_out,pos_4,'single','l');
fwrite(file_out,pos_5,'single','l');
fwrite(file_out,block_size_2_tail,'uint','l');

clearvars pos_1

%read velocities

block_size_3_head=fread(fid, 1, 'uint','l') ;
vel_0=fread(fid, [3 Npart(1+0)], 'single','l') ;
vel_1=fread(fid, [3 Npart(1+1)], 'single','l') ;
vel_2=fread(fid, [3 Npart(1+2)], 'single','l') ;
vel_3=fread(fid, [3 Npart(1+3)], 'single','l') ;
vel_4=fread(fid, [3 Npart(1+4)], 'single','l') ;
vel_5=fread(fid, [3 Npart(1+5)], 'single','l') ;
block_size_3_tail=fread(fid, 1, 'uint','l') ;

 %give hte velocity kick
        
%         kick_to_wake3=-sign(dist_to_wake3(1,:))*vel_pert;
        vel_1(3,:)=vel_1(3,:)-sign(dist_to_wake3(1,:))*vel_pert*Mpc_to_km*(1/(1+z_insert));

%write velocities

fwrite(file_out,block_size_3_head,'uint','l');
fwrite(file_out,vel_0,'single','l');
fwrite(file_out,vel_1,'single','l');
fwrite(file_out,vel_2,'single','l');
fwrite(file_out,vel_3,'single','l');
fwrite(file_out,vel_4,'single','l');
fwrite(file_out,vel_5,'single','l');
fwrite(file_out,block_size_3_tail,'uint','l');

clearvars vel_1 dist_to_wake3

%read particle ids

block_size_4_head=fread(fid, 1, 'uint','l') ;
pid_0=fread(fid, [1 Npart(1+0)], 'uint','l') ;
pid_1=fread(fid, [1 Npart(1+1)], 'uint','l') ;
pid_2=fread(fid, [1 Npart(1+2)], 'uint','l') ;
pid_3=fread(fid, [1 Npart(1+3)], 'uint','l') ;
pid_4=fread(fid, [1 Npart(1+4)], 'uint','l') ;
pid_5=fread(fid, [1 Npart(1+5)], 'uint','l') ;
block_size_4_tail=fread(fid, 1, 'uint','l') ;

%write particle ids

fwrite(file_out,block_size_4_head,'uint','l');
fwrite(file_out,pid_0,'uint','l');
fwrite(file_out,pid_1,'uint','l');
fwrite(file_out,pid_2,'uint','l');
fwrite(file_out,pid_3,'uint','l');
fwrite(file_out,pid_4,'uint','l');
fwrite(file_out,pid_5,'uint','l');
fwrite(file_out,block_size_4_tail,'uint','l');


%from now on a more carefull test must be made, since only the pure dm case
%was tested

%read Variable particle masses

Nm=int8(dot(single(Massarr==0),single(Npart~=0)));

block_size_5_head=fread(fid, int8(Nm~=0), 'uint','l') ;
var_masses=fread(fid, [1 Nm], 'single','l') ;
block_size_5_tail=fread(fid, Nm, 'uint','l') ;

%write Variable particle masses

fwrite(file_out,block_size_5_head,'uint','l');
fwrite(file_out,var_masses,'single','l');
fwrite(file_out,block_size_5_tail,'uint','l');

%read internal energy gas

block_size_6_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
u=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_6_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%write internal energy gas

fwrite(file_out,block_size_6_head,'uint','l');
fwrite(file_out,u,'single','l');
fwrite(file_out,block_size_6_tail,'uint','l');

%read density gas

block_size_7_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
rho=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_7_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%write density gas

fwrite(file_out,block_size_7_head,'uint','l');
fwrite(file_out,rho,'single','l');
fwrite(file_out,block_size_7_tail,'uint','l');


%read smoothing length gas

block_size_8_head=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;
hsml=fread(fid, [1 Npart(1+0)], 'single','l') ;
block_size_8_tail=fread(fid, int8(Npart(1+0)~=0), 'uint','l') ;

%write smoothing length gas

fwrite(file_out,block_size_8_head,'uint','l');
fwrite(file_out,hsml,'single','l');
fwrite(file_out,block_size_8_tail,'uint','l');


%%%%%
%from now one there will be problems if we want to recover those following
%quantities (will not affect a sucesfull copy-paste) since its presence
%(maximum four blocks) is determined in compilation and not know here
%%%%%

%read gravitational potential

block_size_9_head=fread(fid, 1, 'uint','l') ;
pot=fread(fid, [1 block_size_9_head/4], 'single','l') ;
block_size_9_tail=fread(fid, 1, 'uint','l') ;

%write gravitational potential

fwrite(file_out,block_size_9_head,'uint','l');
fwrite(file_out,pot,'single','l');
fwrite(file_out,block_size_9_tail,'uint','l');

%read Accelerations:

block_size_10_head=fread(fid, 1, 'uint','l') ;
acc=fread(fid, [1 block_size_10_head/4], 'single','l') ;
block_size_10_tail=fread(fid, 1, 'uint','l') ;

%write Accelerations:

fwrite(file_out,block_size_10_head,'uint','l');
fwrite(file_out,acc,'single','l');
fwrite(file_out,block_size_10_tail,'uint','l');

%read Rate of entropy production

block_size_11_head=fread(fid, 1, 'uint','l') ;
dAdt=fread(fid, [1 block_size_11_head/4], 'single','l') ;
block_size_11_tail=fread(fid, 1, 'uint','l') ;

%write Rate of entropy production

fwrite(file_out,block_size_11_head,'uint','l');
fwrite(file_out,dAdt,'single','l');
fwrite(file_out,block_size_11_tail,'uint','l');


%read Timesteps of particles

block_size_12_head=fread(fid, 1, 'uint','l') ;
dt=fread(fid, [1 block_size_12_head/4], 'single','l') ;
block_size_12_tail=fread(fid, 1, 'uint','l') ;

%write Timesteps of particles

fwrite(file_out,block_size_12_head,'uint','l');
fwrite(file_out,dt,'single','l');
fwrite(file_out,block_size_12_tail,'uint','l');



fclose(fid);
fclose(file_out);


cd('../wake_insertion/gadget');

end

