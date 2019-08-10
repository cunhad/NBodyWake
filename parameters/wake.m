function [ vSgammaS displacement vel_pert ] = wake( Gmu,z )
%wake computes the relevant wake quantities at redshift z it it was
%produces at zi (in this case the redshift of recombination)
%(example) [ vSgammaS displacement vel_pert] = wake( 1E-6 ,15)

%comoving Mpc, second and solar mass units (if not otherwise specified)

%   Detailed explanation goes here

[ h OmegaBM OmegaCDM OmegaM OmegaL clight zi t_0 Hzero tensor_tilt spectral_indice sigma8 T_cmb_t0 Scalar_amplitude ] = cosmology(  );

% vSgammaS=clight*(sqrt(3))/3; %speed of cosmic string times Lorentz factor in Mpc/h per second units*/

vSgammaS=clight*0.611; %speed of cosmic string times Lorentz factor in Mpc/h per second units (v=0. before it was 0.5)*/

displacement=((12*3.14)/5)*Gmu*t_0*vSgammaS*(sqrt(1+zi))/(1+z);    %displacement in comoving coordinates Mpc/h

vel_pert=((8.*3.14)/5.)*(Gmu)*vSgammaS*(sqrt(1+zi))*(sqrt(1+(z)));  %velocity perturbation in comoving coordinates (Mpc/h)/s 

end