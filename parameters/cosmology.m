    function [ h OmegaBM OmegaCDM OmegaM OmegaL clight zi t_0 Hzero tensor_tilt spectral_indice sigma8 T_cmb_t0 Scalar_amplitude ] = cosmology(  )
%UNTITLED3 Summary of this function goes here

% (example) [ h OmegaBM OmegaCDM OmegaM OmegaL clight zi t_0 Hzero tensor_tilt spectral_indice sigma8 T_cmb_t0 Scalar_amplitude ] = cosmology(  ) 
 
%   Detailed explanation goes here
%Mpc/h, second and solar mass units, matter-dominated period

h=0.7;
OmegaCDM=0.246;
OmegaL=0.7095;
clight=(9.6E-15)*h;  %speed of light in Mpc/h per second units
%zi=1000;   %redshift close to the time of recombination
Hzero=(3.2411E-18)*h; %hubble constant

OmegaBM=0.0445;
OmegaM=OmegaCDM+OmegaBM;
tensor_tilt=1;
spectral_indice = 0.96;
sigma8=0.8628;
T_cmb_t0=2.7255;
Scalar_amplitude=2.49*10^(-9);

t_0=(1/sqrt(OmegaM))*(2/(3*Hzero)); %age of the universe today
zi=-1+2.4*(10^(4))*OmegaM*(h^2);

end

