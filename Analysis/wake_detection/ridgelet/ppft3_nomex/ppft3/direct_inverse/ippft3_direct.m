% function rim = ippftd3d(pp, EPS)
%
% This function inverts ppft transformation and restores the initial image.
%   
% Oren 11/12/2005
function rim = ippft3_direct(pp, EPS)    
    ppx = OnionPeeling3D(pp, EPS);
    rim = invDecimatedFreqs3d(ppx); 
