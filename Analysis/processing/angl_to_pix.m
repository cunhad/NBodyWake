function [ ipix1 ] = angl_to_pix( theta,phi ,nside)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% [ pix_ind ] = angl_to_pix( 1.5708,3.1447 ,256)

piover2 = 0.5*pi;
twopi=2.0*pi;
z0=2.0/3.0;

z = cos(theta);
za = abs(z);
if( phi >= twopi)
    phi = phi - twopi;
end
if (phi < 0.)
    phi = phi + twopi;
end

tt = phi / piover2;  % in [0,4)

nl2 = 2*nside;
nl4 = 4*nside;
ncap  = nl2*(nside-1);    % number of pixels in the north polar cap
npix  = 12*nside*nside;

if( za <= z0 )
    
    jp = floor(nside*(0.5 + tt - z*0.75)); %/*index of ascending edge line*/
    jm = floor(nside*(0.5 + tt + z*0.75)); %/*index of descending edge line*/
    
    ir = nside + 1 + jp - jm;%// ! in {1,2n+1} (ring number counted from z=2/3)
    kshift = 0;
    if (mod(ir,2)==0.)
        kshift = 1;%// ! kshift=1 if ir even, 0 otherwise
    end
    
    ip = floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;%// ! in {1,4n}
    if( ip>nl4 )
        ip = ip - nl4;
    end
    
    ipix1 = ncap + nl4*(ir-1) + ip ;
    
else
    
    tp = tt - floor(tt);	%MOD(tt,1.d0)
    tmp = sqrt( 3.*(1. - za) );
    
    jp = floor( nside * tp * tmp );	% increasing edge line index
    jm = floor( nside * (1. - tp) * tmp );  % decreasing edge line index
    
    ir = jp + jm + 1;	%ring number counted from the closest pole
    ip = floor( tt * ir ) + 1;	% in {1,4*ir}
    if( ip>4*ir )
        ip = ip - 4*ir;
    end
    
    ipix1 = 2*ir*(ir-1) + ip;
    if( z<=0. )
        ipix1 = npix - 2*ir*(ir+1) + ip;
    end
end
  
end

