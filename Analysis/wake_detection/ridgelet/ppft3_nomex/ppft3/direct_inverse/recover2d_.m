% function re=recover2d(fu,u,j,d1,c1,d2, c2,d,c, EPS,a)
%
% this function recovers partially evaluated 2d-image. The substantional work of the algorithm is
% the evaluation of the central part of the image by the given values on its borders
% and by the resampled image fu which covers the central part.
%
%fu (n+1)x(n+1) - the image on pp grid
%u  (n+1)x(n+1) - partially evaluated image on carthesian grid
%j - integer - count of rows filled on u
%d1,d2,d 1..n; c 1..n+1; cl 1..n/2; c2 1..n/2+1 - constants used in recovering
%

function re=recover2d_(fu,u,j,a,EPS);
fu = squeeze(fu);
u  = squeeze(u); 
n = size(u, 1)-1;
alpha = (n/2 -j)/(n/2); 
m = 3*n+1;

    x = a*[-(n/2-j):(n/2-j)]*(-2)*pi/m;
    y = a*[[-n/2:-n/2+j-1]   [-n/2:n/2]*alpha   [n/2-j+1:n/2]]*(-2)*pi/m;
    
    x1 = [-n/2:n/2]*(-2)*pi*a/m*alpha;
    y1 = [-n/2:n/2]*(-2)*pi*a/m;


if j~=0
    pp11 = u;
    ppaa = fu;

        comp1a = zeros(n+1, n+1);
        re1a = zeros(n+1-2*j, n+1);
        re11 = zeros(n+1-2*j, n+1-2*j);


 	for k = 1:j
         comp1a(k, :) = ToeplitzResamp(pp11(k,1:n+1),y1,x1,n, EPS);
     end
    
    for k = 1:j
     %   comp1a(n+2-k, :) = fastresampleI(pp11(n+2-k,1:n+1),y1,x1,n);
      comp1a(n+2-k, :) = ToeplitzResamp(pp11(n+2-k,1:n+1),y1,x1,n, EPS);
    end
    
	for k = 1:n+1
	%	re1a(:,k) = fastresampleI([comp1a(1:j,k); ppaa(:,k); comp1a(n-j+2:n+1,k)], y, x,n);
        re1a(:,k) = ToeplitzResamp([comp1a(1:j,k); ppaa(:,k); comp1a(n-j+2:n+1,k)], y, x,n, EPS);
    end
	
	
	for k = 1:n+1-2*j
   %     re11(k,:) = fastresampleI([pp11(k+j,1:j) re1a(k,:) pp11(k+j,n-j+2:n+1)], y, x,n).';
       re11(k,:) = ToeplitzResamp([pp11(k+j,1:j) re1a(k,:) pp11(k+j,n-j+2:n+1)], y, x,n, EPS).';
    end
    
	re=pp11;
	re(1+j:n+1-j,1+j:n+1-j)=re11;
else
    re = fu;
end 
