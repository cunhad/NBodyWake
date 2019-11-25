% function fim=OnionPeeling3D(pp, EPS)
%
% This function performs the onion-peeling procedure on pp grid in 3d case.
%

function fim=OnionPeeling3D(pp, EPS);
	n=checkInput(pp);
    fim = zeros(n+1, n+1,n+1);


    for j = 1:n/2
        
        fim(j, :,:) = recover2d_( ...
            pp(1, j*3-2,:,:), ...
            fim(j, :,:), ...
            j-1,3,EPS);
        
        fim(n+2-j,:,:) = recover2d_( ...
            pp(1, 3*n+4-j*3,n+1:-1:1,n+1:-1:1), ...
            fim(n+2-j, :,:), ...
            j-1,3,EPS);
        
        fim(:, j,:) = recover2d_( ...
            pp(2, j*3-2,:,:), ...
            fim(:, j,:), ...
            j-1,3,EPS);
        
        fim(:,n+2-j,:) = recover2d_( ...
            pp(2, 3*n+4-j*3,n+1:-1:1,n+1:-1:1), ...
            fim(:,n+2-j,:), ...
            j-1,3,EPS);

        fim(:, :,j) = recover2d_( ...
            pp(3, j*3-2,:,:), ...
            fim(:, :,j), ...
            j-1,3,EPS);
        
        fim(:,:,n+2-j) = recover2d_( ...
            pp(3, 3*n+4-j*3,n+1:-1:1,n+1:-1:1), ...
            fim(:,:,n+2-j), ...
            j-1,3,EPS);     

    end
    fim(n/2+1, n/2+1, n/2+1) = pp(1, n/2*3 +1,n/2+1,n/2+1);
    
    
    
function n=checkInput(pp)
% Check that the array pp is properly structured. pp
% should be of size 3x(3n+1)x(n+1)x(n+1). If everything is ok with pp
% the function returns n.
%
s = size(pp);

% if (s(1)~=3)|(mod(s(2),2)~=1) | (mod(s(3),2)~=1) | (mod(s(4),2)~=1) | (length(s) ~=4)
%    error('pp must be of size 3x(3n+1)x(n+1)x(n+1)');
% end
% 
% 
n = [(s(2)-1)/3;s(3)-1; s(4)-1];
% if (n(1)~=n(2))|(n(1)~=n(3))
%    error('pp must be of size 3x(3n+1)x(n+1)x(n+1)');
% end

n=n(1);
    