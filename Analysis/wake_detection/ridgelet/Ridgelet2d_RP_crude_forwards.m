function [ Radon_hor_, Radon_vert_] = Ridgelet2d_RP_crude_forwards(X)

%  2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)



% implicit variable

[ncx,ncy]=size(X);

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);


%Apply FFT2

for i=1:ncx
    for j=1:ncy
        X_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X(i,j);
    end
end

F_ = fft(fft(X_).').';

if is_odd_ncx==0 F_(:,end+1)=F_(:,1);end
if is_odd_ncy==0 F_(end+1,:)=F_(1,:);end



% figure; imagesc(abs(F_));colorbar;







%Apply Recto-Polar interpolation


%Horizontal lines_centered

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
        y=min(max(round(r),1),ncy+(1-is_odd_ncy));
        x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));    
%         A((t-1)*ncy+r,1)=x;
%         A((t-1)*ncy+r,2)=y;
        RPinterp_hor(t,r)=F_(x,y);
    end
end


for t=1:ncy+(1-is_odd_ncy)
    for r=1:ncx+(1-is_odd_ncx)
        M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
        x=min(max(round(r),1),ncx+(1-is_odd_ncx));
        y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
%         B((t-1)*ncx+r,1)=r;
%         B((t-1)*ncx+r,2)=t+M*(r-1);        
        RPinterp_vert(t,r)=F_(x,y);
    end
end


% 
% figure; imagesc(abs(RPinterp_hor));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure; imagesc(abs(RPinterp_vert));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


%Radon Transformation

RPinterp_hor(:,end)=[];
RPinterp_vert(:,end)=[];

Radon_hor = (ifft(RPinterp_hor.').');
Radon_vert = (ifft(RPinterp_vert.').');

for i=1:ncx+(1-is_odd_ncx)
    for j=1:ncy
        Radon_hor_(i,j)=(exp(-1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*Radon_hor(i,j);
    end
end

for i=1:ncx
    for j=1:ncy+(1-is_odd_ncy)
        Radon_vert_(j,i)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*Radon_vert(j,i);
    end
end

% 
% figure;imagesc(real(Radon_hor_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure;imagesc(real(Radon_vert_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


end
