function [ X_rec_phase ] = Ridgelet2d_RP_crude_backwards( Radon_hor_, Radon_vert_ )

% Inverse 2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)

% implicit variables

ncy=size(Radon_hor_(1,:));
ncy=ncy(2);

ncx=size(Radon_vert_(1,:));
ncx=ncx(2);

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);



for i=1:ncx+(1-is_odd_ncx)
    for j=1:ncy
        Radon_hor_phase_rec(i,j)=(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*Radon_hor_(i,j);
    end
end

for i=1:ncx
    for j=1:ncy+(1-is_odd_ncy)
        Radon_vert_phase_rec(j,i)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*Radon_vert_(j,i);
    end
end


RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');

% Parei aqui

RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1));
RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1));

% figure; imagesc(abs(RPinterp_hor_rec));colorbar;
% figure; imagesc(abs(RPinterp_vert_rec));colorbar;




F_rec=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));
F_rec_count=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));

%Horizontal lines

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
        y=min(max(round(r),1),ncy+(1-is_odd_ncy));
        x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));   
        F_rec(x,y)=F_rec(x,y)+RPinterp_hor_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;
    end
end

% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Vertical lines

for t=1:ncy+(1-is_odd_ncy)
    for r=1:ncx+(1-is_odd_ncx)
        M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
        x=min(max(round(r),1),ncx+(1-is_odd_ncx));
        y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
        F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

% figure; imagesc(abs(F_rec_norm));colorbar;

if is_odd_ncx==0 F_rec_norm(:,end)=[];end
if is_odd_ncy==0 F_rec_norm(end,:)=[];end

X_rec = ifft(ifft(F_rec_norm.').');

for i=1:ncx
    for j=1:ncy
        X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(-1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X_rec(i,j);
    end
end


% figure; imagesc(real(X_rec_phase));colorbar;



end

