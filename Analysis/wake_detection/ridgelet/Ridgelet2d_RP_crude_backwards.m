function [ X_rec_phase ] = Ridgelet2d_RP_crude_backwards( Radon_hor_, Radon_vert_ )

% Inverse 2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)


% implicit variable

[ncx,ncy]=size(Radon_hor_);


for i=1:ncx
    for j=1:ncy
        Radon_hor_phase_rec(i,j)=(exp(1i*pi*(j-1)*(1-1/ncy)))*Radon_hor_(i,j);
    end
end

for i=1:ncx
    for j=1:ncy
        Radon_vert_phase_rec(j,i)=(exp(1i*pi*(i-1)*(1-1/ncx)))*Radon_vert_(j,i);
    end
end


RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');


% figure; imagesc(abs(RPinterp_hor_rec));colorbar;
% figure; imagesc(abs(RPinterp_vert_rec));colorbar;




F_rec=zeros(ncx,ncy);
F_rec_count=zeros(ncx,ncy);

%Horizontal lines

for t=1:ncx
    for r=1:ncy
        M=(ncx-(2*t)+1)/(ncy-1);
        y=min(max(round(r),1),ncy);
        x=min(max(fix(t+M*(r-1)-(ncx+1)/2)+(ncx+1)/2,1),ncx);
        F_rec(x,y)=F_rec(x,y)+RPinterp_hor_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;
    end
end

% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Vertical lines

for t=1:ncy
    for r=1:ncx
        M=(ncy-(2*t)+1)/(ncx-1);
        x=min(max(round(r),1),ncx);
        y=min(max(fix(t+M*(r-1)-(ncy+1)/2)+(ncy+1)/2,1),ncy);
        F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

% figure; imagesc(abs(F_rec_norm));colorbar;

X_rec = ifft(ifft(F_rec_norm.').');

for i=1:ncx
    for j=1:ncy
        X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*(exp(-1i*pi*(j-1)*(1-1/ncy)))*X_rec(i,j);
    end
end


% figure; imagesc(real(X_rec_phase));colorbar;



end

