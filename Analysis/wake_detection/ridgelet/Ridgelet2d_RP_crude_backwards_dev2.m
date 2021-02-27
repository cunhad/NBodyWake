function [ X_rec_phase ] = Ridgelet2d_RP_crude_backwards_dev2( Radon_hor_, Radon_vert_ )

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
 

% %Apply inverse Recto-Polar interpolation


 
 % figure; imagesc(abs(RPinterp_hor_rec));colorbar;
 % figure; imagesc(abs(RPinterp_vert_rec));colorbar;
 
%   
% F_rec=zeros(ncx,ncy);
% F_rec_count=zeros(ncx,ncy);
% 
% RPinterp_vert_rec_=(fft(Radon_vert_phase_rec.').');
% RPinterp_hor_rec_=(fft(Radon_hor_phase_rec.').');
% 
% 
% RPinterp_vert_rec_v=RPinterp_vert_rec_(:,2+(ncy_-is_odd_ncy)/2:end-((ncy_-is_odd_ncy)/2));
% RPinterp_hor_rec_h=RPinterp_hor_rec_(:,2+(ncx_-is_odd_ncx)/2:end-(ncx_-is_odd_ncx)/2);

% 
% % the following is for the not pure real ridgelet
% 
% for i=1:ncx_v_+(1-is_odd_ncx)
%     for j=1:ncy_v_
%         Radon_vert_phase_rec_v(i,j)=(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy_v_)))*Radon_vert_v_(i,j);
%     end
% end
% 
% for i=1:ncx_h_
%     for j=1:ncy_h_+(1-is_odd_ncy)
%         Radon_hor_phase_rec_h(j,i)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx_h_)))*Radon_hor_h_(j,i);
%     end
% end


% % max(abs(Radon_hor_phase_rec_h(:)-Radon_hor_h(:)))
% 
% RPinterp_vert_rec_v=(fft(Radon_vert_phase_rec_v.').');
% RPinterp_hor_rec_h=(fft(Radon_hor_phase_rec_h.').');

%ver isso aqui

% max(abs(RPinterp_hor_rec_h(:)-RPinterp_hor_h(:)))


% RPinterp_vert_rec2=RPinterp_vert_rec_v(:,1:padd:end);
% RPinterp_hor_rec2=RPinterp_hor_rec_h(:,1:padd:end);


% e isso aqui tambem

% RPinterp_vert_rec2=RPinterp_vert_rec_v(:,1+is_odd_ncy:padd:end);
% RPinterp_hor_rec2=RPinterp_hor_rec_h(:,1+is_odd_ncx:padd:end);

% RPinterp_vert_rec2(:,end+1)=conj(RPinterp_vert_rec2(:,1));
% RPinterp_hor_rec2(:,end+1)=conj(RPinterp_hor_rec2(:,1));


% figure; imagesc(abs(RPinterp_hor_rec2));colorbar;
% figure; imagesc(abs(RPinterp_vert_rec2));colorbar;

% RPinterp_vert_rec3=RPinterp_vert_rec_v(:,1:3:end);
% RPinterp_hor_rec3=RPinterp_hor_rec_h(:,1:3:end);



% RPinterp_vert_rec3=RPinterp_vert_rec_v;
% RPinterp_hor_rec3=RPinterp_hor_rec_h;
% 
% RPinterp_vert_rec3(:,end+1)=conj(RPinterp_vert_rec3(:,1));
% RPinterp_hor_rec3(:,end+1)=conj(RPinterp_hor_rec3(:,1));
% 
% figure; imagesc(abs(RPinterp_hor_rec3));colorbar;
% figure; imagesc(abs(RPinterp_vert_rec3));colorbar;



% 
% F_rec=zeros(ncx_v_+(1-is_odd_ncx),ncy_v_+(1-is_odd_ncy));
% F_rec_count=zeros(ncx_v_+(1-is_odd_ncx),ncy_v_+(1-is_odd_ncy));

% F_rec=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));
% F_rec_count=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));







%Apply inverse Recto-Polar interpolation

F_rec=zeros(ncx,ncy);
F_rec_count=zeros(ncx,ncy);


% RPinterp_hor_rec
% RPinterp_hor_hor

%Vertical lines

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
        y=min(max(round(r),1),ncy+(1-is_odd_ncy));
        x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));   
        F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;
    end
end
% 
% for t=1:ncx_v_+(1-is_odd_ncx)
%     for r=1:ncy_v_+(1-is_odd_ncy)
%         M=(ncx_v_+(1-is_odd_ncx)-(2*t)+1)/(ncy_v_-1+(1-is_odd_ncy));
%         y=min(max(round(r),1),ncy_v_+(1-is_odd_ncy));
%         x=min(max(fix(t+M*(r-1)-ceil((ncx_v_+1)/2))+ceil((ncx_v_+1)/2),1),ncx_v_+(1-is_odd_ncx));   
%         F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec3(t,r);
%         F_rec_count(x,y)=F_rec_count(x,y)+1;
%     end
% end


% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Horizontal lines

for t=1:ncy+(1-is_odd_ncy)
    for r=1:ncx+(1-is_odd_ncx)
        M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
        x=min(max(round(r),1),ncx+(1-is_odd_ncx));
        y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
        F_rec(x,y)=F_rec(x,y)+RPinterp_hor_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

% F_rec_norm(isnan(F_rec_norm))=0;

% figure; imagesc(abs(F_rec_norm));colorbar;

if is_odd_ncy==0 F_rec_norm(:,end)=[];end
if is_odd_ncx==0 F_rec_norm(end,:)=[];end

X_rec = ifft(ifft(F_rec_norm.').');

for i=1:ncx
    for j=1:ncy
        X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(-1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X_rec(i,j);
    end
end


% for i=1:ncx
%     for j=1:ncy
%         X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(i-1+ncx)*(1-is_odd_ncx/(3*ncx))))*(exp(-1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*(exp(1i*pi*(j-1+ncy)*(1-is_odd_ncy/(3*ncy))))*X_rec(i,j);
%     end
% end


% figure; imagesc(real(X_rec_phase));colorbar;

% figure; imagesc(imag(X_rec_phase));colorbar;

% X_rec_phase2=X_rec_phase(:,1+ncy:2*ncy);

% figure; imagesc(real(X_rec_phase2));colorbar;

end

