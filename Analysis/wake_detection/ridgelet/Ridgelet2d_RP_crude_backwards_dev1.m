function [ X_rec_phase ] = Ridgelet2d_RP_crude_backwards_dev1( Radon_hor_, Radon_vert_)

% Inverse 2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)

% using real ifft fix
% 
% RPinterp_hor_phase_=Radon_hor_(:,end);
% RPinterp_vert_phase_=Radon_vert_(:,end);
% 
% Radon_hor_(:,end)=[];
% Radon_vert_(:,end)=[];


% 
% 
% for r=1:length(Radon_hor_(1,:))
%     for t=1:length(Radon_hor_(:,1))
%         Radon_hor_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/length(Radon_hor_(1,:)))))*Radon_hor_(t,r);
%     end
% end
% 
% for t=1:length(Radon_vert_(:,1))
%     for r=1:length(Radon_vert_(1,:))
%         Radon_vert_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/length(Radon_vert_(1,:)))))*Radon_vert_(t,r);
%     end
% end
% 
% 
% 
% RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
% RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');
% 
% 
% if sum(RPinterp_vert_rec(:,1)==conj(RPinterp_vert_rec(end:-1:1,1)))/length(RPinterp_vert_rec(:,1))==1
%     ncy_ = length(RPinterp_vert_rec(1,:));
% end
% 
% if  sum(RPinterp_vert_rec(2:end,1)==conj(RPinterp_vert_rec(end:-1:2,1)))/length(RPinterp_vert_rec(:,1))==0
%     ncy_ = length(RPinterp_vert_rec(1,:));
% end
% 
% if sum(RPinterp_hor_rec(2:end,1)==conj(RPinterp_hor_rec(end:-1:2,1)))/length(RPinterp_hor_rec(:,1))==1
%     ncx_ = length(RPinterp_hor_rec(1,:))-1;
% else
%     ncx_ = length(RPinterp_hor_rec(1,:));
% end
% 
% 
% 
% % implicit variable

ncy=size(Radon_hor_(:,1));
ncy = ncy(1)-1;

ncx=size(Radon_vert_(:,1));
ncx = ncx(1)-1;






% 
% % implicit variable
%  
% ncy=size(Radon_vert_(1,:));
% ncy=ncy(2);
% 
% % if mod(ncy_(2)-1,3)==0
% %     ncy_=ncy_(2)-1;
% % else
% %     ncy_=ncy_(2);
% % end
% 
% ncx=size(Radon_hor_(1,:));
% ncx=ncx(2);
% % 
% % if mod(ncx(2)-1,3)==0
% %     ncx=ncx(2)-1;
% % else
% %     ncx=ncx(2);
% % end
% 
% 



% [ncx,ncy]=size(Radon_hor_);
% ncx=ncx-1;

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);



if is_odd_ncy==1 Radon_hor_(end,:)=[];end
if is_odd_ncx==1 Radon_vert_(end,:)=[];end



for r=1:ncx+(1-is_odd_ncx)
    for t=1:ncy+(1-is_odd_ncy)
        Radon_hor_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/(ncx+(1-is_odd_ncx)))))*Radon_hor_(t,r);
    end
end

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        Radon_vert_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/(ncy+(1-is_odd_ncy)))))*Radon_vert_(t,r);
    end
end





RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');
 

% 
% if is_odd_ncx==0 RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1)); end
% if is_odd_ncy==0 RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1)); end


% 
% if is_odd_ncx==0 
% %     RPinterp_hor(:,end)=[];
% %     RPinterp_hor_imag(end)=[];
% %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_hor(:,1)=abs(RPinterp_hor(:,1)).*sign(abs(RPinterp_hor(:,1)));
%     RPinterp_hor_rec(:,1)=abs(RPinterp_hor_rec(:,1)).*exp(+1i*RPinterp_hor_phase_);
%     RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1));
% end
%     
% 
% if is_odd_ncy==0 
% %     RPinterp_hor(:,end)=[];
% %     RPinterp_hor_imag(end)=[];
% %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_hor(:,1)=abs(RPinterp_hor(:,1)).*sign(abs(RPinterp_hor(:,1)));
%     RPinterp_vert_rec(:,1)=abs(RPinterp_vert_rec(:,1)).*exp(+1i*RPinterp_vert_phase_);
%     RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1));
% end
%     


% 
% 
% 
% % RPinterp_hor_phase=angle(RPinterp_hor(:,end));
% 
% if is_odd_ncx==0 
%     RPinterp_hor(:,end)=[];
% %     RPinterp_hor_imag(end)=[];
% %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
%     RPinterp_hor(:,1)=abs(RPinterp_hor(:,1)).*sign(abs(RPinterp_hor(:,1)));
% end
%     
% 
% 
% if is_odd_ncx==0
%     %     RPinterp_hor_rec(:,end)=[];
%     %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
%     RPinterp_hor_toReal_rec=RPinterp_hor_rec(:,1);
%     for t=1:ncy+(1-is_odd_ncy)
%         RPinterp_hor_toReal_rec_(t)=(exp(1i*pi*(t-1)*(1-1/(ncy+(1-is_odd_ncy)))))*RPinterp_hor_toReal_rec(t);
%     end
%     RPinterp_hor_rec(:,1)=fft(RPinterp_hor_toReal_rec_(:));
%     RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1));
%     %     RPinterp_hor_rec(:,1)=RPinterp_hor_toReal_rec_;
%     
% end
% 
% if is_odd_ncy==0
%     %     RPinterp_vert_rec(:,end)=[];
%     %     RPinterp_vert_toReal_rec=ifft(RPinterp_vert_rec(:,1));
%     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
%     RPinterp_vert_toReal_rec=RPinterp_vert_rec(:,1);
%     for t=1:ncx+(1-is_odd_ncx)
%         RPinterp_vert_toReal_rec_(t)=(exp(1i*pi*(t-1)*(1-1/(ncx+(1-is_odd_ncx)))))*RPinterp_vert_toReal_rec(t);
%     end
% %     for t=1:ncx+(1-is_odd_ncx)
% %         RPinterp_vert_toReal__rec(t)=(exp(-1i*pi*(t-1)*(1-is_odd_ncx/ncx)))*RPinterp_vert_toReal_rec(t);
% %     end
%      RPinterp_vert_rec(:,1)=fft(RPinterp_vert_toReal_rec_(:));
%     RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1));
% %     RPinterp_vert_rec(:,1)=RPinterp_vert_toReal_rec_;
% end
% 
% 

% 
% 
% if is_odd_ncx==0 
%     RPinterp_hor(:,end)=[];
%     RPinterp_hor_toReal=ifft(RPinterp_vert(:,1));    
% %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
%     for t=1:ncy+(1-is_odd_ncy)
%         RPinterp_hor_toReal_(t)=(exp(-1i*pi*(t-1)*(1-is_odd_ncy/ncy)))*RPinterp_hor_toReal(t);
%     end
%     RPinterp_hor(:,1)=RPinterp_hor_toReal_;
% end
% 
% if is_odd_ncy==0
%     RPinterp_vert(:,end)=[];
%     RPinterp_vert_toReal=ifft(RPinterp_vert(:,1));
%     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
%     for t=1:ncx+(1-is_odd_ncx)
%         RPinterp_vert_toReal_(t)=(exp(-1i*pi*(t-1)*(1-is_odd_ncx/ncx)))*RPinterp_vert_toReal(t);
%     end
%     RPinterp_vert(:,1)=RPinterp_vert_toReal_;
% end
% % 
% if is_odd_ncx==0 RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1)); end
% if is_odd_ncy==0 RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1)); end


%Apply inverse Recto-Polar interpolation

F_rec=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));
F_rec_count=zeros(ncx+(1-is_odd_ncx),ncy+(1-is_odd_ncy));


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
