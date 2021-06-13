function [ X_rec_phase ] = Ridgelet2d_RP_crude_backwards_dev2( Radon_hor_h_, Radon_vert_v_ )

% Radon_hor_h_ = Radon_hor_;
% Radon_vert_v_ = Radon_vert_;

% Inverse 2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)

% using real ifft fix, with padding

 
% implicit variable

ncy_h_ = size(Radon_hor_h_(:,1));
ncy_h_ = ncy_h_(1)-1;

ncx_v_ = size(Radon_vert_v_(:,1));
ncx_v_ = ncx_v_(1)-1;


ncx_h_ = 3*ncx_v_;
ncy_v_ = 3*ncy_h_;


ncx_= ncx_v_;
ncy_= ncy_h_;


% 
% 
% ncv_=size(Radon_vert_v_(:,:));
% ncy_v_=ncv_(2);
% ncy_h_=ncy_v_/3;
% % ncx_v_=ncv_(1);
% 
% % if mod(ncy_(2)-1,3)==0
% %     ncy_=ncy_(2)-1;
% % else
% %     ncy_=ncy_(2);
% % end
% 
% nch_=size(Radon_hor_h_(:,:));
% ncx_h_=nch_(2);
% ncx_v_=ncx_h_/3;
% % 
% % if mod(ncx(2)-1,3)==0
% %     ncx=ncx(2)-1;
% % else
% %     ncx=ncx(2);
% % end
% 
% 
% 
% ncx_= ncx_h_/3;
% ncy_= ncy_v_/3;



% [ncx,ncy]=size(Radon_hor_);
% ncx=ncx-1;

is_odd_ncx_=mod(ncx_,2);
is_odd_ncy_=mod(ncy_,2);


is_odd_ncx_h_=mod(ncx_h_,2);
is_odd_ncy_h_=mod(ncy_h_,2);

is_odd_ncx_v_=mod(ncx_v_,2);
is_odd_ncy_v_=mod(ncy_v_,2);



% Remove redundant augmentation in "odd" case (done to makd the dimension implicit)



if is_odd_ncy_h_==1 Radon_hor_h_(end,:)=[];end
if is_odd_ncx_v_==1 Radon_vert_v_(end,:)=[];end



% convention to the frequency center fft



for r=1:ncx_h_+(1-is_odd_ncx_h_)
    for t=1:ncy_h_+(1-is_odd_ncy_h_)
        Radon_hor_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/(ncx_h_+(1-is_odd_ncx_h_)))))*Radon_hor_h_(t,r);
    end
end

for t=1:ncx_v_+(1-is_odd_ncx_v_)
    for r=1:ncy_v_+(1-is_odd_ncy_v_)
        Radon_vert_phase_rec(t,r)=(exp(1i*pi*(r-1)*(1-1/(ncy_v_+(1-is_odd_ncy_v_)))))*Radon_vert_v_(t,r);
    end
end

% Obtain Radon by inverse fft of the the RP interpolations


RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');
 

% 
% if is_odd_ncx==0 RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1)); end
% if is_odd_ncy==0 RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1)); end


% 
% if is_odd_ncx_==0
%     %     RPinterp_hor_rec(:,end)=[];
%     %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
%     RPinterp_hor_toReal_rec_=RPinterp_hor_rec(:,1);
%     for t=1:ncy_+(1-is_odd_ncy_)
%         RPinterp_hor_toReal_rec_(t)=(exp(1i*pi*(t-1)*(1-1/(ncy_+(1-is_odd_ncy_)))))*RPinterp_hor_toReal_rec_(t);
%     end
%     RPinterp_hor_rec(:,1)=fft(RPinterp_hor_toReal_rec_(:));
%     RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1));
%     %     RPinterp_hor_rec(:,1)=RPinterp_hor_toReal_rec_;
%     display('the following should be small if inverse of a image')
%     display(imag(RPinterp_hor_rec(1,1)))
%     RPinterp_hor_rec(1,1)=real(RPinterp_hor_rec(1,1));
%     
% end
% 
% if is_odd_ncy_==0
%     %     RPinterp_vert_rec(:,end)=[];
%     %     RPinterp_vert_toReal_rec=ifft(RPinterp_vert_rec(:,1));
%     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
%     RPinterp_vert_toReal_rec_=RPinterp_vert_rec(:,1);
%     for t=1:ncx_+(1-is_odd_ncx_)
%         RPinterp_vert_toReal_rec_(t)=(exp(1i*pi*(t-1)*(1-1/(ncx_+(1-is_odd_ncx_)))))*RPinterp_vert_toReal_rec_(t);
%     end
% %     for t=1:ncx+(1-is_odd_ncx)
% %         RPinterp_vert_toReal__rec(t)=(exp(-1i*pi*(t-1)*(1-is_odd_ncx/ncx)))*RPinterp_vert_toReal_rec(t);
% %     end
%      RPinterp_vert_rec(:,1)=fft(RPinterp_vert_toReal_rec_(:));
%     RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1));
% %     RPinterp_vert_rec(:,1)=RPinterp_vert_toReal_rec_;
%     display('the following should be small if inverse of a image')
%     display(imag(RPinterp_vert_rec(1,1)))
%     RPinterp_vert_rec(1,1)=real(RPinterp_vert_rec(1,1));
% end
% 

% 
% RPinterp_vert_rec(:,end+1)=conj(RPinterp_vert_rec(:,1));
%  RPinterp_hor_rec(:,end+1)=conj(RPinterp_hor_rec(:,1));
%  
 
% RPinterp_hor_rec(1,1)=real(RPinterp_hor_rec(1,1));

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




% Removes the extra degree of freedom introduced by the padding

if is_odd_ncx_
    RPinterp_hor_rec_=RPinterp_hor_rec(:,2:3:end);
else
    RPinterp_hor_rec_=RPinterp_hor_rec(:,1:3:end);
end

if is_odd_ncy_
    RPinterp_vert_rec_=RPinterp_vert_rec(:,2:3:end);
else
    RPinterp_vert_rec_=RPinterp_vert_rec(:,1:3:end);
end


%Apply inverse Recto-Polar interpolation

F_rec=zeros(ncx_+(1-is_odd_ncx_),ncy_+(1-is_odd_ncy_));
F_rec_count=zeros(ncx_+(1-is_odd_ncx_),ncy_+(1-is_odd_ncy_));


% RPinterp_hor_rec
% RPinterp_hor_hor

%Vertical lines

for t=1:ncx_+(1-is_odd_ncx_)
    for r=1:ncy_+(1-is_odd_ncy_)
        M=(ncx_+(1-is_odd_ncx_)-(2*t)+1)/(ncy_-1+(1-is_odd_ncy_));
        y=min(max(round(r),1),ncy_+(1-is_odd_ncy_));
        x=min(max(fix(t+M*(r-1)-ceil((ncx_+1)/2))+ceil((ncx_+1)/2),1),ncx_+(1-is_odd_ncx_));   
        F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec_(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;
    end
end



% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Horizontal lines

for t=1:ncy_+(1-is_odd_ncy_)
    for r=1:ncx_+(1-is_odd_ncx_)
        M=(ncy_+(1-is_odd_ncy_)-(2*t)+1)/(ncx_-1+(1-is_odd_ncx_));
        x=min(max(round(r),1),ncx_+(1-is_odd_ncx_));
        y=min(max(fix(t+M*(r-1)-ceil((ncy_+1)/2))+ceil((ncy_+1)/2),1),ncy_+(1-is_odd_ncy_));
        F_rec(x,y)=F_rec(x,y)+RPinterp_hor_rec_(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

% F_rec_norm(isnan(F_rec_norm))=0;

% figure; imagesc(abs(F_rec_norm));colorbar;






% remove augmentation if even

if is_odd_ncy_==0 F_rec_norm(:,end)=[];end
if is_odd_ncx_==0 F_rec_norm(end,:)=[];end


%Apply iFFT2 (% convention to the frequency center fft)


X_rec = ifft(ifft(F_rec_norm.').');

for i=1:ncx_
    for j=1:ncy_
        X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx_/ncx_)))*(exp(-1i*pi*(j-1)*(1-is_odd_ncy_/ncy_)))*X_rec(i,j);
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

