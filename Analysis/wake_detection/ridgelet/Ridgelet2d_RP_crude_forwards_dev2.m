function [ Radon_hor_h__,Radon_vert_v__] = Ridgelet2d_RP_crude_forwards_dev2(X)


% 2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)

% using real ifft fix, with padding



X_h=[zeros(size(X)); X;zeros(size(X))];
X_v=[zeros(size(X)), X,zeros(size(X))];

% % implicit variables

[ncx_h,ncy_h]=size(X_h);
[ncx_v,ncy_v]=size(X_v);

is_odd_ncx_h=mod(ncx_h,2);
is_odd_ncy_h=mod(ncy_h,2);

is_odd_ncx_v=mod(ncx_v,2);
is_odd_ncy_v=mod(ncy_v,2);



% 
% % implicit variable
% 
% [ncx,ncy]=size(X);
% 
% is_odd_ncx=mod(ncx,2);
% is_odd_ncy=mod(ncy,2);

% 
% %Apply FFT2
% 
% for i=1:ncx
%     for j=1:ncy
%         X_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X(i,j);
%     end
% end


%Apply FFT2  (% convention to the frequency center fft)

for i=1:ncx_h
    for j=1:ncy_h
        X_h_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_h/ncx_h)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_h/ncy_h)))*X_h(i,j);
    end
end


for i=1:ncx_v
    for j=1:ncy_v
        X_v_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_v/ncx_v)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_v/ncy_v)))*X_v(i,j);
    end
end



% F_ = fft(fft(X_).').';



F_h_ = fft(fft(X_h_).').';
F_v_ = fft(fft(X_v_).').';



% Augmentation if even


% 
% if is_odd_ncy==0 F_(:,end+1)=F_(:,1);end
% if is_odd_ncx==0 F_(end+1,:)=F_(1,:);end


if is_odd_ncy_h==0 F_h_(:,end+1)=F_h_(:,1);end
if is_odd_ncx_h==0 F_h_(end+1,:)=F_h_(1,:);end

if is_odd_ncy_v==0 F_v_(:,end+1)=F_v_(:,1);end
if is_odd_ncx_v==0 F_v_(end+1,:)=F_v_(1,:);end



% figure;imagesc(real(F_h_));colorbar;
% figure;imagesc(real(F_v_));colorbar;




%Apply Recto-Polar interpolation



% 
% %Vertical lines_centered
% 
% for t=1:ncx+(1-is_odd_ncx)
%     for r=1:ncy+(1-is_odd_ncy)
%         M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
%         y=min(max(round(r),1),ncy+(1-is_odd_ncy));
%         x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));
%         %         A((t-1)*ncy+r,1)=x;
%         %         A((t-1)*ncy+r,2)=y;
%         RPinterp_vert(t,r)=F_(x,y);
%     end
% end
% % 
% 


%Vertical lines_centered


for t=1:ncx_v+(1-is_odd_ncx_v)
    for r=1:ncy_v+(1-is_odd_ncy_v)
        M = (ncx_v+(1-is_odd_ncx_v)-(2*t)+1)/(ncy_v-1-2*(is_odd_ncy_v)+(1-is_odd_ncy_v));
        y = r;
        x = mod(-1+fix(t+M*(r-1-(is_odd_ncy_v))-ceil((ncx_v+1)/2))+ceil((ncx_v+1)/2),ncx_v+(1-is_odd_ncx_v))+1;    
%         A((t-1)*ncy+r,1)=x;
%         A((t-1)*ncy+r,2)=y;
        RPinterp_vert_v(t,r)=F_v_(x,y);
    end
end

% %Horizontal lines_centered
% 
% 
% for t=1:ncy+(1-is_odd_ncy)
%     for r=1:ncx+(1-is_odd_ncx)
%         M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
%         x=min(max(round(r),1),ncx+(1-is_odd_ncx));
%         y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
%         %         B((t-1)*ncx+r,1)=r;
%         %         B((t-1)*ncx+r,2)=t+M*(r-1);
%         RPinterp_hor(t,r)=F_(x,y);
%     end
% % end

%Horizontal lines_centered

for t=1:ncy_h+(1-is_odd_ncy_h)
    for r=1:ncx_h+(1-is_odd_ncx_h)
        M=(ncy_h+(1-is_odd_ncy_h)-(2*t)+1)/(ncx_h-1-2*(is_odd_ncx_h)+(1-is_odd_ncx_h));
        x = r ;
        y = mod(-1+fix(t+M*(r-1-(is_odd_ncx_h))-ceil((ncy_h+1)/2))+ceil((ncy_h+1)/2),ncy_h+(1-is_odd_ncy_h))+1;
%         x=min(max(round(r),1),ncx_h+(1-is_odd_ncx_h));
%         y=min(max(fix(t+M*(r-1)-ceil((ncy_h+1)/2))+ceil((ncy_h+1)/2),1),ncy_h+(1-is_odd_ncy_h));
%         B_((t-1)*ncx_h+(1-is_odd_ncx_h)+r,1)=r;
%         B_((t-1)*ncx_h+(1-is_odd_ncx_h)+r,2)=mod(fix(t+M*(r-1-(is_odd_ncx_h))-ceil((ncy_h+1)/2))+ceil((ncy_h+1)/2),ncy_h+(1-is_odd_ncy_h));   
        B_((t-1)*ncx_h+(1-is_odd_ncx_h)+r,1)=x;
        B_((t-1)*ncx_h+(1-is_odd_ncx_h)+r,2)=y;   
        
        RPinterp_hor_h(t,r)=F_h_(x,y);
    end
end

% figure;imagesc(real(RPinterp_hor_h));colorbar;
% figure;imagesc(real(RPinterp_vert_v));colorbar;







% RPinterp_hor_=[sum(RPinterp_hor')',RPinterp_hor];
% RPinterp_vert_=[sum(RPinterp_vert')',RPinterp_vert];

% Radon_hor = (ifft(RPinterp_hor_.').');
% Radon_vert = (ifft(RPinterp_vert_.').');
% 


% if is_odd_ncx==0 RPinterp_hor(:,end)=[]; end
% if is_odd_ncy==0 RPinterp_vert(:,end)=[]; end

% if is_odd_ncx==0 
%     RPinterp_hor(:,end)=[];
%     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% end
% if is_odd_ncy==0 
%     RPinterp_vert(:,end)=[];
%     RPinterp_vert(:,1)=real(RPinterp_vert(:,1));
% end

% if is_odd_ncx_h==0 
%     RPinterp_hor_h(:,end)=[];
%     RPinterp_hor_h_toReal=ifft(RPinterp_hor_h(:,1));    
% %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
%     for t=1:ncy_h+(1-is_odd_ncy_h)
%         RPinterp_hor_h_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncy_h+(1-is_odd_ncy_h)))))*RPinterp_hor_h_toReal(t);
%     end
%     RPinterp_hor_h(:,1)=RPinterp_hor_h_toReal_;
% end
% 
% if is_odd_ncy_v==0
%     RPinterp_vert_v(:,end)=[];
%     RPinterp_vert_v_toReal=ifft(RPinterp_vert_v(:,1));
%     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
%     for t=1:ncx_v+(1-is_odd_ncx_v)
%         RPinterp_vert_v_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncx_v+(1-is_odd_ncx_v)))))*RPinterp_vert_v_toReal(t);
%     end
%     RPinterp_vert_v(:,1)=RPinterp_vert_v_toReal_;
% end




% Obtain Radon by inverse fft of the the RP interpolations

Radon_hor_h = (ifft(RPinterp_hor_h.').');
Radon_vert_v = (ifft(RPinterp_vert_v.').');





% convention to the frequency center fft


for r=1:ncx_h+(1-is_odd_ncx_h)
    for t=1:ncy_h+(1-is_odd_ncy_h)
        Radon_hor_h__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncx_h+(1-is_odd_ncx_h)))))*Radon_hor_h(t,r);
    end
end

for t=1:ncx_v+(1-is_odd_ncx_v)
    for r=1:ncy_v+(1-is_odd_ncy_v)
        Radon_vert_v__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncy_v+(1-is_odd_ncy_v)))))*Radon_vert_v(t,r);
    end
end



% Radon_hor__=Radon_hor;
% Radon_vert__=Radon_vert;
% 
% figure;imagesc(real(Radon_hor__)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% 
% figure;imagesc(real(Radon_vert__)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);




% Redundant augmentation in "odd" case (to make the dimension implicit)



if is_odd_ncy_h==1 Radon_hor_h__(end+1,:)=Radon_hor_h__(1,:);end

if is_odd_ncx_v==1 Radon_vert_v__(end+1,:)=Radon_vert_v__(1,:);end


end
