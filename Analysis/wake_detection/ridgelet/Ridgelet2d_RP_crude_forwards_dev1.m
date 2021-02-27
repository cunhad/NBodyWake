function [ Radon_hor_h__,Radon_vert_v__] = Ridgelet2d_RP_crude_forwards_dev1(X)

%  2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)
% try zero-padding


% implicit variable

% Choose this if padding

X_h=[zeros(size(X)); X;zeros(size(X))];
X_v=[zeros(size(X)), X,zeros(size(X))];

% Or this is not padding

% X_h=X;
% X_v=X;


[ncx_h,ncy_h]=size(X_h);
[ncx_v,ncy_v]=size(X_v);

is_odd_ncx_h=mod(ncx_h,2);
is_odd_ncy_h=mod(ncy_h,2);

is_odd_ncx_v=mod(ncx_v,2);
is_odd_ncy_v=mod(ncy_v,2);


%Apply FFT2

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



F_h_ = fft(fft(X_h_).').';
F_v_ = fft(fft(X_v_).').';

if is_odd_ncy_h==0 F_h_(:,end+1)=F_h_(:,1);end
if is_odd_ncx_h==0 F_h_(end+1,:)=F_h_(1,:);end

if is_odd_ncy_v==0 F_v_(:,end+1)=F_v_(:,1);end
if is_odd_ncx_v==0 F_v_(end+1,:)=F_v_(1,:);end


% figure; imagesc(abs(F_h_));colorbar;
% figure; imagesc(abs(F_v_));colorbar;


%see this

% if is_odd_ncx_h==1 F_h_(1,:)=[];F_h_(end,:)=[];ncx_h=ncx_h-2; end
% 
% if is_odd_ncy_v==1 F_v_(:,1)=[];F_v_(:,end)=[];ncy_v=ncy_v-2; end






%Apply Recto-Polar interpolation


%Vertical lines_centered


for t=1:ncx_v+(1-is_odd_ncx_v)
    for r=1:ncy_v+(1-is_odd_ncy_v)
        M=(ncx_v+(1-is_odd_ncx_v)-(2*t)+1)/(ncy_v-1+(1-is_odd_ncy_v));
        y=min(max(round(r),1),ncy_v+(1-is_odd_ncy_v));
        x=min(max(fix(t+M*(r-1)-ceil((ncx_v+1)/2))+ceil((ncx_v+1)/2),1),ncx_v+(1-is_odd_ncx_v));    
%         A((t-1)*ncy+r,1)=x;
%         A((t-1)*ncy+r,2)=y;
        RPinterp_vert_v(t,r)=F_v_(x,y);
    end
end

%Horizontal lines_centered

for t=1:ncy_h+(1-is_odd_ncy_h)
    for r=1:ncx_h+(1-is_odd_ncx_h)
        M=(ncy_h+(1-is_odd_ncy_h)-(2*t)+1)/(ncx_h-1+(1-is_odd_ncx_h));
        x=min(max(round(r),1),ncx_h+(1-is_odd_ncx_h));
        y=min(max(fix(t+M*(r-1)-ceil((ncy_h+1)/2))+ceil((ncy_h+1)/2),1),ncy_h+(1-is_odd_ncy_h));
%         B((t-1)*ncx+r,1)=r;
%         B((t-1)*ncx+r,2)=t+M*(r-1);        
        RPinterp_hor_h(t,r)=F_h_(x,y);
    end
end



% figure; imagesc(abs(RPinterp_hor_h));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure; imagesc(abs(RPinterp_vert_v));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


%Radon Transformation

%gives a pure real result

RPinterp_hor_h_a=cat(2,RPinterp_hor_h(:,1+(ncx_h-is_odd_ncx_h)/2:end),RPinterp_hor_h,RPinterp_hor_h(:,1:(ncx_h-is_odd_ncx_h)/2));
Radon_hor_h_b = (ifft(RPinterp_hor_h_a.').');
Radon_hor_h__=Radon_hor_h_b(:,1:2:end);
% Radon_hor_h__(:,end)=[];

RPinterp_vert_v_a=cat(2,RPinterp_vert_v(:,1+(ncy_v-is_odd_ncy_v)/2:end),RPinterp_vert_v,RPinterp_vert_v(:,1:(ncy_v-is_odd_ncy_v)/2));
Radon_vert_v_b = (ifft(RPinterp_vert_v_a.').');
Radon_vert_v__=Radon_vert_v_b(:,1:2:end);
% Radon_vert_v__(:,end)=[];

% figure;imagesc(real(Radon_hor_h__));colorbar;
% figure;imagesc(real(Radon_vert_v__));colorbar;



% the following does not gives a pure real result

% RPinterp_vert_v(:,end)=[];
% RPinterp_hor_h(:,end)=[];
% 
% 
% 
% Radon_vert_v = (ifft(RPinterp_vert_v.').');
% Radon_hor_h = (ifft(RPinterp_hor_h.').');
% 
% for i=1:ncx_v+(1-is_odd_ncx_v)
%     for j=1:ncy_v
%         Radon_vert_v_(i,j)=(exp(-1i*pi*(j-1)*(1-is_odd_ncy_v/ncy_v)))*Radon_vert_v(i,j);
%     end
% end
% 
% for i=1:ncx_h
%     for j=1:ncy_h+(1-is_odd_ncy_h)
%         Radon_hor_h_(j,i)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx_h/ncx_h)))*Radon_hor_h(j,i);
%     end
% end

% 
% 
% figure;imagesc(real(Radon_hor_h_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure;imagesc(real(Radon_vert_v_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


end

