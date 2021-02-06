function [ output_args ] = Ridgelet2d_PL( input_args )
%using PolarLab

addpath(genpath('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/Analysis/wake_detection/slices_curv_2d'));

clearvars;
%Construct the sample image

N=64;
nc=64;
ncx=nc;ncy=nc;

% nc=128;
% ncx=nc;ncy=nc/4;
% % ncx=nc;ncy=nc;
% 
% % X=ones(ncx,ncy);
% X=zeros(ncx,ncy);
% X(ncx/4:3*ncx/4,ncy/4:3*ncy/4)=1;
% % X(:,nc/8)=1;
% % X(:,nc/2)=1;
% % for i=1:ncx
% % %     X(i,min(max(round(i/4),1),ncy))=1;
% %     X(i,10+min(max(round(i/8),1),ncy))=1;
% % 
% % end
% % X=X+0.5*randn([ncx,ncy]);
% % X=randn([ncx,ncy]);
% 
% % 
% % %odd image test
% % 


% 
% nc=1+2*64;
% % nc=1+2*16;
% % nc=1+2*8;
% ncx=nc;ncy=1+(nc-1)/2;
% % ncx=nc;ncy=nc;
X=zeros(ncx,ncy);
X(1+(ncx-1)/4:end-(ncx-1)/4,1+(ncy-1)/4:end-(ncy-1)/4)=1;
% % X(1+3*(ncx-1)/8:end-3*(ncx-1)/8,1+3*(ncy-1)/8:end-3*(ncy-1)/8)=1;
% % X=ones(ncx,ncy);
% % 
% % nc=129;
% % ncx=nc;ncy=nc;
% % X=ones(ncx,ncy);

% X=zeros(N); X(N/2-2:N/2+3,N/2-3:N-1)=1;
 

figure; imagesc(X);colorbar;

%build in matlab rad transf
figure;
% theta = 0:180;
% theta-linspace()
theta = 180/nc:180/nc:180;
[R,xp] = radon(X,theta);
% imagesc(theta,xp,R);colorbar;
% imagesc(R(1+(ncx-1)/4:end-(ncx-1)/4,1:(end-1)/2)); colorbar;
imagesc(R(-(ncx-1)/4+(end+1)/2:(ncx-1)/4+(end+1)/2,1:(end-1)/2)); colorbar;


% figure; imagesc(X); colormap(gray(256)); axis image; colorbar;

% The ground truth computed brute-force
[Xc,Yc]=Create_Grid('X',[N,pi],'');
% Yref=Brute_Force_Transform(X,Xc,Yc); 
OS=20;
% Y=AFTUSF_Spline_1(X,Xc,Yc,OS); 
Y=XPolar_Transform(X,3,8,1,OS); 
F_=Y;

nc=2*nc;
ncx=nc;ncy=nc;

figure; imagesc(abs(F_));colorbar;

% 
% 
% %Apply Recto-Polar interpolation
% 
% %Horizontal lines
% 
% for t=1:ncx
%     for r=1:ncy
%         M=(ncx-(2*t)+1)/(ncy-1);
%         Y=min(max(round(r),1),ncy);
%         X=min(max(round(t+M*(r-1)),1),ncx);
%         RPinterp_hor(t,r)=F_(X,Y);
%     end
% end
% 
% % figure; imagesc(abs(RPinterp_hor));colorbar;
% 
% 
% %Vertical lines
% 
% for t=1:ncy
%     for r=1:ncx
%         M=(ncy-(2*t)+1)/(ncx-1);
%         X=min(max(round(r),1),ncx);
%         Y=min(max(round(t+M*(r-1)),1),ncy);
%         RPinterp_vert(t,r)=F_(X,Y);
%     end
% end
% 
% % figure; imagesc(abs(RPinterp_hor));colorbar;
% % figure; imagesc(abs(RPinterp_vert));colorbar;
% 
% 
% % figure; plot(abs(RPinterp_vert(1,:)))
% % figure; plot(abs(ifft(RPinterp_vert(1,:))))
% 
% 


%Radon Transformation

% Radon = ifftshift(ifft(RPinterp_hor).',2);
% Radon = ifftshift(ifft(RPinterp_vert).',2);
Radon_hor = (ifft(F_.').');
figure;imagesc(abs(Radon_hor)); colorbar;



end

