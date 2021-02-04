function [ output_args ] = Ridgelet2d_PL( input_args )
%using PolarLab

addpath(genpath('/home/asus/Dropbox/Disrael/Work/Research/NBodyWake/production/Analysis/wake_detection/slices_curv_2d'));

clearvars;
%Construct the sample image

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



nc=1+2*64;
% nc=1+2*16;
% nc=1+2*8;
ncx=nc;ncy=1+(nc-1)/2;
% ncx=nc;ncy=nc;
X=zeros(ncx,ncy);
X(1+(ncx-1)/4:end-(ncx-1)/4,1+(ncy-1)/4:end-(ncy-1)/4)=1;
% X(1+3*(ncx-1)/8:end-3*(ncx-1)/8,1+3*(ncy-1)/8:end-3*(ncy-1)/8)=1;
% X=ones(ncx,ncy);
% 
% nc=129;
% ncx=nc;ncy=nc;
% X=ones(ncx,ncy);


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




end

