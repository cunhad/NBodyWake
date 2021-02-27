function [ ] = Script_Ridgelet2d_RP_crude()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

clearvars;
%Construct the sample image


nc=128;
% nc=8;
% ncx=nc;ncy=nc/2;
% ncx=89;ncy=5;
ncx=nc;ncy=nc;

% X=ones(ncx,ncy);
% X=zeros(ncx,ncy);
% X(1+ncx/4:3*ncx/4,1+ncy/4:3*ncy/4)=1;
% X(:,nc/8)=1;
% X(:,nc/2)=1;
% for i=1:ncx
% %     X(i,min(max(round(i/4),1),ncy))=1;
%     X(i,10+min(max(round(i/8),1),ncy))=1;
% 
% end
% X=X+0.5*randn([ncx,ncy]);
X=randn([ncx,ncy]);

figure; imagesc(X);colorbar;

% 
% 
% nc=1+2*64;
% % nc=1+2*16;
% % nc=1+2*8;
% % nc=1+2*4;
% ncx=nc;ncy=1+(nc-1)/2;
% % ncx=nc;ncy=nc;
% X=zeros(ncx,ncy);
% X(1+(ncx-1)/4:end-(ncx-1)/4,1+(ncy-1)/4:end-(ncy-1)/4)=1;
% % X(1+3*(ncx-1)/8:end-3*(ncx-1)/8,1+3*(ncy-1)/8:end-3*(ncy-1)/8)=1;
% % X=ones(ncx,ncy);
% % 
% % nc=129;
% % ncx=nc;ncy=nc;
% % X=ones(ncx,ncy);
% 

% figure; imagesc(X);colorbar;

% figure; histogram(X(:));

%build in matlab rad transf
figure;
% theta = 0:180;
% theta-linspace()
theta = 180/nc:180/nc:180;
[R,xp] = radon(X,theta);
% imagesc(theta,xp,R);colorbar;
% imagesc(R(1+(ncx-1)/4:end-(ncx-1)/4,1:(end-1)/2)); colorbar;
imagesc(R(-(ncx-1)/4+(end+1)/2:(ncx-1)/4+(end+1)/2,1:(end-1)/2)); colorbar;



%Radon Transformation

[ Radon_vert_,Radon_hor_] = Ridgelet2d_RP_crude_forwards(X);


figure;imagesc(real(Radon_vert_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

figure;imagesc(real(Radon_hor_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

figure;imagesc(imag(Radon_vert_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

figure;imagesc(imag(Radon_hor_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

% 
% figure; histogram(real(Radon_vert_(:)));
% figure; histogram(real(Radon_hor_(:)));


[ X_rec_phase ] = Ridgelet2d_RP_crude_backwards(Radon_vert_, Radon_hor_  );

figure; imagesc(real(X_rec_phase));colorbar;

% max(abs((X_rec_phase(:)-X(:))))


Radon_hor=zeros(size(Radon_hor_));
Radon_vert=zeros(size(Radon_vert_));

% Radon_hor(65,1)=1;
% Radon_hor(34,1)=1;
Radon_hor(1,15)=1;

% figure;imagesc(real(Radon_hor)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure;imagesc(real(Radon_vert)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


[ X_rec_phase ] = Ridgelet2d_RP_crude_backwards( Radon_hor, Radon_vert);

figure; imagesc(real(X_rec_phase));colorbar;


end

