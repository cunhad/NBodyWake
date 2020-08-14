function [ output_args ] = Ridgelet2d( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


nc=128;
ncx=nc;ncy=nc/2;

X=zeros(ncx,ncy);
% X(end/2,:)=1;
% X(end/2,end/2)=1;
% X(nc/4:3*nc/4,nc/4:3*nc/4)=1;
% X(nc/2,:)=1;
X(:,nc/4:1+nc/4)=1;
% X(nc/2:1+nc/2,:)=1;
% X=ones(nc);

figure; imagesc(X);


% XP=sum(D,2);

% figure; plot(XP);

% Y = fftn(X);
% Y = fftshift(X,2);
% Y2 = fftshift(X);
% 
%  Y = fft(fft(X).').'
% Y2 = fft(fft(X).').';
% figure; imagesc(abs(Y2))

% F = fftshift(fftshift(X).').';
% figure; imagesc(abs(F));
% 
% 
% F = fftshift(fft(fftshift(fft(X)).',2),1).';
% figure; imagesc(abs(F));


 F = fftshift(fft(fftshift(fft(X).',2)).',2);
%F=fftshift(fft(X),1);
figure; imagesc(abs(F));

% figure; imagesc(abs(F));
% figure; imagesc(real(F));
% figure; imagesc(imag(F));



%X(:,nc/4:1+nc/4)=1;
% Proj=ifft(Y);
% % figure; plot(Proj);
%  
% Proj=ifft(Y)';
% figure; plot(Proj);


theta = (0:1/4:2)*pi;
 rad = 0:floor(sqrt((ncx^2)+ncy^2));
 [T, R] = meshgrid(theta, rad);
 [X, Y] = pol2cart(T, R);
 Z = X.^2+Y.^2;
 
 
 theta = (0:1/4:2)*pi;
 rad = 0:floor(sqrt((ncx^2)+ncy^2));
 [T, R] = meshgrid(theta, rad);
 [X, Y] = pol2cart(T, R);
 Z = X.^2+Y.^2;
 
 
 
 





























% F = circshift(Y,[1+nc/2 1+nc/2]);
% F=Y;


% figure; imagesc(abs(F));
% figure; imagesc(real(F));
% figure; imagesc(imag(F));

% xi=[1:nc,nc*ones(1,nc)];
% xf=[nc:-1:1,ones(1,nc)];
% yi=[ones(1,nc),1:nc];
% yf=[nc*ones(1,nc),nc:-1:1];

xi=[ones(1,2*nc)];
xf=[nc*ones(1,nc),nc:-1:1];
yi=[ones(1,2*nc)];
yf=[1:nc,nc*ones(1,nc)];

radon=zeros(2*nc,nc);

for i=1:2*nc
[x y]=bresenham(xi(i),yi(i),xf(i),yf(i));
ind = sub2ind([nc nc],x,y);
% radon(i,:)=ifft(F(ind));
radon(i,:)=ifftn(F(ind));
% radon(i,:)=ifftshift(F(ind));

end

figure; imagesc(abs(radon));
figure; imagesc(real(radon));
figure; imagesc(imag(radon));



end

