function [ output_args ] = Ridgelet2d( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


nc=128;
X=zeros(nc);
% X(end/2,:)=1;
% X(end/2,end/2)=1;
X(nc/4:3*nc/4,nc/4:3*nc/4)=1;
% X=ones(nc);


% figure; imagesc(X);

Y = fftn(X);
% Y = fftshift(X);

% figure; imagesc(abs(Y));
% figure; imagesc(real(Y));
% figure; imagesc(imag(Y));

% F = circshift(Y,[1+nc/2 1+nc/2]);
F=Y;


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

