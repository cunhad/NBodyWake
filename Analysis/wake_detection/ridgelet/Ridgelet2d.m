function [ output_args ] = Ridgelet2d( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nc=128;
ncx=nc;ncy=nc/2;

X=zeros(ncx,ncy);
X(:,nc/4)=1;
X=X+randn([ncx,ncy]);
% X=randn([ncx,ncy]);
figure; imagesc(X);


% y = 2*randn([1,256]);
% Y = fft(y,256);Pyy = Y.*conj(Y)/256;
% plot(Pyy(1:128))
% 
% y = 2*randn([1,256]);
% Y = fft(y,256);Pyy = Y.*conj(Y)/256;
% plot(Pyy(1:128))
% 
% X=1+randn([1,256]);
% F=fft(X,256); Pyy = F.*conj(F);
% figure;plot(Pyy);



% X=randn([1,256]);F = fft(fft(X,256).',256).';figure; imagesc(abs(F));
% X=randn([ncx,ncy]);F = fft(fft(X,256).',256).';figure; imagesc(abs(F));
% F = fft2(X);
F = fftshift(fft(fftshift(fft(X).',2)).',2);
figure; imagesc(abs(F));

% F = fftshift(fft(fftshift(fft(X).',256)).',256);
% figure; imagesc(abs(F));

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


 
 [ncx_F,ncy_F]=size(F);
nc_r=max(ncx_F,ncx_F);
% theta = (0:1/(2*nc_r):-(1/(2*nc_r))+pi);
theta = linspace(0,pi,1+2*nc_r);
theta(end) = [];
rad = -ceil(sqrt((ncx_F^2)+(ncy_F^2))):ceil(sqrt((ncx_F^2)+(ncy_F^2)));
[T, R] = meshgrid(theta, rad);
% [x, y] = pol2cart(T, R);
x=R.*cos(T);
y=R.*sin(T);
ncx_centre=floor((ncx_F+1)/2);
ncy_centre=floor((ncy_F+1)/2);
lin_ncx_centred=(0:ncx_F-1)-ncx_centre;
lin_ncy_centred=(0:ncy_F-1)-ncy_centre;
[x_q,y_q] = meshgrid(lin_ncx_centred,lin_ncy_centred);
% F_intp = griddata(x_q(:),y_q(:),F(:),x,y,'cubic');
F2=F';
F_intp = griddata(x_q(:),y_q(:),F2(:),x,y,'cubic');

figure;imagesc(abs(F_intp)); colorbar;

figure;surf(abs(F_intp)); colorbar;
figure;surf(x,y,abs(F_intp)); colorbar;

F_intp_nonan=F_intp;

F_intp_nonan(isnan(F_intp))=0;


Radon = ifft(F_intp_nonan').';
figure;imagesc(real(Radon)); colorbar;
figure;imagesc(abs(Radon)); colorbar;

Radon = ifftshift(ifft(F_intp_nonan).',2);figure;imagesc(abs(Radon)); colorbar;


end

