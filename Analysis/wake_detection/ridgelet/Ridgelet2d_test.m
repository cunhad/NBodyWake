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
% X(:,nc/4:1+nc/4)=1;
X(:,nc/4)=1;
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
 
 
 
%  lin_spc_z_aux=linspace(min(lin_spc_z),max(lin_spc_z),length(lin_spc_z)*upsample_factor);
%  lin_spc_r_aux=linspace(min(lin_spc_r),max(lin_spc_r),1+length(lin_spc_r)*upsample_factor);
%  [lin_spc_zq,lin_spc_rq] = meshgrid(lin_spc_z,lin_spc_r);
%  [lin_spc_z_auxq,lin_spc_r_auxq] = meshgrid(lin_spc_z_aux,lin_spc_r_aux);
%  Cq_cut_extended_aux=Cq_cut_extended';
%  Cq_cut_extended_aux_grid = griddata(lin_spc_rq(:)',lin_spc_zq(:)',Cq_cut_extended_aux(:)',lin_spc_r_auxq,lin_spc_z_auxq,'cubic');
 
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

% F_intp_real = griddata(x_q(:),y_q(:),real(F(:)),x,y,'cubic');
% F_intp_imag = griddata(x_q(:),y_q(:),imag(F(:)),x,y,'cubic');
% lin_F_ncx_centred=(1:);
% lin_F_ncy_centred=(0:ncy-1)-ncy_centre;

% figure;surf(abs(F_intp)); colorbar;
figure;imagesc(abs(F_intp)); colorbar;

figure;surf(abs(F_intp)); colorbar;
figure;surf(x,y,abs(F_intp)); colorbar;

F_intp_nonan=F_intp;

F_intp_nonan(isnan(F_intp))=0;


% Radon=ifft2(F_intp_nonan);
% figure;surf(real(Radon)); colorbar;


% Radon = ifftshift(ifft(ifftshift(ifft(F_intp_nonan).',2)).',2);


% Radon = ifft(ifftshift(ifft(ifftshift(F_intp_nonan).',2)).',2);
% Radon = ifft(ifftshift(ifft(ifftshift(F_intp_nonan),1).'),1).';

% figure;imagesc(abs(Radon)); colorbar;
% 
% 
% 
% 
% Radon = ifftshift(ifft(F_intp_nonan).',2);
% figure;imagesc(abs(Radon)); colorbar;


Radon = ifft(F_intp_nonan').';
figure;imagesc(real(Radon)); colorbar;



% figure;surf(x,y,abs(F_intp_real)); colorbar;
% figure;surf(x,y,abs(F_intp)); colorbar;
% figure;surf(R,T,abs(F_intp)); colorbar;
% figure;surf(x_q,y_q,abs(F_intp)); colorbar;



























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




nc=128;
X=zeros(nc);
% X(end/2,:)=1;
% X(end/2,end/2)=1;
% X(nc/4:3*nc/4,nc/4:3*nc/4)=1;
X(:,nc/2)=1;
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

