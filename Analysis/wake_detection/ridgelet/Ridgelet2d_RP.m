function [ output_args ] = Ridgelet2d_RP( input_args )
% 2d Ridgelet transformation with recto-polar interpolation, not cetered

clearvars;
%Construct the sample image

nc=128;
ncx=nc;ncy=nc/4;
% ncx=nc;ncy=nc;

% X=ones(ncx,ncy);
X=zeros(ncx,ncy);
X(ncx/4:3*ncx/4,ncy/4:3*ncy/4)=1;
% X(:,nc/8)=1;
% X(:,nc/2)=1;
% for i=1:ncx
% %     X(i,min(max(round(i/4),1),ncy))=1;
%     X(i,10+min(max(round(i/8),1),ncy))=1;
% 
% end
% X=X+0.5*randn([ncx,ncy]);
% X=randn([ncx,ncy]);

% 
% %odd image test
% 



% nc=129;
% ncx=nc;ncy=1+(nc-1)/2;
% X=zeros(ncx,ncy);
% X((ncx-1)/4:1+end-(ncx-1)/4,(ncy-1)/4:1+end-(ncy-1)/4)=1;
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
imagesc(theta,xp,R);colorbar;

%Apply FFT2

F = fftshift(fft(fftshift(fft(X).',2)).',2);
% F = fft(fft(X).').';
% F = fftshift(fft(fftshift(fft(X).',2).'),2);
% F = fftshift(fft((fft(X).').'));

%use zero=mode-at-the-center convention

for i=1:ncx
    for j=1:ncy
        F_(i,j)=((-1)^(i+j))*F(i,j);
    end
end



figure; imagesc(abs(F_));colorbar;


% %recenter
% F_ext=zeros(size(F_)+1);
% F_ext(1:end-1,1:end-1)=F_;
% F_ext(:,end)=[F_(1,1); F_(end:-1:1,1)];
% F_ext(end,:)=[F_(1,1) F_(1,end:-1:1)];
% 
% figure; imagesc(abs(F_ext));colorbar;
% 
% [X,Y] = meshgrid(1:(ncy+1),1:(ncx+1));
% [Xq,Yq] = meshgrid(1:(ncy),1:(ncx));
% Xq=Xq+0.5;
% Yq=Yq+0.5;
% Fq = interp2(X,Y,F_ext,Xq,Yq,'cubic');
% 
% figure; imagesc(abs(Fq));colorbar;
% 
% F_=Fq;



%Apply Recto-Polar interpolation

%Horizontal lines

for t=1:ncx
    for r=1:ncy
        M=(ncx-(2*t)+1)/(ncy-1);
        Y=min(max(round(r),1),ncy);
        X=min(max(round(t+M*(r-1)),1),ncx);
        RPinterp_hor(t,r)=F_(X,Y);
    end
end

% figure; imagesc(abs(RPinterp_hor));colorbar;


%Vertical lines

for t=1:ncy
    for r=1:ncx
        M=(ncy-(2*t)+1)/(ncx-1);
        X=min(max(round(r),1),ncx);
        Y=min(max(round(t+M*(r-1)),1),ncy);
        RPinterp_vert(t,r)=F_(X,Y);
    end
end

% figure; imagesc(abs(RPinterp_hor));colorbar;
% figure; imagesc(abs(RPinterp_vert));colorbar;


% figure; plot(abs(RPinterp_vert(1,:)))
% figure; plot(abs(ifft(RPinterp_vert(1,:))))




%Radon Transformation

% Radon = ifftshift(ifft(RPinterp_hor).',2);
% Radon = ifftshift(ifft(RPinterp_vert).',2);
Radon_hor = (ifft(RPinterp_hor.').');
Radon_vert = (ifft(RPinterp_vert.').');

figure;imagesc(abs(Radon_hor)); colorbar;
figure;imagesc(abs(Radon_vert)); colorbar;


% figure;subplot(2,1,1)
% imagesc(abs(Radon_hor));
% subplot(2,1,2)
% imagesc(abs(Radon_vert));

% image1=imagesc(abs(Radon_hor));
% image2=imagesc(abs(Radon_vert));
% imshowpair(abs(Radon_hor),abs(Radon_vert),'montage');

















%inverse Radon Transformation

RPinterp_hor_rec=(fft(Radon_hor.').');
RPinterp_vert_rec=(fft(Radon_vert.').');

F_rec=zeros(ncx,ncy);
F_rec_count=zeros(ncx,ncy);

%Horizontal lines

for t=1:ncx
    for r=1:ncy
        M=(ncx-(2*t)+1)/(ncy-1);
        Y=min(max(round(r),1),ncy);
        X=min(max(round(t+M*(r-1)),1),ncx);
        F_rec(X,Y)=F_rec(X,Y)+RPinterp_hor_rec(t,r);
        F_rec_count(X,Y)=F_rec_count(X,Y)+1;
    end
end

% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Vertical lines

for t=1:ncy
    for r=1:ncx
        M=(ncy-(2*t)+1)/(ncx-1);
        X=min(max(round(r),1),ncx);
        Y=min(max(round(t+M*(r-1)),1),ncy);
        F_rec(X,Y)=F_rec(X,Y)+RPinterp_vert_rec(t,r);
        F_rec_count(X,Y)=F_rec_count(X,Y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

% figure; imagesc(abs(F_rec_norm));colorbar;




%get rid of zero-mode-at-the-center convention

for i=1:ncx
    for j=1:ncy
        F_rec_norm(i,j)=((-1)^(i+j))*F_rec_norm(i,j);
    end
end

figure; imagesc(abs(F_rec_norm));colorbar;

%un re-center


% F_rec_norm_ext=zeros(size(F_rec_norm)+1);
% 
% F_rec_norm_ext(2:end,2:end)=F_rec_norm;
% F_rec_norm_ext(:,1)=[F_rec_norm(end,end); F_rec_norm(end:-1:1,end)];
% F_rec_norm_ext(1,:)=[F_rec_norm(end,end) F_rec_norm(end,end:-1:1)];
% 
% % figure; imagesc(abs(F_rec_norm_ext));colorbar;
% 
% [X,Y] = meshgrid(0:(ncy),0:(ncx));
% [Xq,Yq] = meshgrid(1:(ncy),1:(ncx));
% Xq=Xq-0.5;
% Yq=Yq-0.5;
% F_rec_normq = interp2(X,Y,F_rec_norm_ext,Xq,Yq,'cubic');
% figure; imagesc(abs(F_rec_normq));colorbar;
% F_rec_norm=F_rec_normq;

%un re-center (2nd try)
% 
% F_rec_norm_ext=zeros(size(F_rec_norm));
% F_rec_norm_ext(2:end,2:end)=F_rec_norm(1:end-1,1:end-1);
% F_rec_norm_ext(1,2:end)=F_rec_norm(1,1:end-1);
% F_rec_norm_ext(2:end,1)=F_rec_norm(1:end-1,1);
% F_rec_norm_ext(1,1)=F_rec_norm(end,end);
% F_rec_norm=F_rec_norm_ext;

%Apply iFFT2

X_rec = ifft(ifftshift(ifft(ifftshift(F_rec_norm,2).'),2).');



figure; imagesc(abs(X_rec));colorbar;

% figure; imagesc(real(X_rec));colorbar;





diff=RPinterp_hor_inv-RPinterp_hor;
max(abs(diff(:)))


end