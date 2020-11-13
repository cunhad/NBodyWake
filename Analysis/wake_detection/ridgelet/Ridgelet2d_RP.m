function [ output_args ] = Ridgelet2d_RP( input_args )
% 2d Ridgelet transformation with recto-polar interpolation

clearvars;
%Construct the sample image

nc=128;
% ncx=nc;ncy=nc/4;
ncx=nc;ncy=nc;


X=zeros(ncx,ncy);
% X(:,nc/8)=1;
X(:,nc/2)=1;
% for i=1:ncx
% %     X(i,min(max(round(i/4),1),ncy))=1;
%     X(i,min(max(round(i),1),ncy))=1;
% 
% end
% X=X+0.5*randn([ncx,ncy]);
% X=randn([ncx,ncy]);
figure; imagesc(X);colorbar;

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






end