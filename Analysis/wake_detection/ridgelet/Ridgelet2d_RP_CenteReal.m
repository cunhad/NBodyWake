function [ output_args ] = Ridgelet2d_RP_CenteReal( input_args )
% 2d Ridgelet transformation with recto-polar interpolation, try to ceter
% by upsampling beforehand

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
% nc=1+2*4;
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
%Apply FFT2

% F = fftshift(fft(fftshift(fft(X).',2)).',2);
% use zero=mode-at-the-center convention
for i=1:ncx
    for j=1:ncy
%         X_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*(exp(-1i*pi*(j-1)*(1-1/ncy)))*X(i,j);
%         X_(i,j)=(exp(-1i*pi*(i-1)*(1+1/ncx)))*(exp(-1i*pi*(j-1)*(1+1/ncy)))*X(i,j);
        X_(i,j)=(exp(1i*pi*(i-1)*(1-1/ncx)))*(exp(1i*pi*(j-1)*(1-1/ncy)))*X(i,j);
    end
end

F = fft(fft(X_).').';
figure; imagesc(abs(F));colorbar;

% F = fftshift(fft(fftshift(fft(X).',2)).',2);
% figure; imagesc(abs(F));colorbar;

% F = fft(fft(X).').';
% F = fftshift(fft(fftshift(fft(X).',2).'),2);
% F = fftshift(fft((fft(X).').'));

% %use zero=mode-at-the-center convention
% 
% for i=1:ncx
%     for j=1:ncy
%         F_(i,j)=((-exp(1i*pi/ncx))^(i-1))*((-exp(1i*pi/ncy))^(j-1))*F(i,j);
%     end
% end

% figure; imagesc(abs(F_));colorbar;


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

F_=F;

% %Horizontal lines
% for t=1:ncx
%     for r=1:ncy
%         M=(ncx-(2*t)+1)/(ncy-1);
%         y=min(max(round(r),1),ncy);
%         x=min(max(round(t+M*(r-1)),1),ncx);
% %         Y=min(max(floor(r),1),ncy);
% %         X=min(max(floor(t+M*(r-1)),1),ncx);
% %         A((t-1)*ncy+r,1)=X;
% %         A((t-1)*ncy+r,2)=Y;
%         Yh_q(t,r)=r;
%         Xh_q(t,r)=t+M*(r-1);
% %         B((t-1)*ncy+r,2)=r;
% %         B((t-1)*ncy+r,1)=t+M*(r-1);
%         RPinterp_hor(t,r)=F_(x,y);
% %         RPinterp_hor(t,r)=((-1)^(1+t))*F_(X,Y);
%     end
% end
% 

%Horizontal lines_centered

for t=1:ncx
    for r=1:ncy
        M=(ncx-(2*t)+1)/(ncy-1);
        y=min(max(round(r),1),ncy);
%         x=min(max(round(t+M*(r-1)),1),ncx);    
        x=min(max(fix(t+M*(r-1)-(ncx+1)/2)+(ncx+1)/2,1),ncx);    
%         Y=min(max(floor(r),1),ncy);
%         X=min(max(floor(t+M*(r-1)),1),ncx);
        A((t-1)*ncy+r,1)=x;
        A((t-1)*ncy+r,2)=y;
        Yh_q(t,r)=r;
        Xh_q(t,r)=t+M*(r-1);
%         B((t-1)*ncy+r,2)=r;
%         B((t-1)*ncy+r,1)=t+M*(r-1);
        RPinterp_hor(t,r)=F_(x,y);
%         RPinterp_hor(t,r)=((-1)^(1+t))*F_(X,Y);
    end
end


% x=min(max(round(t+M*([1:ncy]-1-((ncy-1)/2))),1),ncx)

% figure; imagesc(abs(RPinterp_hor));colorbar;

% [grid_x,grid_y] = meshgrid(1:ncx,1:ncy);
% gr_tr_x=grid_x';
% gr_tr_y=grid_y';
% RPinterp_hor_q = griddata(gr_tr_x(:),gr_tr_y(:),F_(:),Xh_q,Yh_q,'linear');
% figure; imagesc(abs(RPinterp_hor_q));colorbar;


% figure; imagesc(abs(RPinterp_hor));colorbar;
% figure; imagesc(abs(RPinterp_hor));colorbar;

%Vertical lines


for t=1:ncy
    for r=1:ncx
        M=(ncy-(2*t)+1)/(ncx-1);
        x=min(max(round(r),1),ncx);
%         y=min(max(round(t+M*(r-1)),1),ncy);
        y=min(max(fix(t+M*(r-1)-(ncy+1)/2)+(ncy+1)/2,1),ncy);
%         A((t-1)*ncx+r,1)=x;
%         A((t-1)*ncx+r,2)=y;
%         X=min(max(floor(r),1),ncx);
%         Y=min(max(floor(t+M*(r-1)),1),ncy);
        Xv_q(t,r)=r;
        Yv_q(t,r)=t+M*(r-1);
        B((t-1)*ncx+r,1)=r;
        B((t-1)*ncx+r,2)=t+M*(r-1);        
        RPinterp_vert(t,r)=F_(x,y);
%         RPinterp_vert(t,r)=((-1)^(1+t))*F_(X,Y);
    end
end

% figure; imagesc(abs(RPinterp_vert));colorbar;

% [grid_x,grid_y] = meshgrid(1:ncx,1:ncy);
% gr_tr_x=grid_x';
% gr_tr_y=grid_y';
% 
% RPinterp_vert_q = griddata(gr_tr_x(:),gr_tr_y(:),F_(:),Xv_q,Yv_q,'linear');
% figure; imagesc(abs(RPinterp_vert_q));colorbar;


figure; imagesc(abs(RPinterp_hor));colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

figure; imagesc(abs(RPinterp_vert));colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


%Centered interpolation

[grid_x,grid_y] = meshgrid(1:ncx,1:ncy);
gr_tr_x=grid_x';
gr_tr_y=grid_y';

RPinterp_hor_q = griddata(gr_tr_x(:),gr_tr_y(:),F_(:),Xh_q,Yh_q,'linear');
figure; imagesc(abs(RPinterp_hor_q));colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

RPinterp_vert_q = griddata(gr_tr_x(:),gr_tr_y(:),F_(:),Xv_q,Yv_q,'linear');
figure; imagesc(abs(RPinterp_vert_q));colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


% figure; plot(abs(RPinterp_vert(1,:)))
% figure; plot(abs(ifft(RPinterp_vert(1,:))))




%Radon Transformation

% Radon_hor = ifftshift(ifft(RPinterp_hor).',2);
% Radon_vert = ifftshift(ifft(RPinterp_vert).',2);
Radon_hor = (ifft(RPinterp_hor.').');
Radon_vert = (ifft(RPinterp_vert.').');

for i=1:ncx
    for j=1:ncy
        Radon_hor_(i,j)=(exp(-1i*pi*(j-1)*(1-1/ncy)))*Radon_hor(i,j);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

for i=1:ncx
    for j=1:ncy
        Radon_vert_(j,i)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*Radon_vert(j,i);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

figure;imagesc(real(Radon_hor_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

figure;imagesc(real(Radon_vert_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);

% Radon_hor_q = ifftshift(ifft(RPinterp_hor_q).',2);
% Radon_vert_q = ifftshift(ifft(RPinterp_vert_q).',2);
% 
% for i=1:ncx
%     for j=1:ncy
% %         F_rec_norm(i,j)=((-1)^(i+j))*F_rec_norm(i,j);
%         RPinterp_hor_q_norm(i,j)=((-1)^(j+0))*RPinterp_hor_q(i,j);
%     end
% end
% Radon_hor_q = (ifft(RPinterp_hor_q_norm.').');
% figure;imagesc(abs(Radon_hor_q)); colorbar;

% for i=1:ncx
%     for j=1:ncy
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(j-1)*(1-1/ncy)))*RPinterp_hor_q(i,j);
% %         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
%     end
% end

Radon_hor_q = (ifft(RPinterp_hor_q.').');
% figure;imagesc(abs(Radon_hor_q)); colorbar;

Radon_vert_q = (ifft(RPinterp_vert_q.').');
% figure;imagesc(abs(Radon_vert_q)); colorbar;


for i=1:ncx
    for j=1:ncy
        Radon_hor_q_(i,j)=(exp(-1i*pi*(j-1)*(1-1/ncy)))*Radon_hor_q(i,j);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

for i=1:ncx
    for j=1:ncy
        Radon_vert_q_(j,i)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*Radon_vert_q(j,i);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

figure;imagesc(real(Radon_hor_q_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


figure;imagesc(real(Radon_vert_q_)); colorbar;
xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);



% figure;subplot(2,1,1)
% imagesc(abs(Radon_hor));
% subplot(2,1,2)
% imagesc(abs(Radon_vert));

% image1=imagesc(abs(Radon_hor));
% image2=imagesc(abs(Radon_vert));
% imshowpair(abs(Radon_hor),abs(Radon_vert),'montage');

















%inverse Radon Transformation


for i=1:ncx
    for j=1:ncy
        Radon_hor_phase_rec(i,j)=(exp(1i*pi*(j-1)*(1-1/ncy)))*Radon_hor_(i,j);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

for i=1:ncx
    for j=1:ncy
        Radon_vert_phase_rec(j,i)=(exp(1i*pi*(i-1)*(1-1/ncx)))*Radon_vert_(j,i);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end





RPinterp_hor_rec=(fft(Radon_hor_phase_rec.').');
RPinterp_vert_rec=(fft(Radon_vert_phase_rec.').');



figure; imagesc(abs(RPinterp_hor_rec));colorbar;
figure; imagesc(abs(RPinterp_vert_rec));colorbar;






for i=1:ncx
    for j=1:ncy
        Radon_hor_phase_rec_q(i,j)=(exp(1i*pi*(j-1)*(1-1/ncy)))*Radon_hor_q_(i,j);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end

for i=1:ncx
    for j=1:ncy
        Radon_vert_phase_rec_q(j,i)=(exp(1i*pi*(i-1)*(1-1/ncx)))*Radon_vert_q_(j,i);
%         RPinterp_hor_q_(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*RPinterp_hor_q(i,j);
    end
end



RPinterp_hor_rec_q=(fft(Radon_hor_phase_rec_q.').');
RPinterp_vert_rec_q=(fft(Radon_vert_phase_rec_q.').');


figure; imagesc(abs(RPinterp_hor_rec_q));colorbar;
figure; imagesc(abs(RPinterp_vert_rec_q));colorbar;






F_rec=zeros(ncx,ncy);
F_rec_count=zeros(ncx,ncy);

%Horizontal lines

for t=1:ncx
    for r=1:ncy
        M=(ncx-(2*t)+1)/(ncy-1);
        y=min(max(round(r),1),ncy);
%         x=min(max(round(t+M*(r-1)),1),ncx);
        x=min(max(fix(t+M*(r-1)-(ncx+1)/2)+(ncx+1)/2,1),ncx);
%         Y=min(max(floor(r),1),ncy);
%         X=min(max(floor(t+M*(r-1)),1),ncx);
        F_rec(x,y)=F_rec(x,y)+RPinterp_hor_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;
    end
end

% figure; imagesc(F_rec_count);colorbar;
% figure; imagesc(abs(F_rec));colorbar;


%Vertical lines

for t=1:ncy
    for r=1:ncx
        M=(ncy-(2*t)+1)/(ncx-1);
        x=min(max(round(r),1),ncx);
%         y=min(max(round(t+M*(r-1)),1),ncy);
        y=min(max(fix(t+M*(r-1)-(ncy+1)/2)+(ncy+1)/2,1),ncy);
%         X=min(max(floor(r),1),ncx);
%         Y=min(max(floor(t+M*(r-1)),1),ncy);
        F_rec(x,y)=F_rec(x,y)+RPinterp_vert_rec(t,r);
        F_rec_count(x,y)=F_rec_count(x,y)+1;        
    end
end

%normalization

F_rec_norm=F_rec./F_rec_count;

figure; imagesc(abs(F_rec_norm));colorbar;
% 
% X_rec = ifft((ifft((F_rec_norm_phase.')).'));
% F = fft(fft(X_).').';


X_rec = ifft(ifft(F_rec_norm.').');

% X_rec = ifft(ifft(F_rec_norm_phase).').';

for i=1:ncx
    for j=1:ncy
%         F_rec_norm(i,j)=((-1)^(i+j))*F_rec_norm(i,j);
        X_rec_phase(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*(exp(-1i*pi*(j-1)*(1-1/ncy)))*X_rec(i,j);
    end
end




figure; imagesc(real(X_rec_phase));colorbar;






RPinterp=[RPinterp_hor_rec_q; RPinterp_vert_rec_q'];
Xinterp=[Xh_q;Xv_q'];
Yinterp=[Yh_q;Yv_q'];

% F_rec_norm_q = griddata(Xinterp(:),Yinterp(:),RPinterp(:),grid_y,grid_x,'linear')';

F_rec_norm_q = griddata(Xinterp(:),Yinterp(:),RPinterp(:),grid_x,grid_y,'linear')';


% figure; imagesc(abs(F_rec_norm));colorbar;


%get rid of zero-mode-at-the-center convention
% 
% for i=1:ncx
%     for j=1:ncy
% %         F_rec_norm(i,j)=((-1)^(i+j))*F_rec_norm(i,j);
%         F_rec_norm_phase_q(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*(exp(-1i*pi*(j-1)*(1-1/ncy)))*F_rec_norm_q(i,j);
%     end
% end

figure; imagesc(abs(F_rec_norm_q));colorbar;

X_rec_q = ifft(ifft(F_rec_norm_q.').');

for i=1:ncx
    for j=1:ncy
%         F_rec_norm(i,j)=((-1)^(i+j))*F_rec_norm(i,j);
        X_rec_phase_q(i,j)=(exp(-1i*pi*(i-1)*(1-1/ncx)))*(exp(-1i*pi*(j-1)*(1-1/ncy)))*X_rec_q(i,j);
    end
end


figure; imagesc(abs(X_rec_phase_q));colorbar;




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
% 
% X_rec = ifft(ifftshift(ifft(ifftshift(F_rec_norm_phase,2).'),2).');
% 
% 
% 
% figure; imagesc(abs(X_rec));colorbar;
% 
% % figure; imagesc(real(X_rec));colorbar;
% 
% 
% 
% 
% 
% diff=RPinterp_hor_inv-RPinterp_hor;
% max(abs(diff(:)))
% 

end