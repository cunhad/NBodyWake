function [ Radon_hor_h_,Radon_vert_v_] = Ridgelet3d_RP_crude_forwards(X)

%  2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)
% try zero-padding


% (A,B,C) refers to the padded versions of X,Y, and Z dimensions

% implicit variables

% Choose this if padding

% X_A=[zeros(size(X)); X;zeros(size(X))];
% X_B=[zeros(size(X)), X,zeros(size(X))];
% X_C=cat(3,zeros(size(X)), X,zeros(size(X)));

% Or this is not padding
% 
% X_A=X;
% X_B=X;
% X_C=X;


[ncx_A,ncy_A,ncz_A]=size(X_A);
[ncx_B,ncy_B,ncz_B]=size(X_B);
[ncx_C,ncy_C,ncz_C]=size(X_C);

is_odd_ncx_A=mod(ncx_A,2);
is_odd_ncy_A=mod(ncy_A,2);
is_odd_ncz_A=mod(ncz_A,2);

is_odd_ncx_B=mod(ncx_B,2);
is_odd_ncy_B=mod(ncy_B,2);
is_odd_ncz_B=mod(ncz_B,2);

is_odd_ncx_C=mod(ncx_C,2);
is_odd_ncy_C=mod(ncy_C,2);
is_odd_ncz_C=mod(ncz_C,2);

%Apply FFT2

for i=1:ncx_A
    for j=1:ncy_A
        for k=1:ncz_A
            X_A_(i,j,k)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_A/ncx_A)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_A/ncy_A)))*(exp(1i*pi*(k-1)*(1-is_odd_ncz_A/ncz_A)))*X_A(i,j,k);
        end
    end
end


for i=1:ncx_B
    for j=1:ncy_B
        for k=1:ncz_B
            X_B_(i,j,k)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_B/ncx_B)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_B/ncy_B)))*(exp(1i*pi*(k-1)*(1-is_odd_ncz_B/ncz_B)))*X_B(i,j,k);
        end
    end
end

for i=1:ncx_C
    for j=1:ncy_C
        for k=1:ncz_C
            X_C_(i,j,k)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_C/ncx_C)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_C/ncy_C)))*(exp(1i*pi*(k-1)*(1-is_odd_ncz_C/ncz_C)))*X_C(i,j,k);
        end
    end
end

F_A_ = fftn(X_A_);
F_B_ = fftn(X_B_);
F_C_ = fftn(X_C_);

if is_odd_ncz_A==0 F_A_(:,:,end+1)=F_A_(:,:,1);end
if is_odd_ncy_A==0 F_A_(:,end+1,:)=F_A_(:,1,:);end
if is_odd_ncx_A==0 F_A_(end+1,:,:)=F_A_(1,:,:);end

if is_odd_ncz_B==0 F_B_(:,:,end+1)=F_B_(:,:,1);end
if is_odd_ncy_B==0 F_B_(:,end+1,:)=F_B_(:,1,:);end
if is_odd_ncx_B==0 F_B_(end+1,:,:)=F_B_(1,:,:);end

if is_odd_ncz_C==0 F_C_(:,:,end+1)=F_C_(:,:,1);end
if is_odd_ncy_C==0 F_C_(:,end+1,:)=F_C_(:,1,:);end
if is_odd_ncx_C==0 F_C_(end+1,:,:)=F_C_(1,:,:);end


% figure; imagesc(abs(F_h_));colorbar;
% figure; imagesc(abs(F_v_));colorbar;







%Apply Recto-Polar interpolation


% X_A lines_centered
% 
% aux_count=1;

for r = 1:ncx_A+(1-is_odd_ncx_A)
    for ty = 1:ncy_A+(1-is_odd_ncy_A)
        for tz = 1:ncz_A+(1-is_odd_ncz_A)
            My=(ncy_A+(1-is_odd_ncy_A)-(2*ty)+1)/(ncx_A-1+(1-is_odd_ncx_A));
            Mz=(ncz_A+(1-is_odd_ncz_A)-(2*tz)+1)/(ncx_A-1+(1-is_odd_ncx_A));
            x=min(max(round(r),1),ncx_A+(1-is_odd_ncx_A));
            y=min(max(fix(ty+My*(r-1)-ceil((ncy_A+1)/2))+ceil((ncy_A+1)/2),1),ncy_A+(1-is_odd_ncy_A));
            z=min(max(fix(tz+Mz*(r-1)-ceil((ncz_A+1)/2))+ceil((ncz_A+1)/2),1),ncz_A+(1-is_odd_ncz_A));
%             A(aux_count,1)=r;
%             A(aux_count,2)=ty;
%             A(aux_count,3)=tz;
%             A(aux_count,4)=x;
%             A(aux_count,5)=y;
%             A(aux_count,6)=z;
%             A(aux_count,7)=F_A_(x,y,z);
            RPinterp_A(r,ty,tz)=F_A_(x,y,z);
%             aux_count=aux_count+1;
        end
    end
end

% X_B lines_centered

for tx = 1:ncx_B+(1-is_odd_ncx_B)
    for r = 1:ncy_B+(1-is_odd_ncy_B)
        for tz = 1:ncz_B+(1-is_odd_ncz_B)
            Mx=(ncx_B+(1-is_odd_ncx_B)-(2*tx)+1)/(ncy_B-1+(1-is_odd_ncy_B));
            Mz=(ncz_B+(1-is_odd_ncz_B)-(2*tz)+1)/(ncy_B-1+(1-is_odd_ncy_B));
            y=min(max(round(r),1),ncy_B+(1-is_odd_ncy_B));
            x=min(max(fix(tx+Mx*(r-1)-ceil((ncx_B+1)/2))+ceil((ncx_B+1)/2),1),ncx_B+(1-is_odd_ncx_B));
            z=min(max(fix(tz+Mz*(r-1)-ceil((ncz_B+1)/2))+ceil((ncz_B+1)/2),1),ncz_B+(1-is_odd_ncz_B));
            RPinterp_B(r,tx,tz)=F_B_(x,y,z);
        end
    end
end

% X_C lines_centered

for tx = 1:ncx_C+(1-is_odd_ncx_C)
    for ty = 1:ncy_C+(1-is_odd_ncy_C)
        for r = 1:ncz_C+(1-is_odd_ncz_C)
            Mx=(ncx_C+(1-is_odd_ncx_C)-(2*tx)+1)/(ncz_C-1+(1-is_odd_ncz_C));
            My=(ncy_C+(1-is_odd_ncy_C)-(2*ty)+1)/(ncz_C-1+(1-is_odd_ncz_C));
            z=min(max(round(r),1),ncz_C+(1-is_odd_ncz_C));
            x=min(max(fix(tx+Mx*(r-1)-ceil((ncx_C+1)/2))+ceil((ncx_C+1)/2),1),ncx_C+(1-is_odd_ncx_C));
            y=min(max(fix(ty+My*(r-1)-ceil((ncy_C+1)/2))+ceil((ncy_C+1)/2),1),ncy_C+(1-is_odd_ncy_C));
            RPinterp_C(r,tx,ty)=F_C_(x,y,z);
        end
    end
end



% figure; imagesc(abs(RPinterp_hor_h));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure; imagesc(abs(RPinterp_vert_v));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


%Radon Transformation




RPinterp_A(end,:,:)=[];
RPinterp_B(end,:,:)=[];
RPinterp_C(end,:,:)=[];

Radon_A = ifft(RPinterp_A);
Radon_B = ifft(RPinterp_B);
Radon_C = ifft(RPinterp_C);



for r=1:ncx_A
    for j=1:ncy_A+(1-is_odd_ncy_A)
        for k=1:ncz_A+(1-is_odd_ncz_A)
%             Radon_A_(r,j,k)=(exp(-1i*pi*(r-1)*((1-is_odd_ncx_A)/ncx_A)))*Radon_A(r,j,k);
            Radon_A_(r,j,k)=(exp(-1i*pi*(r-1)*(1-is_odd_ncx_A/ncx_A)))*Radon_A(r,j,k);
%             Radon_A_(r,j,k)=(exp(-1i*pi*(j-1)*((1-is_odd_ncy_A)/ncy_A)))*(exp(-1i*pi*(k-1)*((1-is_odd_ncz_A)/ncz_A)))*Radon_A(r,j,k);
        end
    end
end

for i=1:ncx_B+(1-is_odd_ncx_B)
    for j=1:ncy_B
        Radon_vert_v_(i,j)=(exp(-1i*pi*(j-1)*(1-is_odd_ncy_B/ncy_B)))*Radon_A(i,j);
    end
end


for i=1:ncx_B+(1-is_odd_ncx_B)
    for j=1:ncy_B
        Radon_vert_v_(i,j)=(exp(-1i*pi*(j-1)*(1-is_odd_ncy_B/ncy_B)))*Radon_vert_v(i,j);
    end
end

for i=1:ncx_A
    for j=1:ncy_A+(1-is_odd_ncy_A)
        Radon_hor_h_(j,i)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx_A/ncx_A)))*Radon_hor_h(j,i);
    end
end

% 
% figure;imagesc(real(Radon_hor_h_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure;imagesc(real(Radon_vert_v_)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


end
