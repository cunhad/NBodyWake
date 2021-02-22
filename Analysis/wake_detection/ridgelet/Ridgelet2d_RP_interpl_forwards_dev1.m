function [ Radon_hor_h_,Radon_vert_v_] = Ridgelet2d_RP_interpl_forwards_dev1(X)

%  2-D Radon transformation, recto-polar crude interpolation (nearest, close to center)
% try zero-padding


% implicit variable

% Choose this if padding

% X_h=[zeros(size(X)); X;zeros(size(X))];
% X_v=[zeros(size(X)), X,zeros(size(X))];

% Or this is not padding

X_h=X;
X_v=X;


[ncx_h,ncy_h]=size(X_h);
[ncx_v,ncy_v]=size(X_v);

is_odd_ncx_h=mod(ncx_h,2);
is_odd_ncy_h=mod(ncy_h,2);

is_odd_ncx_v=mod(ncx_v,2);
is_odd_ncy_v=mod(ncy_v,2);


%Apply FFT2

for i=1:ncx_h
    for j=1:ncy_h
        X_h_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_h/ncx_h)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_h/ncy_h)))*X_h(i,j);
    end
end


for i=1:ncx_v
    for j=1:ncy_v
        X_v_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx_v/ncx_v)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy_v/ncy_v)))*X_v(i,j);
    end
end



F_h_ = fft(fft(X_h_).').';
F_v_ = fft(fft(X_v_).').';

% clearvars F_h_A F_h_B F_h_C

% if is_odd_ncx_h==0 F_h_A(:,end+1)=F_h_A(:,1);end
% if is_odd_ncy_h==0 F_h_A(end+1,:)=F_h_A(1,:);end

if is_odd_ncx_h==0    
    F_h_(:,end+1)=F_h_(:,1);
    F_h_(:,end+1)=F_h_(:,2);
    F_h_(:,end+1)=F_h_(:,3);
    F_h_(:,end+1)=F_h_(:,end-4);
    F_h_(:,end+1)=F_h_(:,end-4);
    F_h_ = circshift(F_h_,[0 2]);    
else    
    F_h_(:,end+1)=F_h_(:,1);
    F_h_(:,end+1)=F_h_(:,2);
    F_h_(:,end+1)=F_h_(:,end-3);
    F_h_(:,end+1)=F_h_(:,end-3);
    F_h_ = circshift(F_h_,[0 2]);
    
end

if is_odd_ncy_h==0   
    F_h_(end+1,:)=F_h_(1,:);
    F_h_(end+1,:)=F_h_(2,:);
    F_h_(end+1,:)=F_h_(3,:);
    F_h_(end+1,:)=F_h_(end-4,:);
    F_h_(end+1,:)=F_h_(end-4,:);
    F_h_ = circshift(F_h_,[2 0]);  
else    
    F_h_(end+1,:)=F_h_(1,:);
    F_h_(end+1,:)=F_h_(2,:);
    F_h_(end+1,:)=F_h_(end-3,:);
    F_h_(end+1,:)=F_h_(end-3,:);
    F_h_ = circshift(F_h_,[2 0]);
end



if is_odd_ncx_v==0    
    F_v_(:,end+1)=F_v_(:,1);
    F_v_(:,end+1)=F_v_(:,2);
    F_v_(:,end+1)=F_v_(:,3);
    F_v_(:,end+1)=F_v_(:,end-4);
    F_v_(:,end+1)=F_v_(:,end-4);
    F_v_ = circshift(F_v_,[0 2]);    
else    
    F_v_(:,end+1)=F_v_(:,1);
    F_v_(:,end+1)=F_v_(:,2);
    F_v_(:,end+1)=F_v_(:,end-3);
    F_v_(:,end+1)=F_v_(:,end-3);
    F_v_ = circshift(F_v_,[0 2]);
    
end

if is_odd_ncy_v==0   
    F_v_(end+1,:)=F_v_(1,:);
    F_v_(end+1,:)=F_v_(2,:);
    F_v_(end+1,:)=F_v_(3,:);
    F_v_(end+1,:)=F_v_(end-4,:);
    F_v_(end+1,:)=F_v_(end-4,:);
    F_v_ = circshift(F_v_,[2 0]);  
else    
    F_v_(end+1,:)=F_v_(1,:);
    F_v_(end+1,:)=F_v_(2,:);
    F_v_(end+1,:)=F_v_(end-3,:);
    F_v_(end+1,:)=F_v_(end-3,:);
    F_v_ = circshift(F_v_,[2 0]);
end





% F_h_ = circshift(F_h_,[1 1]);


% clearvars F_h_A F_h_B F_h_C
% 
% F_h_A=F_h_;
% F_h_B=F_h_A;
% 
% % F_h_B(:,end+1)=F_h_A
% 
% F_h_B(:,end+1)=F_h_A(:,1);
% F_h_B(:,end+1)=F_h_A(:,2);
% F_h_B(:,end+1)=F_h_A(:,end);
% F_h_B(end+1,:)=F_h_B(1,:);
% F_h_B(end+1,:)=F_h_B(2,:);
% F_h_B(end+1,:)=F_h_B(end-2,:);
% 
% F_h_C = circshift(F_h_B,[1 1]);

% Augment to n+1 boundaries for interpolation

% F_h_(:,end+1)=F_h_(:,1);
% F_h_(:,end+1)=F_h_(:,2);
% F_h_(:,end+1)=F_h_(:,end-2);
% F_h_(end+1,:)=F_h_(1,:);
% F_h_(end+1,:)=F_h_(2,:);
% F_h_(end+1,:)=F_h_(end-2,:);
% F_h_ = circshift(F_h_,[1 1]);
% 
% F_v_(:,end+1)=F_v_(:,1);
% F_v_(:,end+1)=F_v_(:,2);
% F_v_(:,end+1)=F_v_(:,end-2);
% F_v_(end+1,:)=F_v_(1,:);
% F_v_(end+1,:)=F_v_(2,:);
% F_v_(end+1,:)=F_v_(end-2,:);
% F_v_ = circshift(F_v_,[1 1]);

% 

% if is_odd_ncx_h==0 F_h_A(:,end+1)=F_h_A(:,1);end
% if is_odd_ncy_h==0 F_h_A(end+1,:)=F_h_A(1,:);end
% 
% 
% F_h_B=F_h_A;
% 
% F_h_B
% 
% F_h_B(:,end+1)=F_h_B(:,2);
% F_h_B(:,end+1)=F_h_B(:,3);
% F_h_B(end+1,:)=F_h_B(2,:);
% F_h_B(end+1,:)=F_h_B(3,:);
% 
% 
% F_h_C = circshift(F_h_B,[1 1]);





% if is_odd_ncx_h==0 F_h_(:,end+1)=F_h_(:,1);end
% if is_odd_ncy_h==0 F_h_(end+1,:)=F_h_(1,:);end
% 
% if is_odd_ncx_v==0 F_v_(:,end+1)=F_v_(:,1);end
% if is_odd_ncy_v==0 F_v_(end+1,:)=F_v_(1,:);end


% figure; imagesc(abs(F_h_));colorbar;
% figure; imagesc(abs(F_v_));colorbar;


% Apply interpolation

%Vertical and horizontal lines (double interpolating)


n_thetas=2*(2*(ncx_v+ncy_h));
d_theta=2*pi/n_thetas;

thetas=0:d_theta:2*pi-d_theta;
theta_x_half=atan(floor(ncy_h/2)/floor(ncx_v/2));

aux_count_vert=1;
aux_count_hor=1;
clearvars sample_points_vert sample_points_hor

RPinterp_vert_v_x_cell={};
RPinterp_vert_v_y_cell={};
RPinterp_hor_h_x_cell={};
RPinterp_hor_h_y_cell={};

aux_count_interp_vert=1;
aux_count_interp_vert_interp=1;
aux_count_interp_vert_cell=1;
aux_count_interp_hor=1;
aux_count_interp_hor_interp=1;
aux_count_interp_hor_cell=1;

%para armazenar, ao invez de celulcas, faz uma entrada auxiliar labelling
%de 1 a lenght(t) (um numero para cada angulo)

for t=thetas
    if (t<theta_x_half)|((t>=pi-theta_x_half)&(t<pi+theta_x_half))|(t>=2*pi-theta_x_half)
        for r=-floor(abs(2+ncx_h/(2*cos(t)))):1/2:floor(abs(2+ncx_v/(2*cos(t))))
            sample_points_vert(aux_count_vert,1)=r;
            sample_points_vert(aux_count_vert,2)=t;
            sample_points_vert(aux_count_vert,3)=r*cos(t);
            sample_points_vert(aux_count_vert,4)=r*sin(t);                
            sample_points_vert(aux_count_vert,5)=aux_count_interp_vert_interp;                
            aux_count_vert=aux_count_vert+1;
            RPinterp_vert_v_x(aux_count_interp_vert)=r*cos(t);
            RPinterp_vert_v_y(aux_count_interp_vert)=r*sin(t);
            aux_count_interp_vert=aux_count_interp_vert+1;
        end 
        RPinterp_vert_v_x_cell{aux_count_interp_vert_cell}=RPinterp_vert_v_x;
        RPinterp_vert_v_y_cell{aux_count_interp_vert_cell}=RPinterp_vert_v_y;
        aux_count_interp_vert_cell=aux_count_interp_vert_cell+1;
        aux_count_interp_vert=1;
        clearvars RPinterp_vert_v_x RPinterp_vert_v_y
        aux_count_interp_vert_interp=aux_count_interp_vert_interp+1;
    else
        for r=-floor(2+abs(ncy_h/(2*sin(t)))):1/2:floor(2+abs(ncy_h/(2*sin(t))))
            sample_points_hor(aux_count_hor,1)=r;
            sample_points_hor(aux_count_hor,2)=t;
            sample_points_hor(aux_count_hor,3)=r*cos(t);
            sample_points_hor(aux_count_hor,4)=r*sin(t);            
            sample_points_hor(aux_count_hor,5)=aux_count_interp_hor_interp;            
            aux_count_hor=aux_count_hor+1;
            RPinterp_hor_h_x(aux_count_interp_hor)=r*cos(t);
            RPinterp_hor_h_y(aux_count_interp_hor)=r*sin(t);
            aux_count_interp_hor=aux_count_interp_hor+1;
        end 
        RPinterp_hor_h_x_cell{aux_count_interp_hor_cell}=RPinterp_hor_h_x;
        RPinterp_hor_h_y_cell{aux_count_interp_hor_cell}=RPinterp_hor_h_y;
        aux_count_interp_hor_cell=aux_count_interp_hor_cell+1;
        aux_count_interp_hor=1;
        clearvars RPinterp_hor_h_x RPinterp_hor_h_y
        aux_count_interp_hor_interp=aux_count_interp_hor_interp+1;
    end
end



figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))

figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))


figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
hold on;
scatter(sample_points_hor(:,3),sample_points_hor(:,4))

% %Vertical lines_centered
% 
% 
% for t=1:ncx_v+(1-is_odd_ncx_v)
%     for r=1:ncy_v+(1-is_odd_ncy_v)
%         M=(ncx_v+(1-is_odd_ncx_v)-(2*t)+1)/(ncy_v-1+(1-is_odd_ncy_v));
%         y=min(max(round(r),1),ncy_v+(1-is_odd_ncy_v));
%         x=min(max(fix(t+M*(r-1)-ceil((ncx_v+1)/2))+ceil((ncx_v+1)/2),1),ncx_v+(1-is_odd_ncx_v));    
% %         A((t-1)*ncy+r,1)=x;
% %         A((t-1)*ncy+r,2)=y;
%         RPinterp_vert_v(t,r)=F_v_(x,y);
%     end
% end
% 
% %Horizontal lines_centered
% 
% for t=1:ncy_h+(1-is_odd_ncy_h)
%     for r=1:ncx_h+(1-is_odd_ncx_h)
%         M=(ncy_h+(1-is_odd_ncy_h)-(2*t)+1)/(ncx_h-1+(1-is_odd_ncx_h));
%         x=min(max(round(r),1),ncx_h+(1-is_odd_ncx_h));
%         y=min(max(fix(t+M*(r-1)-ceil((ncy_h+1)/2))+ceil((ncy_h+1)/2),1),ncy_h+(1-is_odd_ncy_h));
% %         B((t-1)*ncx+r,1)=r;
% %         B((t-1)*ncx+r,2)=t+M*(r-1);        
%         RPinterp_hor_h(t,r)=F_h_(x,y);
%     end
% end



% figure; imagesc(abs(RPinterp_hor_h));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% figure; imagesc(abs(RPinterp_vert_v));colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


%Radon Transformation
% 
% RPinterp_vert_v(:,end)=[];
% RPinterp_hor_h(:,end)=[];

% [X,Y] = meshgrid(1:(ncy+1),1:(ncx+1));
% [Xq,Yq] = meshgrid(1:(ncy),1:(ncx));
% Xq=Xq+0.5;
% Yq=Yq+0.5;
% Fq = interp2(X,Y,F_ext,Xq,Yq,'cubic');
% 
% figure; imagesc(abs(Fq));colorbar;
% 
% F_=Fq;

% [x_h,y_h] = meshgrid(-(-1+length(F_h_(1,:)))/2:(-1+length(F_h_(1,:)))/2,-(-1+length(F_h_(:,1)))/2:(-1+length(F_h_(:,1)))/2);
[x_h,y_h] = meshgrid(-(-1+length(F_h_(:,1)))/2:(-1+length(F_h_(:,1)))/2,-(-1+length(F_h_(1,:)))/2:(-1+length(F_h_(1,:)))/2);
RPinterp_hor_h = interp2(x_h,y_h,F_h_',sample_points_hor(:,3),sample_points_hor(:,4),'cubic');

[x_v,y_v] = meshgrid(-(-1+length(F_v_(:,1)))/2:(-1+length(F_v_(:,1)))/2,-(-1+length(F_v_(1,:)))/2:(-1+length(F_v_(1,:)))/2);
RPinterp_vert_v = interp2(x_v,y_v,F_v_',sample_points_vert(:,3),sample_points_vert(:,4),'cubic');


aux_RPinterp_count=1;
Radon_vert_v_aux=[];
Radon_vert_v=[];
for i=1:length(sample_points_vert(:,1))
    if (sample_points_vert(i,5)==aux_RPinterp_count)
        Radon_vert_v_aux(end+1)=RPinterp_vert_v(i);
    else
        Radon_vert_v=[aux_RPinterp_count*ones(size(Radon_vert_v_aux));(ifft(Radon_vert_v_aux.').')];
        display(aux_RPinterp_count)
        
        Radon_vert_v_aux=[];
        aux_RPinterp_count=aux_RPinterp_count+1;
        Radon_vert_v_aux(end+1)=RPinterp_vert_v(i);
    end
end
Radon_vert_v=[aux_RPinterp_count*ones(size(Radon_vert_v_aux));(ifft(Radon_vert_v_aux.').')];




%parei aqui


% Radon_vert_v_aux_=Radon_vert_v_aux;
% Radon_vert_v_aux_(end)=[];
% Radon_vert_v_aux_(end)=[];
% Radon_vert_v_aux_(end)=[];
% Radon_vert_v_aux_(1)=[];
% Radon_vert_v_aux_(1)=[];
% Radon_vert_v=[(ifft(Radon_vert_v_aux_.').')];
% figure; plot(abs(Radon_vert_v))
% figure; plot(real(Radon_vert_v))
% figure; plot(imag(Radon_vert_v))

% 
%     Radon_vert_v
% 
% RPinterp_hor_h



Radon_vert_v = (ifft(RPinterp_vert_v.').');
Radon_hor_h = (ifft(RPinterp_hor_h.').');

for i=1:ncx_v+(1-is_odd_ncx_v)
    for j=1:ncy_v
        Radon_vert_v_(i,j)=(exp(-1i*pi*(j-1)*(1-is_odd_ncy_v/ncy_v)))*Radon_vert_v(i,j);
    end
end

for i=1:ncx_h
    for j=1:ncy_h+(1-is_odd_ncy_h)
        Radon_hor_h_(j,i)=(exp(-1i*pi*(i-1)*(1-is_odd_ncx_h/ncx_h)))*Radon_hor_h(j,i);
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

