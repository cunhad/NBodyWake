function [ Radon_X__,Radon_Y__,Radon_Z__,interp_X_info,interp_Y_info,interp_Z_info,ncx,ncy,ncz] = Ridgelet3d_interp_forwards_dev3_old(X)



% Forward Radon Recto-polar interpolation


% implicit variable

[ncx,ncy,ncz]=size(X);

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);
is_odd_ncz=mod(ncz,2);

%Apply FFT3 (% convention to the frequency center fft)

for i=1:ncx
    for j=1:ncy
        for k=1:ncz
            X_(i,j,k)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*(exp(1i*pi*(k-1)*(1-is_odd_ncz/ncz)))*X(i,j,k);
        end
    end
end


% F_ = fft(fft(X_).').';
% 
% F_ = fft(X_,[],3);

F_ = fftn(X_);


% Augmentation if even, plus four axtra colums (to increase the cover of
% the crude intrerpolation

if is_odd_ncx==0   
    F_(end+1,:,:)=F_(1,:,:);
    F_(end+1,:,:)=F_(2,:,:);
    F_(end+1,:,:)=F_(end-2,:,:);
    F_ = circshift(F_,[1 0 0]);  
else    
    F_(end+1,:,:)=F_(1,:,:);
    F_(end+1,:,:)=F_(end-1,:,:);
    F_ = circshift(F_,[1 0 0]);
end


if is_odd_ncy==0    
    F_(:,end+1,:)=F_(:,1,:);
    F_(:,end+1,:)=F_(:,2,:);
    F_(:,end+1,:)=F_(:,end-2,:);
    F_ = circshift(F_,[0 1 0]);    
else    
    F_(:,end+1,:)=F_(:,1,:);
    F_(:,end+1,:)=F_(:,end-1,:);
    F_ = circshift(F_,[0 1 0]);
    
end

if is_odd_ncz==0    
    F_(:,:,end+1)=F_(:,:,1);
    F_(:,:,end+1)=F_(:,:,2);
    F_(:,:,end+1)=F_(:,:,end-2);
    F_ = circshift(F_,[0 0 1]);    
else    
    F_(:,:,end+1)=F_(:,:,1);
    F_(:,:,end+1)=F_(:,:,end-1);
    F_ = circshift(F_,[0 0 1]);
    
end



% figure; imagesc(real(F_')); colorbar;
% figure; imagesc(imag(F_')); colorbar;

% Choose the points to be interpolated


% Piramid Z

D_theta_xy=2*atan((ncx/2)/(ncz/2));
D_theta_yx=2*atan((ncy/2)/(ncz/2));

n_thetas_xy=2*(2*(ncx));
n_thetas_yx=2*(2*(ncy));

d_theta_xy=D_theta_xy/n_thetas_xy;
d_theta_yx=D_theta_yx/n_thetas_yx;

thetas_xy=0:d_theta_xy:D_theta_xy-d_theta_xy*0;
thetas_xy=thetas_xy-D_theta_xy/2;

thetas_yx=0:d_theta_yx:D_theta_yx-d_theta_yx*0;
thetas_yx=thetas_yx-D_theta_yx/2;


aux_count_Z=1;
aux_count_interp_Z_interp=1;
aux_count_interp_Z=1;
% aux_count_interp_Z_cell=1;

aux_count_t1 = 1;
aux_count_t2 = 1;


for t1 = thetas_xy
    for t2 = thetas_yx        
        half_diag = floor((ncz/2)*sqrt((tan(t1)^2)+(tan(t2)^2)+1));
        spherical_theta = atan(sqrt((tan(t1)^2)+(tan(t2)^2)));
        spherical_phi = atan(tan(t2)/tan(t1))+double(t1<0)*pi;
        if (t1==0&t2==0) 
            spherical_phi = 0;
        end
        for r = -half_diag:half_diag
            
%             sample_points_Z(aux_count_Z,1)=spherical_theta;
%             sample_points_Z(aux_count_Z,2)=spherical_phi;
%             sample_points_Z(aux_count_Z,3)=r;
%             sample_points_Z(aux_count_Z,4)=r*sin(spherical_theta)*cos(spherical_phi);
%             sample_points_Z(aux_count_Z,5)=r*sin(spherical_theta)*sin(spherical_phi);
%             sample_points_Z(aux_count_Z,6)=r*cos(spherical_theta);
%             sample_points_Z(aux_count_Z,7)=aux_count_interp_Z_interp;
            
            sample_points_Z(1,aux_count_Z)=spherical_theta;
            sample_points_Z(2,aux_count_Z)=spherical_phi;
            sample_points_Z(3,aux_count_Z)=r;
            sample_points_Z(4,aux_count_Z)=r*sin(spherical_theta)*cos(spherical_phi);
            sample_points_Z(5,aux_count_Z)=r*sin(spherical_theta)*sin(spherical_phi);
            sample_points_Z(6,aux_count_Z)=r*cos(spherical_theta);
            sample_points_Z(7,aux_count_Z)=aux_count_interp_Z_interp;
            
            
            aux_count_Z=aux_count_Z+1;   
            
%             interp_Z_x(aux_count_interp_Z)=r*sin(spherical_theta)*cos(spherical_phi);
%             interp_Z_y(aux_count_interp_Z)=r*sin(spherical_theta)*sin(spherical_phi);
%             interp_Z_z(aux_count_interp_Z)=r*cos(spherical_theta);
            aux_count_interp_Z=aux_count_interp_Z+1;
        end
        
        aux_count_interp_Z_interp = aux_count_interp_Z_interp +1;
        
        interp_Z_info(aux_count_t1,aux_count_t2)=aux_count_interp_Z-1;        
        aux_count_t2 = aux_count_t2 +1;

        
%         interp_Z_x_cell{aux_count_interp_Z_cell}=interp_Z_x;
%         interp_Z_y_cell{aux_count_interp_Z_cell}=interp_Z_y;
%         interp_Z_z_cell{aux_count_interp_Z_cell}=interp_Z_z;%         
%         aux_count_interp_Z_cell=aux_count_interp_Z_cell+1;
        
        aux_count_interp_Z = 1;
        
%         clearvars interp_Z_x interp_Z_y interp_Z_z
        
        
    end
    
    aux_count_t2 = 1;
    aux_count_t1 = aux_count_t1+1;
end


% sample_points_Z = sample_points_Z';


% 
% figure;
% scatter3(sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6))








% % Piramid Y  (Z->Y,Y->-Z with respect to piramid Z)
% 
% D_theta_xz=2*atan((ncx/2)/(ncy/2));
% D_theta_zx=2*atan((ncz/2)/(ncy/2));
% 
% n_thetas_xz=2*(2*(ncx));
% n_thetas_zx=2*(2*(ncz));
% 
% d_theta_xz=D_theta_xz/n_thetas_xz;
% d_theta_zx=D_theta_zx/n_thetas_zx;
% 
% thetas_xz=0:d_theta_xz:D_theta_xz-d_theta_xz*0;
% thetas_xz=thetas_xz-D_theta_xz/2;
% 
% thetas_zx=0:d_theta_zx:D_theta_zx-d_theta_zx*0;
% thetas_zx=thetas_zx-D_theta_zx/2;
% 
% 
% aux_count_Y=1;
% aux_count_interp_Y_interp=1;
% aux_count_interp_Y=1;
% % aux_count_interp_Y_cell=1;
% 
% 
% aux_count_t1 = 1;
% aux_count_t2 = 1;
% 
% for t1 = thetas_xz
%     for t2 = thetas_zx        
%         half_diag = floor((ncy/2)*sqrt((tan(t1)^2)+(tan(t2)^2)+1));
%         spherical_theta = atan(sqrt((tan(t1)^2)+(tan(t2)^2)));
%         spherical_phi = atan(tan(t2)/tan(t1))+double(t1<0)*pi;
%         if (t1==0&t2==0) 
%             spherical_phi = 0;
%         end        
%         for r = -half_diag:half_diag
%             sample_points_Y(aux_count_Y,1)=spherical_theta;
%             sample_points_Y(aux_count_Y,2)=spherical_phi;
%             sample_points_Y(aux_count_Y,3)=r;
%             sample_points_Y(aux_count_Y,4)=r*sin(spherical_theta)*cos(spherical_phi);
%             sample_points_Y(aux_count_Y,5)=r*cos(spherical_theta);
%             sample_points_Y(aux_count_Y,6)=-r*sin(spherical_theta)*sin(spherical_phi);
%             sample_points_Y(aux_count_Y,7)=aux_count_interp_Y_interp;
%             aux_count_Y=aux_count_Y+1;   
%             
% %             interp_Y_x(aux_count_interp_Y)=r*sin(spherical_theta)*cos(spherical_phi);
% %             interp_Y_y(aux_count_interp_Y)=r*cos(spherical_theta);
% %             interp_Y_z(aux_count_interp_Y)=-r*sin(spherical_theta)*sin(spherical_phi);
%             aux_count_interp_Y=aux_count_interp_Y+1;
%         end
%         
%         aux_count_interp_Y_interp = aux_count_interp_Y_interp +1;      
%         
%         interp_Y_info(aux_count_t1,aux_count_t2)=aux_count_interp_Y-1;        
%         aux_count_t2 = aux_count_t2 +1;
%         
% %         interp_Y_x_cell{aux_count_interp_Y_cell}=interp_Y_x;
% %         interp_Y_y_cell{aux_count_interp_Y_cell}=interp_Y_y;
% %         interp_Y_z_cell{aux_count_interp_Y_cell}=interp_Y_z;
% %         aux_count_interp_Y_cell=aux_count_interp_Y_cell+1;
%         
%         aux_count_interp_Y = 1;
%         
% %         clearvars interp_Y_x interp_Y_y interp_Y_z
%         
%     end
%     
%     aux_count_t2 = 1;
%     aux_count_t1 = aux_count_t1+1;
% end
% 
% % 
% % figure;
% % scatter3(sample_points_Y(:,4),sample_points_Y(:,5),sample_points_Y(:,6))
% 
% 
% 
% 
% % Piramid X  (Z->X,X->-Z with respect to piramid Z)
% 
% 
% D_theta_zy=2*atan((ncz/2)/(ncx/2));
% D_theta_yz=2*atan((ncy/2)/(ncx/2));
% 
% n_thetas_zy=2*(2*(ncz));
% n_thetas_yz=2*(2*(ncy));
% 
% d_theta_zy=D_theta_zy/n_thetas_zy;
% d_theta_yz=D_theta_yz/n_thetas_yz;
% 
% thetas_zy=0:d_theta_zy:D_theta_zy-d_theta_zy*0;
% thetas_zy=thetas_zy-D_theta_zy/2;
% 
% thetas_yz=0:d_theta_yz:D_theta_yz-d_theta_yz*0;
% thetas_yz=thetas_yz-D_theta_yz/2;
% 
% 
% aux_count_X=1;
% aux_count_interp_X_interp=1;
% aux_count_interp_X=1;
% % aux_count_interp_X_cell=1;
% 
% aux_count_t1 = 1;
% aux_count_t2 = 1;
% 
% for t1 = thetas_zy
%     for t2 = thetas_yz        
%         half_diag = floor((ncx/2)*sqrt((tan(t1)^2)+(tan(t2)^2)+1));
%         spherical_theta = atan(sqrt((tan(t1)^2)+(tan(t2)^2)));
%         spherical_phi = atan(tan(t2)/tan(t1))+double(t1<0)*pi;
%         if (t1==0&t2==0) 
%             spherical_phi = 0;
%         end
%         for r = -half_diag:half_diag
%             sample_points_X(aux_count_X,1)=spherical_theta;
%             sample_points_X(aux_count_X,2)=spherical_phi;
%             sample_points_X(aux_count_X,3)=r;
%             sample_points_X(aux_count_X,4)=r*cos(spherical_theta);
%             sample_points_X(aux_count_X,5)=r*sin(spherical_theta)*sin(spherical_phi);
%             sample_points_X(aux_count_X,6)=-r*sin(spherical_theta)*cos(spherical_phi);
%             sample_points_X(aux_count_X,7)=aux_count_interp_X_interp;
%             aux_count_X=aux_count_X+1;   
%             
% %             interp_X_x(aux_count_interp_X)=r*cos(spherical_theta);
% %             interp_X_y(aux_count_interp_X)=r*sin(spherical_theta)*sin(spherical_phi);
% %             interp_X_z(aux_count_interp_X)=-r*sin(spherical_theta)*cos(spherical_phi);
%             aux_count_interp_X=aux_count_interp_X+1;
%         end
%         
%         aux_count_interp_X_interp = aux_count_interp_X_interp +1;
%         
%         interp_X_info(aux_count_t1,aux_count_t2)=aux_count_interp_X-1;        
%         aux_count_t2 = aux_count_t2 +1;     
%         
% %         interp_X_x_cell{aux_count_interp_X_cell}=interp_X_x;
% %         interp_X_y_cell{aux_count_interp_X_cell}=interp_X_y;
% %         interp_X_z_cell{aux_count_interp_X_cell}=interp_X_z;
% %         aux_count_interp_X_cell=aux_count_interp_X_cell+1;
%         
%         aux_count_interp_X = 1;
%         
% %         clearvars interp_X_x interp_X_y interp_X_z
%         
%     end
%     aux_count_t2 = 1;
%     aux_count_t1 = aux_count_t1+1;
% end
% 
% % 
% % figure;
% % scatter3(sample_points_X(:,4),sample_points_X(:,5),sample_points_X(:,6))



% 
% 
% % figure;
% % scatter3(sample_points_X(:,4),sample_points_X(:,5),sample_points_X(:,6))
% % hold on 
% % scatter3(sample_points_Y(:,4),sample_points_Y(:,5),sample_points_Y(:,6))
% % scatter3(sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6))
% 
% 
% 
% 
% % 
% % % Record number on points in each line
% % 
% % for t=1:length(interp_X_x_cell)
% %     interp_Z_cell_info(t)=length(interp_X_x_cell{t});
% % end
% % 
% % for t=1:length(interp_Y_x_cell)
% %     interp_Y_cell_info(t)=length(interp_Y_x_cell{t});
% % end
% % 
% % for t=1:length(interp_X_x_cell)
% %     interp_X_cell_info(t)=length(interp_X_x_cell{t});
% % end




%Apply (cubic) interpolation


[x_,y_,z_] = meshgrid(-(-1+length(F_(:,1,1)))/2:(-1+length(F_(:,1,1)))/2,-(-1+length(F_(1,:,1)))/2:(-1+length(F_(1,:,1)))/2,-(-1+length(F_(1,1,:)))/2:(-1+length(F_(1,1,:)))/2);
% [y_,x_,z_] = meshgrid(-(-1+length(F_(1,:,1)))/2:(-1+length(F_(1,:,1)))/2,-(-1+length(F_(:,1,1)))/2:(-1+length(F_(:,1,1)))/2,-(-1+length(F_(1,1,:)))/2:(-1+length(F_(1,1,:)))/2);

% x_ = permute(x_, [2 1 3]);
% y_ = permute(y_, [2 1 3]);
% z_ = permute(z_, [2 1 3]);

% interp_X = interp3(x_,y_,z_,permute(F_, [2 1 3]),sample_points_X(:,4),sample_points_X(:,5),sample_points_X(:,6),'cubic');
% interp_Y = interp3(x_,y_,z_,permute(F_, [2 1 3]),sample_points_Y(:,4),sample_points_Y(:,5),sample_points_Y(:,6),'cubic');
% interp_Z = interp3(x_,y_,z_,permute(F_, [2 1 3]),sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'cubic');
interp_Z = interp3(x_,y_,z_,permute(F_, [2 1 3]),sample_points_Z(4,:),sample_points_Z(5,:),sample_points_Z(6,:),'cubic');


% interp_X = interp3(x_,y_,z_,F_,sample_points_X(:,4),sample_points_X(:,5),sample_points_X(:,6),'cubic');
% interp_Y = interp3(x_,y_,z_,F_,sample_points_Y(:,4),sample_points_Y(:,5),sample_points_Y(:,6),'cubic');
% interp_Z = interp3(x_,y_,z_,F_,sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'cubic');



% 
% interp_Z = interp3(x_,y_,z_,F_,sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'spline');
% interp_Z = interp3(x_,y_,z_,F_,sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'makima');
% interp_Z = interp3(x_,y_,z_,F_,sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'nearest');
% interp_Z = interp3(x_,y_,z_,F_,sample_points_Z(:,4),sample_points_Z(:,5),sample_points_Z(:,6),'linear');



% Obtain Radon by inverse fft of the the interpolations



% aux_interp_X_count=1;
% Radon_X_aux=[];
% Radon_X=[];
% 
% for i=1:length(sample_points_X(:,1))
%     if (sample_points_X(i,7)==aux_interp_X_count)
%         Radon_X_aux(end+1)=interp_X(i);
%     else
%         
%         Radon_X=[Radon_X,(ifft(Radon_X_aux.').')];
%         Radon_X_aux=[];       
%         aux_interp_X_count=aux_interp_X_count+1;      
%         Radon_X_aux(end+1)=interp_X(i);
% 
%     end
% end
% 
% Radon_X=[Radon_X,(ifft(Radon_X_aux.').')];
% 
% 
% aux_interp_Y_count=1;
% Radon_Y_aux=[];
% Radon_Y=[];
% 
% for i=1:length(sample_points_Y(:,1))
%     if (sample_points_Y(i,7)==aux_interp_Y_count)
%         Radon_Y_aux(end+1)=interp_Y(i);
%     else
%         
%         Radon_Y=[Radon_Y,(ifft(Radon_Y_aux.').')];
%         Radon_Y_aux=[];       
%         aux_interp_Y_count=aux_interp_Y_count+1;      
%         Radon_Y_aux(end+1)=interp_Y(i);
% 
%     end
% end
% 
% Radon_Y=[Radon_Y,(ifft(Radon_Y_aux.').')];

% aux_interp_Z_count=1;
% Radon_Z_aux=[];
% Radon_Z=[];
% 
% for i=1:length(sample_points_Z(:,1))
%     if (sample_points_Z(i,7)==aux_interp_Z_count)
%         Radon_Z_aux(end+1)=interp_Z(i);
%     else
%         
%         Radon_Z=[Radon_Z,(ifft(Radon_Z_aux.').')];
%         Radon_Z_aux=[];       
%         aux_interp_Z_count=aux_interp_Z_count+1;      
%         Radon_Z_aux(end+1)=interp_Z(i);
% 
%     end
% end
% 
% Radon_Z=[Radon_Z,(ifft(Radon_Z_aux.').')];
% 
aux_interp_Z_count=1;
Radon_Z_aux=[];
Radon_Z=[];



for i=1:length(sample_points_Z(1,:))
    if (sample_points_Z(7,i)==aux_interp_Z_count)
        Radon_Z_aux(end+1)=interp_Z(i);
    else
        
        Radon_Z=[Radon_Z,(ifft(Radon_Z_aux.').')];
        Radon_Z_aux=[];       
        aux_interp_Z_count=aux_interp_Z_count+1;      
        Radon_Z_aux(end+1)=interp_Z(i);

    end
end

Radon_Z=[Radon_Z,(ifft(Radon_Z_aux.').')];

% convention to the frequency center fft




% aux_interp_X_count__=1;
% Radon_X_aux__=[];
% Radon_X__=[];
% for i=1:length(sample_points_X(:,1))
%     if (sample_points_X(i,7)==aux_interp_X_count__)
%         Radon_X_aux__(end+1)=Radon_X(i);
%     else
%         ncr = length(Radon_X_aux__);
%         for r=1:ncr
%             Radon_X_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_X_aux__(r);
%         end
%         
%         Radon_X__=[Radon_X__,Radon_X_aux__];
%         Radon_X_aux__=[];  
%         aux_interp_X_count__=aux_interp_X_count__+1;
%         Radon_X_aux__(end+1)=Radon_X(i);
%     end
% end
% ncr = length(Radon_X_aux__);
% for r=1:ncr
%     Radon_X_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_X_aux__(r);
% end
% Radon_X__=[Radon_X__,Radon_X_aux__];
% 
% 
% aux_interp_Y_count__=1;
% Radon_Y_aux__=[];
% Radon_Y__=[];
% for i=1:length(sample_points_Y(:,1))
%     if (sample_points_Y(i,7)==aux_interp_Y_count__)
%         Radon_Y_aux__(end+1)=Radon_Y(i);
%     else
%         ncr = length(Radon_Y_aux__);
%         for r=1:ncr
%             Radon_Y_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Y_aux__(r);
%         end
%         
%         Radon_Y__=[Radon_Y__,Radon_Y_aux__];
%         Radon_Y_aux__=[];  
%         aux_interp_Y_count__=aux_interp_Y_count__+1;
%         Radon_Y_aux__(end+1)=Radon_Y(i);
%     end
% end
% ncr = length(Radon_Y_aux__);
% for r=1:ncr
%     Radon_Y_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Y_aux__(r);
% end
% Radon_Y__=[Radon_Y__,Radon_Y_aux__];




% aux_interp_Z_count__=1;
% Radon_Z_aux__=[];
% Radon_Z__=[];
% for i=1:length(sample_points_Z(:,1))
%     if (sample_points_Z(i,7)==aux_interp_Z_count__)
%         Radon_Z_aux__(end+1)=Radon_Z(i);
%     else
%         ncr = length(Radon_Z_aux__);
%         for r=1:ncr
%             Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
%         end
%         
%         Radon_Z__=[Radon_Z__,Radon_Z_aux__];
%         Radon_Z_aux__=[];  
%         aux_interp_Z_count__=aux_interp_Z_count__+1;
%         Radon_Z_aux__(end+1)=Radon_Z(i);
%     end
% end
% ncr = length(Radon_Z_aux__);
% for r=1:ncr
%     Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
% end
% Radon_Z__=[Radon_Z__,Radon_Z_aux__];


aux_interp_Z_count__=1;
Radon_Z_aux__=[];
Radon_Z__=[];
for i=1:length(sample_points_Z(1,:))
    if (sample_points_Z(7,i)==aux_interp_Z_count__)
        Radon_Z_aux__(end+1)=Radon_Z(i);
    else
        ncr = length(Radon_Z_aux__);
        for r=1:ncr
            Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
        end
        
        Radon_Z__=[Radon_Z__,Radon_Z_aux__];
        Radon_Z_aux__=[];  
        aux_interp_Z_count__=aux_interp_Z_count__+1;
        Radon_Z_aux__(end+1)=Radon_Z(i);
    end
end
ncr = length(Radon_Z_aux__);
for r=1:ncr
    Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
end
Radon_Z__=[Radon_Z__,Radon_Z_aux__];



% Store also info related to the interpolation, so the Radon Transformation
% can be properly plotted



Radon_X__ = Radon_Z__;
Radon_Y__ = Radon_Z__;

% Radon_X__(2,:)=sample_points_X(:,7);
% 
% Radon_Y__(2,:)=sample_points_Y(:,7);


% Radon_X__(2,:)=sample_points_Z(:,7);
% 
% Radon_Y__(2,:)=sample_points_Z(:,7);
% 
% Radon_Z__(2,:)=sample_points_Z(:,7);


Radon_X__(2,:)=sample_points_Z(7,:);

Radon_Y__(2,:)=sample_points_Z(7,:);

Radon_Z__(2,:)=sample_points_Z(7,:);



interp_X_info = interp_Z_info;
interp_Y_info = interp_Z_info;

end