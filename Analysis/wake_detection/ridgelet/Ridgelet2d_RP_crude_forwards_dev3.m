function [ Radon_hor__,Radon_vert__,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy] = Ridgelet2d_RP_crude_forwards_dev3(X)



% Forward Radon Recto-polar interpolation


% implicit variable

[ncx,ncy]=size(X);

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);


%Apply FFT2 (% convention to the frequency center fft)

for i=1:ncx
    for j=1:ncy
        X_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X(i,j);
    end
end


F_ = fft(fft(X_).').';



% Augmentation if even, plus four axtra colums (to increase the cover of
% the crude intrerpolation

% 
% if is_odd_ncy==0    
%     F_(:,end+1)=F_(:,1);
%     F_(:,end+1)=F_(:,2);
%     F_(:,end+1)=F_(:,3);
%     F_(:,end+1)=F_(:,end-4);
%     F_(:,end+1)=F_(:,end-4);
%     F_ = circshift(F_,[0 2]);    
% else    
%     F_(:,end+1)=F_(:,1);
%     F_(:,end+1)=F_(:,2);
%     F_(:,end+1)=F_(:,end-3);
%     F_(:,end+1)=F_(:,end-3);
%     F_ = circshift(F_,[0 2]);
%     
% end
% 
% 
% 
% if is_odd_ncx==0   
%     F_(end+1,:)=F_(1,:);
%     F_(end+1,:)=F_(2,:);
%     F_(end+1,:)=F_(3,:);
%     F_(end+1,:)=F_(end-4,:);
%     F_(end+1,:)=F_(end-4,:);
%     F_ = circshift(F_,[2 0]);  
% else    
%     F_(end+1,:)=F_(1,:);
%     F_(end+1,:)=F_(2,:);
%     F_(end+1,:)=F_(end-3,:);
%     F_(end+1,:)=F_(end-3,:);
%     F_ = circshift(F_,[2 0]);
% end



if is_odd_ncy==0    
    F_(:,end+1)=F_(:,1);
    F_(:,end+1)=F_(:,2);
%     F_(:,end+1)=F_(:,3);
    F_(:,end+1)=F_(:,end-2);
%     F_(:,end+1)=F_(:,end-4);
    F_ = circshift(F_,[0 1]);    
else    
    F_(:,end+1)=F_(:,1);
%     F_(:,end+1)=F_(:,2);
%     F_(:,end+1)=F_(:,end-3);
    F_(:,end+1)=F_(:,end-1);
    F_ = circshift(F_,[0 1]);
    
end



if is_odd_ncx==0   
    F_(end+1,:)=F_(1,:);
    F_(end+1,:)=F_(2,:);
%     F_(end+1,:)=F_(3,:);
%     F_(end+1,:)=F_(end-4,:);
    F_(end+1,:)=F_(end-2,:);
    F_ = circshift(F_,[1 0]);  
else    
    F_(end+1,:)=F_(1,:);
%     F_(end+1,:)=F_(2,:);
%     F_(end+1,:)=F_(end-3,:);
    F_(end+1,:)=F_(end-1,:);
    F_ = circshift(F_,[1 0]);
end


% figure; imagesc(real(F_')); colorbar;
% figure; imagesc(imag(F_')); colorbar;

% Choose the points to be interpolated

%Vertical and horizontal lines 



D_theta_x=2*atan((ncx/2)/(ncy/2));
D_theta_y=2*atan((ncy/2)/(ncx/2));

n_thetas_x=2*(2*(ncx));
n_thetas_y=2*(2*(ncy));

d_theta_x=D_theta_x/n_thetas_x;
d_theta_y=D_theta_y/n_thetas_y;


thetas_y=0:d_theta_y:D_theta_y-d_theta_y;
thetas_y=thetas_y-D_theta_y/2;

thetas_x=0:d_theta_x:D_theta_x-d_theta_x;
thetas_x=thetas_x+D_theta_y/2;


thetas=[thetas_y,thetas_x];




aux_count_vert=1;
aux_count_hor=1;
aux_count_interp_hor_interp=1;
aux_count_interp_vert_interp=1;

aux_count_interp_hor=1;
aux_count_interp_hor_cell=1;
aux_count_interp_vert=1;
aux_count_interp_vert_cell=1;

for t=thetas
    
    if aux_count_interp_hor_interp<=n_thetas_y
        
        
        for r=-floor(abs(ncx/(2*cos(t)))):1:floor(abs(ncx/(2*cos(t))))
            
            sample_points_hor(aux_count_hor,1)=r;
            sample_points_hor(aux_count_hor,2)=t;
            sample_points_hor(aux_count_hor,3)=r*cos(t);
            sample_points_hor(aux_count_hor,4)=r*sin(t);
            sample_points_hor(aux_count_hor,5)=aux_count_interp_hor_interp;
            aux_count_hor=aux_count_hor+1;   
            
            RPinterp_hor_x(aux_count_interp_hor)=r*cos(t);
            RPinterp_hor_y(aux_count_interp_hor)=r*sin(t);
            aux_count_interp_hor=aux_count_interp_hor+1;
        end
        
        aux_count_interp_hor_interp=aux_count_interp_hor_interp+1;
        
        RPinterp_hor_x_cell{aux_count_interp_hor_cell}=RPinterp_hor_x;
        RPinterp_hor_y_cell{aux_count_interp_hor_cell}=RPinterp_hor_y;
        aux_count_interp_hor_cell=aux_count_interp_hor_cell+1;
        aux_count_interp_hor=1;
        clearvars RPinterp_hor_x RPinterp_hor_y
        
    else
        
        for r=-floor(abs(ncy/(2*sin(t)))):1:floor(abs(ncy/(2*sin(t))))
            sample_points_vert(aux_count_vert,1)=r;
            sample_points_vert(aux_count_vert,2)=t;
            sample_points_vert(aux_count_vert,3)=r*cos(t);
            sample_points_vert(aux_count_vert,4)=r*sin(t);
            sample_points_vert(aux_count_vert,5)=aux_count_interp_vert_interp;
            aux_count_vert=aux_count_vert+1;
            
            RPinterp_vert_x(aux_count_interp_vert)=r*cos(t);
            RPinterp_vert_y(aux_count_interp_vert)=r*sin(t);
            aux_count_interp_vert=aux_count_interp_vert+1; 
            
        end

        aux_count_interp_vert_interp=aux_count_interp_vert_interp+1;
        
        RPinterp_vert_x_cell{aux_count_interp_vert_cell}=RPinterp_vert_x;
        RPinterp_vert_y_cell{aux_count_interp_vert_cell}=RPinterp_vert_y;
        aux_count_interp_vert_cell=aux_count_interp_vert_cell+1;
        aux_count_interp_vert=1;
        clearvars RPinterp_vert_x RPinterp_vert_y

    end
end

%  figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
%  
%  xlim([-(ncx)/2 (ncx)/2])
%  ylim([-(ncy)/2 (ncy)/2])
%  
%  
%  hold on
%   scatter(sample_points_hor(:,3),sample_points_hor(:,4))
%  colorbar;
%  
%  

% Record number on points in each line

for t=1:length(RPinterp_vert_x_cell)
    RPinterp_vert_cell_info(t)=length(RPinterp_vert_x_cell{t});
end

for t=1:length(RPinterp_hor_x_cell)
    RPinterp_hor_cell_info(t)=length(RPinterp_hor_x_cell{t});
end


%Apply (cubic) interpolation


[x_,y_] = meshgrid(-(-1+length(F_(:,1)))/2:(-1+length(F_(:,1)))/2,-(-1+length(F_(1,:)))/2:(-1+length(F_(1,:)))/2);

RPinterp_hor = interp2(x_,y_,F_',sample_points_hor(:,3),sample_points_hor(:,4),'cubic');
RPinterp_vert = interp2(x_,y_,F_',sample_points_vert(:,3),sample_points_vert(:,4),'cubic');

% RPinterp_hor = interp2(x_,y_,F_',sample_points_hor(:,3),sample_points_hor(:,4),'spline');
% RPinterp_vert = interp2(x_,y_,F_',sample_points_vert(:,3),sample_points_vert(:,4),'spline');

% RPinterp_hor = interp2(x_,y_,F_',sample_points_hor(:,3),sample_points_hor(:,4),'makima');
% RPinterp_vert = interp2(x_,y_,F_',sample_points_vert(:,3),sample_points_vert(:,4),'makima');

% RPinterp_hor = interp2(x_,y_,F_',sample_points_hor(:,3),sample_points_hor(:,4),'nearest');
% RPinterp_vert = interp2(x_,y_,F_',sample_points_vert(:,3),sample_points_vert(:,4),'nearest');

% RPinterp_hor = interp2(x_,y_,F_',sample_points_hor(:,3),sample_points_hor(:,4),'linear');
% RPinterp_vert = interp2(x_,y_,F_',sample_points_vert(:,3),sample_points_vert(:,4),'linear');




% Obtain Radon by inverse fft of the the interpolations



aux_RPinterp_vert_count=1;
Radon_vert_aux=[];
% Radon_vert_aux2=[];
Radon_vert=[];
% Radon_vert2=[];


for i=1:length(sample_points_vert(:,1))
    if (sample_points_vert(i,5)==aux_RPinterp_vert_count)
        Radon_vert_aux(end+1)=RPinterp_vert(i);
%         Radon_vert_aux2(end+1)=RPinterp_vert(i);
    else%    
%         Radon_vert_aux(end+1)=RPinterp_vert(i);  
%         Radon_vert_aux2(end+1)=RPinterp_vert(i);
        
        Radon_vert=[Radon_vert,(ifft(Radon_vert_aux.').')];
%         Radon_vert2=[Radon_vert2,Radon_vert_aux2];
%         Radon_vert=[aux_RPinterp_count*ones(size(Radon_vert_aux));(ifft(Radon_vert_aux.').')];
%         display(aux_RPinterp_vert_count);
        
        Radon_vert_aux=[];  
%         Radon_vert_aux2=[];
        
        aux_RPinterp_vert_count=aux_RPinterp_vert_count+1;
        
        Radon_vert_aux(end+1)=RPinterp_vert(i);
%         Radon_vert_aux2(end+1)=RPinterp_vert(i);

    end
end

Radon_vert=[Radon_vert,(ifft(Radon_vert_aux.').')];
% Radon_vert2=[Radon_vert2,Radon_vert_aux2];
% display(aux_RPinterp_vert_count);



aux_RPinterp_hor_count=1;
Radon_hor_aux=[];
Radon_hor=[];
for i=1:length(sample_points_hor(:,1))
    if (sample_points_hor(i,5)==aux_RPinterp_hor_count)
        Radon_hor_aux(end+1)=RPinterp_hor(i);
    else
       
        Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
%         display(aux_RPinterp_hor_count);
        
        Radon_hor_aux=[];  
        
        aux_RPinterp_hor_count=aux_RPinterp_hor_count+1;
        
        Radon_hor_aux(end+1)=RPinterp_hor(i);

    end
end

Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
% display(aux_RPinterp_hor_count);



% convention to the frequency center fft




aux_RPinterp_hor_count__=1;
Radon_hor_aux__=[];
Radon_hor__=[];
for i=1:length(sample_points_hor(:,1))
    if (sample_points_hor(i,5)==aux_RPinterp_hor_count__)
        Radon_hor_aux__(end+1)=Radon_hor(i);
    else
        ncr = length(Radon_hor_aux__);
        for r=1:ncr
            Radon_hor_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_hor_aux__(r);
        end
        
        Radon_hor__=[Radon_hor__,Radon_hor_aux__];
%         display(aux_RPinterp_hor_count__);
        Radon_hor_aux__=[];  
        aux_RPinterp_hor_count__=aux_RPinterp_hor_count__+1;
        Radon_hor_aux__(end+1)=Radon_hor(i);
    end
end


ncr = length(Radon_hor_aux__);
for r=1:ncr
    Radon_hor_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_hor_aux__(r);
end
Radon_hor__=[Radon_hor__,Radon_hor_aux__];
% display(aux_RPinterp_hor_count__);

% ncr = length(Radon_hor_test);
% for r=1:ncr
%     Radon_hor_test__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_hor_test(r);
% end



aux_RPinterp_vert_count__=1;
Radon_vert_aux__=[];
Radon_vert__=[];
for i=1:length(sample_points_vert(:,1))
    if (sample_points_vert(i,5)==aux_RPinterp_vert_count__)
        Radon_vert_aux__(end+1)=Radon_vert(i);
    else
        ncr = length(Radon_vert_aux__);
        for r=1:ncr
            Radon_vert_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_vert_aux__(r);
        end
        
        Radon_vert__=[Radon_vert__,Radon_vert_aux__];
%         display(aux_RPinterp_vert_count__);
        Radon_vert_aux__=[];  
        aux_RPinterp_vert_count__=aux_RPinterp_vert_count__+1;
        Radon_vert_aux__(end+1)=Radon_vert(i);
    end
end


ncr = length(Radon_vert_aux__);
for r=1:ncr
    Radon_vert_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_vert_aux__(r);
end
Radon_vert__=[Radon_vert__,Radon_vert_aux__];
% display(aux_RPinterp_hor_count__);





% Store also info related to the interpolation, so the Radon Transformation
% can be properly plotted


Radon_hor__(2,:)=sample_points_hor(:,5);

Radon_vert__(2,:)=sample_points_vert(:,5);



end
