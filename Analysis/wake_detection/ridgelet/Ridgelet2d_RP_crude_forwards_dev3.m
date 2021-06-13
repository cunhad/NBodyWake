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

% if is_odd_ncy==0 F_(:,end+1)=F_(:,1);end
% if is_odd_ncx==0 F_(end+1,:)=F_(1,:);end

if is_odd_ncy==0    
    F_(:,end+1)=F_(:,1);
    F_(:,end+1)=F_(:,2);
    F_(:,end+1)=F_(:,3);
    F_(:,end+1)=F_(:,end-4);
    F_(:,end+1)=F_(:,end-4);
    F_ = circshift(F_,[0 2]);    
else    
    F_(:,end+1)=F_(:,1);
    F_(:,end+1)=F_(:,2);
    F_(:,end+1)=F_(:,end-3);
    F_(:,end+1)=F_(:,end-3);
    F_ = circshift(F_,[0 2]);
    
end



if is_odd_ncx==0   
    F_(end+1,:)=F_(1,:);
    F_(end+1,:)=F_(2,:);
    F_(end+1,:)=F_(3,:);
    F_(end+1,:)=F_(end-4,:);
    F_(end+1,:)=F_(end-4,:);
    F_ = circshift(F_,[2 0]);  
else    
    F_(end+1,:)=F_(1,:);
    F_(end+1,:)=F_(2,:);
    F_(end+1,:)=F_(end-3,:);
    F_(end+1,:)=F_(end-3,:);
    F_ = circshift(F_,[2 0]);
end

% figure; imagesc(real(F_')); colorbar;
% figure; imagesc(imag(F_')); colorbar;

% Apply interpolation

%Vertical and horizontal lines (double interpolating)




D_theta_x=2*atan((ncx/2)/(ncy/2));
D_theta_y=2*atan((ncy/2)/(ncx/2));

n_thetas_x=2*(2*(ncx));
n_thetas_y=2*(2*(ncy));

d_theta_x=D_theta_x/n_thetas_x;
d_theta_y=D_theta_y/n_thetas_y;


% thetas_x=0:d_theta_x:D_theta_x-d_theta_x;
% thetas_x=thetas_x-D_theta_x/2;

% thetas_y=0:d_theta_y:D_theta_y-d_theta_y;
% thetas_y=thetas_y+D_theta_x/2;

% thetas=[thetas_x,thetas_y];


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
            
%             //going much more than should on r

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

for t=1:length(RPinterp_vert_x_cell)
    RPinterp_vert_cell_info(t)=length(RPinterp_vert_x_cell{t});
end

for t=1:length(RPinterp_hor_x_cell)
    RPinterp_hor_cell_info(t)=length(RPinterp_hor_x_cell{t});
end




% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% 
% figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% xlim([-(ncx+1)/2 (ncx+1)/2])
% ylim([-(ncy+1)/2 (ncy+1)/2])
% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% xlim([-(ncx+1)/2 (ncx+1)/2])
% ylim([-(ncy+1)/2 (ncy+1)/2])
% 
% % 
% 
% figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% xlim([-(ncx-1)/2 (ncx-1)/2])
% ylim([-(ncy-1)/2 (ncy-1)/2])
% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% xlim([-(ncx)/2 (ncx)/2])
% ylim([-(ncy)/2 (ncy)/2])
% 
% 
% hold on
%  scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% colorbar;





% 
% 
% 
% 
% %Vertical and horizontal lines (double interpolating)
% 
% 
% n_thetas=2*(2*(ncx+ncy));
% d_theta=2*pi/n_thetas;
% 
% thetas=0:d_theta:2*pi-d_theta;
% theta_x_half=atan(floor(ncy/2)/floor(ncx/2));
% 
% aux_count_vert=1;
% aux_count_hor=1;
% clearvars sample_points_vert sample_points_hor
% 
% RPinterp_vert_x_cell={};
% RPinterp_vert_y_cell={};
% RPinterp_hor_x_cell={};
% RPinterp_hor_y_cell={};
% 
% aux_count_interp_vert=1;
% aux_count_interp_vert_interp=1;
% aux_count_interp_vert_cell=1;
% aux_count_interp_hor=1;
% aux_count_interp_hor_interp=1;
% aux_count_interp_hor_cell=1;
% 
% 
% %para armazenar, ao invez de celulas, faz uma entrada auxiliar labelling
% %de 1 a lenght(t) (um numero para cada angulo)
% 
% for t=thetas
%     if (t<theta_x_half)|((t>=pi-theta_x_half)&(t<pi+theta_x_half))|(t>=2*pi-theta_x_half)
% %     if (t<theta_x_half)|((t>=pi-theta_x_half)&(t<=pi+theta_x_half))|(t>2*pi-theta_x_half)
%         for r=-floor(abs(2+ncx/(2*cos(t)))):1/2:floor(abs(2+ncx/(2*cos(t))))            
%             sample_points_hor(aux_count_hor,1)=r;
%             sample_points_hor(aux_count_hor,2)=t;
%             sample_points_hor(aux_count_hor,3)=r*cos(t);
%             sample_points_hor(aux_count_hor,4)=r*sin(t);            
%             sample_points_hor(aux_count_hor,5)=aux_count_interp_hor_interp;            
%             aux_count_hor=aux_count_hor+1;
%             RPinterp_hor_x(aux_count_interp_hor)=r*cos(t);
%             RPinterp_hor_y(aux_count_interp_hor)=r*sin(t);
%             aux_count_interp_hor=aux_count_interp_hor+1;
%         end 
%         RPinterp_hor_x_cell{aux_count_interp_hor_cell}=RPinterp_hor_x;
%         RPinterp_hor_y_cell{aux_count_interp_hor_cell}=RPinterp_hor_y;
%         aux_count_interp_hor_cell=aux_count_interp_hor_cell+1;
%         aux_count_interp_hor=1;
%         clearvars RPinterp_hor_x RPinterp_hor_y
%         aux_count_interp_hor_interp=aux_count_interp_hor_interp+1;
%     else
%         for r=-floor(2+abs(ncy/(2*sin(t)))):1/2:floor(2+abs(ncy/(2*sin(t))))
%             sample_points_vert(aux_count_vert,1)=r;
%             sample_points_vert(aux_count_vert,2)=t;
%             sample_points_vert(aux_count_vert,3)=r*cos(t);
%             sample_points_vert(aux_count_vert,4)=r*sin(t);                
%             sample_points_vert(aux_count_vert,5)=aux_count_interp_vert_interp;                
%             aux_count_vert=aux_count_vert+1;
%             RPinterp_vert_x(aux_count_interp_vert)=r*cos(t);
%             RPinterp_vert_y(aux_count_interp_vert)=r*sin(t);
%             aux_count_interp_vert=aux_count_interp_vert+1;            
%         end 
%         RPinterp_vert_x_cell{aux_count_interp_vert_cell}=RPinterp_vert_x;
%         RPinterp_vert_y_cell{aux_count_interp_vert_cell}=RPinterp_vert_y;
%         aux_count_interp_vert_cell=aux_count_interp_vert_cell+1;
%         aux_count_interp_vert=1;
%         clearvars RPinterp_vert_x RPinterp_vert_y
%         aux_count_interp_vert_interp=aux_count_interp_vert_interp+1;        
%     end
% end
% 
% 
% 
% for t=1:length(RPinterp_vert_x_cell)
%     RPinterp_vert_cell_info(t)=length(RPinterp_vert_x_cell{t});
% end
% 
% for t=1:length(RPinterp_hor_x_cell)
%     RPinterp_hor_cell_info(t)=length(RPinterp_hor_x_cell{t});
% end
% 


% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% 
% ncx_new = ncx;
% ncy_new = ncy;
% 
% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% xlim([-(ncx_new-1)/2 (ncx_new-1)/2])
% ylim([-(ncy_new-1)/2 (ncy_new-1)/2])
% 
% 
% figure; scatter(sample_points_hor(:,3),sample_points_hor(:,4))
% xlim([-(ncx_new-1)/2 (ncx_new-1)/2])
% ylim([-(ncy_new-1)/2 (ncy_new-1)/2])
% 
% 
% % 
% % ncx+(1-is_odd_ncx)
% % 
% % 
% % ncx_new=ncx+(1-is_odd_ncx)+4;
% % ncy_new=ncy+(1-is_odd_ncy)+4;
% 
% 
% 
% 
% figure; scatter(sample_points_vert(:,3),sample_points_vert(:,4))
% 
% xlim([-(ncx_new-1)/2 (ncx_new-1)/2])
% ylim([-(ncy_new-1)/2 (ncy_new-1)/2])
% 
% 
% hold on
%  scatter(sample_points_hor(:,3),sample_points_hor(:,4))



% [x_h,y_h] = meshgrid(-(-1+length(F_h_(1,:)))/2:(-1+length(F_h_(1,:)))/2,-(-1+length(F_h_(:,1)))/2:(-1+length(F_h_(:,1)))/2);

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


% [x_,y_] = meshgrid(-(-1+length(F_(1,:)))/2:(-1+length(F_(1,:)))/2,-(-1+length(F_(:,1)))/2:(-1+length(F_(:,1)))/2);
% RPinterp_hor = interp2(x_,y_,F_,sample_points_hor(:,3),sample_points_hor(:,4),'cubic');
% RPinterp_vert = interp2(x_,y_,F_,sample_points_vert(:,3),sample_points_vert(:,4),'cubic');


% sample_points_hor_x_test=sample_points_hor(1:17,3);
% sample_points_hor_y_test=sample_points_hor(1:17,4);
% 
% RPinterp_hor_test = interp2(x_,y_,F_',sample_points_hor_x_test,sample_points_hor_y_test,'cubic');



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




% Radon_hor_test = ifft(RPinterp_hor_test.').';








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




Radon_hor__(2,:)=sample_points_hor(:,5);

Radon_vert__(2,:)=sample_points_vert(:,5);




% Radon_hor__(2:6,:)=sample_points_hor(:,:)';
% 
% Radon_vert__(2:6,:)=sample_points_vert(:,:)';


























% 
% 
% 
% %Apply Recto-Polar interpolation
% 
% 
% 
% 
% %Vertical lines_centered
% 
% for t=1:ncx+(1-is_odd_ncx)
%     for r=1:ncy+(1-is_odd_ncy)
%         M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
%         y=min(max(round(r),1),ncy+(1-is_odd_ncy));
%         x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));
%         %         A((t-1)*ncy+r,1)=x;
%         %         A((t-1)*ncy+r,2)=y;
%         RPinterp_vert(t,r)=F_(x,y);
%     end
% end
% 
% 
% %Horizontal lines_centered
% 
% 
% for t=1:ncy+(1-is_odd_ncy)
%     for r=1:ncx+(1-is_odd_ncx)
%         M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
%         x=min(max(round(r),1),ncx+(1-is_odd_ncx));
%         y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
%         %         B((t-1)*ncx+r,1)=r;
%         %         B((t-1)*ncx+r,2)=t+M*(r-1);
%         RPinterp_hor(t,r)=F_(x,y);
%     end
% end
% 
% 
% 
% % RPinterp_hor_=[sum(RPinterp_hor')',RPinterp_hor];
% % RPinterp_vert_=[sum(RPinterp_vert')',RPinterp_vert];
% 
% % 
% % figure;imagesc(real(RPinterp_hor));colorbar;
% % figure;imagesc(real(RPinterp_vert));colorbar;
% 
% % if is_odd_ncx==0 RPinterp_vert(end,:)=[]; end
% % if is_odd_ncy==0 RPinterp_hor(end,:)=[]; end
% 
% % Radon_hor = (ifft(RPinterp_hor_.').');
% % Radon_vert = (ifft(RPinterp_vert_.').');
% % 
% 
% 
% % if is_odd_ncx==0 RPinterp_hor(:,end)=[]; end
% % if is_odd_ncy==0 RPinterp_vert(:,end)=[]; end
% 
% % % if is_odd_ncx==0 
% % %     RPinterp_hor(:,end)=[];
% % %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% % % end
% % % if is_odd_ncy==0 
% % %     RPinterp_vert(:,end)=[];
% % %     RPinterp_vert(:,1)=real(RPinterp_vert(:,1));
% % % end
% % 
% % 
% % % //store just the imaginary part,for either odd and even imaginary parts
% % % (in the even case, only removes the first colum imag
% % % (in the odd case, do not remove imaginary -> set the extra vector to sero)
% % % 
% % % RPinterp_hor_imag=imag(RPinterp_hor(:,1));
% % % 
% % % if is_odd_ncx==0 
% % %     RPinterp_hor(:,end)=[];
% % % %     RPinterp_hor_imag(end)=[];
% % % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% % %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% % % end
% % %     
% % % RPinterp_vert_imag=imag(RPinterp_vert(:,1));
% % % 
% % % if is_odd_ncy==0
% % %     RPinterp_vert(:,end)=[];
% % % %     RPinterp_vert_imag(end)=[];
% % % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% % %     RPinterp_vert(:,1)=real(RPinterp_vert(:,1));
% % % end
% % % 
% % % 
% % % 
% % 
% % 
% % % with the pseudo modulus idea, 
% % % since the phase in matlab depends on [-pi,pi], 
% % % abs(phase)-pi/2>0 
% % % is the "negativeness" condition
% % 
% % 
% % RPinterp_hor_phase=angle(RPinterp_hor(:,1));
% % 
% % if is_odd_ncx==0 
% %     RPinterp_hor(:,end)=[];
% % %     RPinterp_hor_imag(end)=[];
% % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_hor(:,1)=abs(RPinterp_hor(:,1)).*sign(real(RPinterp_hor(:,1)));
% % end
% %     
% % RPinterp_vert_phase=angle(RPinterp_vert(:,1));
% % 
% % if is_odd_ncy==0
% %     RPinterp_vert(:,end)=[];
% % %     RPinterp_vert_imag(end)=[];
% % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_vert(:,1)=abs(RPinterp_vert(:,1)).*sign(real(RPinterp_vert(:,1)));
% % end
% % 
% %     
% %     
% %     
% % % 
% % % if is_odd_ncx==0 
% % %     RPinterp_hor(:,end)=[];
% % %     RPinterp_hor_toReal=ifft(RPinterp_hor(:,1));    
% % % %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% % %     for t=1:ncy+(1-is_odd_ncy)
% % %         RPinterp_hor_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncy+(1-is_odd_ncy)))))*RPinterp_hor_toReal(t);
% % %     end
% % %     RPinterp_hor(:,1)=RPinterp_hor_toReal_;
% % % end
% % % 
% % % 
% % % if is_odd_ncy==0
% % %     RPinterp_vert(:,end)=[];
% % %     RPinterp_vert_toReal=ifft(RPinterp_vert(:,1));
% % %     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
% % %     for t=1:ncx+(1-is_odd_ncx)
% % %         RPinterp_vert_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncx+(1-is_odd_ncx)))))*RPinterp_vert_toReal(t);
% % %     end
% % %     RPinterp_vert(:,1)=RPinterp_vert_toReal_;
% % % end
% % % 
% 
% 
% 
% 
% 
% % Obtain Radon by inverse fft of the the RP interpolations
% 
% 
% Radon_hor = (ifft(RPinterp_hor.').');
% Radon_vert = (ifft(RPinterp_vert.').');
% 
% 
% 
% 
% 
% 
% % convention to the frequency center fft
% 
% 
% % 
% % for r=1:ncx
% %     for t=1:ncy+(1-is_odd_ncy)
% %         Radon_hor__(t,r)=(exp(-1i*pi*(r-1)*(1-is_odd_ncx/ncx)))*Radon_hor(t,r);
% %     end
% % end
% % 
% % for t=1:ncx+(1-is_odd_ncx)
% %     for r=1:ncy
% %         Radon_vert__(t,r)=(exp(-1i*pi*(r-1)*(1-is_odd_ncy/ncy)))*Radon_vert(t,r);
% %     end
% % end
% 
% 
% for r=1:ncx+(1-is_odd_ncx)
%     for t=1:ncy+(1-is_odd_ncy)
%         Radon_hor__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncx+(1-is_odd_ncx)))))*Radon_hor(t,r);
%     end
% end
% 
% for t=1:ncx+(1-is_odd_ncx)
%     for r=1:ncy+(1-is_odd_ncy)
%         Radon_vert__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncy+(1-is_odd_ncy)))))*Radon_vert(t,r);
%     end
% end
% 
% % % RPinterp_hor_imag=imag(RPinterp_hor(:,1));    
% % % RPinterp_hor_imag=imag(RPinterp_vert(:,1));
% % % 
% % % Radon_hor__(:,end+1)=RPinterp_hor_imag;
% % % Radon_vert__(:,end+1)=RPinterp_vert_imag;
% % 
% % 
% % Radon_hor__(:,end+1)=RPinterp_hor_phase;
% % Radon_vert__(:,end+1)=RPinterp_vert_phase;
% 
% 
% % Radon_hor__=Radon_hor;
% % Radon_vert__=Radon_vert;
% % 
% % figure;imagesc(real(Radon_hor__)); colorbar;
% % xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% % ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% % 
% % 
% % figure;imagesc(real(Radon_vert__)); colorbar;
% % xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% % ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Redundant augmentation in "odd" case (to make the dimension implicit)
% 
% 
% if is_odd_ncy==1 Radon_hor__(end+1,:)=Radon_hor__(1,:);end
% if is_odd_ncx==1 Radon_vert__(end+1,:)=Radon_vert__(1,:);end



end
