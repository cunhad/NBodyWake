function [ Radon_hor__,Radon_vert__] = Ridgelet2d_RP_crude_forwards_dev1(X)

% using real ifft fix, 

% sometimes it is better to pass, not the imag part, but just the phase (mwith signed modulo)





% implicit variable

[ncx,ncy]=size(X);

is_odd_ncx=mod(ncx,2);
is_odd_ncy=mod(ncy,2);


%Apply FFT2

for i=1:ncx
    for j=1:ncy
        X_(i,j)=(exp(1i*pi*(i-1)*(1-is_odd_ncx/ncx)))*(exp(1i*pi*(j-1)*(1-is_odd_ncy/ncy)))*X(i,j);
    end
end


F_ = fft(fft(X_).').';

if is_odd_ncy==0 F_(:,end+1)=F_(:,1);end
if is_odd_ncx==0 F_(end+1,:)=F_(1,:);end





%Apply Recto-Polar interpolation




%Vertical lines_centered

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        M=(ncx+(1-is_odd_ncx)-(2*t)+1)/(ncy-1+(1-is_odd_ncy));
        y=min(max(round(r),1),ncy+(1-is_odd_ncy));
        x=min(max(fix(t+M*(r-1)-ceil((ncx+1)/2))+ceil((ncx+1)/2),1),ncx+(1-is_odd_ncx));
        %         A((t-1)*ncy+r,1)=x;
        %         A((t-1)*ncy+r,2)=y;
        RPinterp_vert(t,r)=F_(x,y);
    end
end


%Horizontal lines_centered


for t=1:ncy+(1-is_odd_ncy)
    for r=1:ncx+(1-is_odd_ncx)
        M=(ncy+(1-is_odd_ncy)-(2*t)+1)/(ncx-1+(1-is_odd_ncx));
        x=min(max(round(r),1),ncx+(1-is_odd_ncx));
        y=min(max(fix(t+M*(r-1)-ceil((ncy+1)/2))+ceil((ncy+1)/2),1),ncy+(1-is_odd_ncy));
        %         B((t-1)*ncx+r,1)=r;
        %         B((t-1)*ncx+r,2)=t+M*(r-1);
        RPinterp_hor(t,r)=F_(x,y);
    end
end



% RPinterp_hor_=[sum(RPinterp_hor')',RPinterp_hor];
% RPinterp_vert_=[sum(RPinterp_vert')',RPinterp_vert];

% 
% figure;imagesc(real(RPinterp_hor));colorbar;
% figure;imagesc(real(RPinterp_vert));colorbar;

% if is_odd_ncx==0 RPinterp_vert(end,:)=[]; end
% if is_odd_ncy==0 RPinterp_hor(end,:)=[]; end

% Radon_hor = (ifft(RPinterp_hor_.').');
% Radon_vert = (ifft(RPinterp_vert_.').');
% 


% if is_odd_ncx==0 RPinterp_hor(:,end)=[]; end
% if is_odd_ncy==0 RPinterp_vert(:,end)=[]; end

% % if is_odd_ncx==0 
% %     RPinterp_hor(:,end)=[];
% %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% % end
% % if is_odd_ncy==0 
% %     RPinterp_vert(:,end)=[];
% %     RPinterp_vert(:,1)=real(RPinterp_vert(:,1));
% % end
% 
% 
% % //store just the imaginary part,for either odd and even imaginary parts
% % (in the even case, only removes the first colum imag
% % (in the odd case, do not remove imaginary -> set the extra vector to sero)
% % 
% % RPinterp_hor_imag=imag(RPinterp_hor(:,1));
% % 
% % if is_odd_ncx==0 
% %     RPinterp_hor(:,end)=[];
% % %     RPinterp_hor_imag(end)=[];
% % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% % end
% %     
% % RPinterp_vert_imag=imag(RPinterp_vert(:,1));
% % 
% % if is_odd_ncy==0
% %     RPinterp_vert(:,end)=[];
% % %     RPinterp_vert_imag(end)=[];
% % %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
% %     RPinterp_vert(:,1)=real(RPinterp_vert(:,1));
% % end
% % 
% % 
% % 
% 
% 
% % with the pseudo modulus idea, 
% % since the phase in matlab depends on [-pi,pi], 
% % abs(phase)-pi/2>0 
% % is the "negativeness" condition
% 
% 
% RPinterp_hor_phase=angle(RPinterp_hor(:,1));
% 
% if is_odd_ncx==0 
%     RPinterp_hor(:,end)=[];
% %     RPinterp_hor_imag(end)=[];
% %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
%     RPinterp_hor(:,1)=abs(RPinterp_hor(:,1)).*sign(real(RPinterp_hor(:,1)));
% end
%     
% RPinterp_vert_phase=angle(RPinterp_vert(:,1));
% 
% if is_odd_ncy==0
%     RPinterp_vert(:,end)=[];
% %     RPinterp_vert_imag(end)=[];
% %     RPinterp_hor_imag=RPinterp_hor(:,1)-real(RPinterp_hor(:,1));
%     RPinterp_vert(:,1)=abs(RPinterp_vert(:,1)).*sign(real(RPinterp_vert(:,1)));
% end
% 
%     
%     
%     
% % 
% % if is_odd_ncx==0 
% %     RPinterp_hor(:,end)=[];
% %     RPinterp_hor_toReal=ifft(RPinterp_hor(:,1));    
% % %     RPinterp_hor(:,1)=real(RPinterp_hor(:,1));
% %     for t=1:ncy+(1-is_odd_ncy)
% %         RPinterp_hor_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncy+(1-is_odd_ncy)))))*RPinterp_hor_toReal(t);
% %     end
% %     RPinterp_hor(:,1)=RPinterp_hor_toReal_;
% % end
% % 
% % 
% % if is_odd_ncy==0
% %     RPinterp_vert(:,end)=[];
% %     RPinterp_vert_toReal=ifft(RPinterp_vert(:,1));
% %     %     RPinterp_vert(:,1)=fft(circshift(RPinterp_vert(:,1),65));
% %     for t=1:ncx+(1-is_odd_ncx)
% %         RPinterp_vert_toReal_(t)=(exp(-1i*pi*(t-1)*(1-1/(ncx+(1-is_odd_ncx)))))*RPinterp_vert_toReal(t);
% %     end
% %     RPinterp_vert(:,1)=RPinterp_vert_toReal_;
% % end
% % 


Radon_hor = (ifft(RPinterp_hor.').');
Radon_vert = (ifft(RPinterp_vert.').');

% 
% for r=1:ncx
%     for t=1:ncy+(1-is_odd_ncy)
%         Radon_hor__(t,r)=(exp(-1i*pi*(r-1)*(1-is_odd_ncx/ncx)))*Radon_hor(t,r);
%     end
% end
% 
% for t=1:ncx+(1-is_odd_ncx)
%     for r=1:ncy
%         Radon_vert__(t,r)=(exp(-1i*pi*(r-1)*(1-is_odd_ncy/ncy)))*Radon_vert(t,r);
%     end
% end


for r=1:ncx+(1-is_odd_ncx)
    for t=1:ncy+(1-is_odd_ncy)
        Radon_hor__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncx+(1-is_odd_ncx)))))*Radon_hor(t,r);
    end
end

for t=1:ncx+(1-is_odd_ncx)
    for r=1:ncy+(1-is_odd_ncy)
        Radon_vert__(t,r)=(exp(-1i*pi*(r-1)*(1-1/(ncy+(1-is_odd_ncy)))))*Radon_vert(t,r);
    end
end

% % RPinterp_hor_imag=imag(RPinterp_hor(:,1));    
% % RPinterp_hor_imag=imag(RPinterp_vert(:,1));
% % 
% % Radon_hor__(:,end+1)=RPinterp_hor_imag;
% % Radon_vert__(:,end+1)=RPinterp_vert_imag;
% 
% 
% Radon_hor__(:,end+1)=RPinterp_hor_phase;
% Radon_vert__(:,end+1)=RPinterp_vert_phase;


% Radon_hor__=Radon_hor;
% Radon_vert__=Radon_vert;
% 
% figure;imagesc(real(Radon_hor__)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);
% 
% 
% figure;imagesc(real(Radon_vert__)); colorbar;
% xlabel('displacement', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('angle', 'interpreter', 'latex', 'fontsize', 20);


% 

if is_odd_ncy==1 Radon_hor__(end+1,:)=Radon_hor__(1,:);end
if is_odd_ncx==1 Radon_vert__(end+1,:)=Radon_vert__(1,:);end



end
