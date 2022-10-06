% function [Rid_Z] = Ridgelet3d_fromRadon_dev3(Radon_Z__,interp_Z_info,lev_3drig)
% function [Ridgelet_Z__] = Ridgelet3d_fromRadon_dev3(Radon_Z_,sample_points_Z_id,interp_Z_info,lev_3drig)
function [ Ridgelet_Z__ ] = Ridgelet3d_fromRadon_dev3( Radon_Z__,interp_Z_info,lev_3drig)

% new version

ncZ_tx = length(interp_Z_info(:,1));
ncZ_ty = length(interp_Z_info(1,:));

count_aux=1;
Radon_Z_aux__=[];
Ridgelet_Z__=[];

for t1 = 1:ncZ_tx
    for t2 = 1:ncZ_ty
        
        for r=1:interp_Z_info(t1,t2)
%             r_loc=(ncZ_r-interp_Z_info(t1,t2))/2+r;
            Radon_Z_aux__(end+1)=Radon_Z__(1,count_aux);
            count_aux = count_aux +1;            
        end
        
        n_levels=floor(log2(length(Radon_Z_aux__)));
        
        [dc_dwt,levels] = wavedec(real(Radon_Z_aux__),n_levels,'db1');
        [Ridgelet_Z_aux_] = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);        
%         Ridgelet_Z_aux_ = [D,D];
%         Ridgelet_Z_aux_(end) = [];
        
%         Radon_Z_aux2__ = Radon_Z_aux__;
        Ridgelet_Z__=[Ridgelet_Z__,Ridgelet_Z_aux_];
        Radon_Z_aux__=[];
        
    end
end



% % old version 3
% 
% 
% n_min = min(interp_Z_info(:));
% ha= floor((n_min-1)/4); %half amplitude for allowable noNan
% 
% 
% % Radon_Z_ = Radon_Z__;
% 
% n_levels=floor(log2(length(Radon_Z_)));
% [dc_dwt,levels] = wavedec(real(Radon_Z_),n_levels,'db1');
% D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
% Ridgelet_Z = D;
% 
% 
% 
% aux_interp_Z_count__=1;
% Ridgelet_Z_aux__=[];
% Ridgelet_Z__=[];
% for i=1:length(sample_points_Z_id)
%     if (sample_points_Z_id(i)==aux_interp_Z_count__)
%         Ridgelet_Z_aux__(end+1)=Ridgelet_Z(i);
%     else
%         
% %         Ridgelet_Z_aux__(1)= NaN;
% %         Ridgelet_Z_aux__(end)= NaN;
%         Ridgelet_Z_aux__(1:ha) = NaN;
%         Ridgelet_Z_aux__(end-ha:end) = NaN;
% 
%         Ridgelet_Z__=[Ridgelet_Z__,Ridgelet_Z_aux__];
%         
%         Ridgelet_Z_aux__=[];  
%         aux_interp_Z_count__=aux_interp_Z_count__+1;
%         Ridgelet_Z_aux__(end+1)=Ridgelet_Z(i);
%     end
% end
% 
% % Ridgelet_Z_aux__(1)= NaN;
% % Ridgelet_Z_aux__(end)= NaN;
% Ridgelet_Z_aux__(1:ha) = NaN;
% Ridgelet_Z_aux__(end-ha:end) = NaN;
%         
%         
% Ridgelet_Z__=[Ridgelet_Z__,Ridgelet_Z_aux__];






%Old Version 2





% % Radon_Z = 1;
% % 
% % sample_points_Z = 1;
% 
% % Radon_Z_ = 
% 
% aux_interp_Z_count__=1;
% Radon_Z_aux__=[];
% Ridgelet_Z__=[];
% for i=1:length(sample_points_Z_id)
%     if (sample_points_Z_id(i)==aux_interp_Z_count__)
%         Radon_Z_aux__(end+1)=Radon_Z_(i);
%     else
% %         ncr = length(Radon_Z_aux__);
% %         for r=1:ncr
% %             Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
% %         end
%         
%         n_levels=floor(log2(length(Radon_Z_aux__)));
%         [dc_dwt,levels] = wavedec(real(Radon_Z_aux__),n_levels,'db1');
%         D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
%         Ridgelet_Z_aux2__ = D;
%         
% %         [c,l] = wavedec(real(Radon_Z_aux__),lev_3drig,'db1');
% %         [D] = detcoef(c,l,[1]);        
% %         Ridgelet_Z_aux2__ = [D,D];
% %         Ridgelet_Z_aux2__(end) = [];
%         
% %         Radon_Z_aux2__ = Radon_Z_aux__;
%         Ridgelet_Z__=[Ridgelet_Z__,Ridgelet_Z_aux2__];
%         
%         Radon_Z_aux__=[];  
%         aux_interp_Z_count__=aux_interp_Z_count__+1;
%         Radon_Z_aux__(end+1)=Radon_Z_(i);
%     end
% end
% % ncr = length(Radon_Z_aux__);
% % for r=1:ncr
% %     Radon_Z_aux__(r)=(exp(-1i*pi*(r-1)*(1-1/ncr)))*Radon_Z_aux__(r);
% % end
% 
%         n_levels=floor(log2(length(Radon_Z_aux__)));
%         [dc_dwt,levels] = wavedec(real(Radon_Z_aux__),n_levels,'db1');
%         D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
%         Ridgelet_Z_aux2__ = D;
% 
% %         [c,l] = wavedec(real(Radon_Z_aux__),lev_3drig,'db1');
% %         [D] = detcoef(c,l,[1]);
% %         Ridgelet_Z_aux2__ = [D,D];
% %         Ridgelet_Z_aux2__(end) = [];
% 
% % Radon_Z_aux2__ = Radon_Z_aux__;
% 
% Ridgelet_Z__=[Ridgelet_Z__,Ridgelet_Z_aux2__];







%old version

% 
% 
% % Radon_Z_(1,end+1) = 0;
% % 
% % Radon_Z__ = Radon_Z_;
% 
% % lev_3drig = 1;
% 
% Radon_Z__(1,end+1) = 0;
% 
% % ncZ_r = max(interp_Z_info(:));
% ncZ_tx = length(interp_Z_info(:,1));
% ncZ_ty = length(interp_Z_info(1,:));
% 
% % Radon_Z__vis = zeros(ncZ_tx,ncZ_ty,ncZ_r);
% Rid_Z = [];
% 
% count_aux=1;
% 
% for t1 = 1:ncZ_tx
%     for t2 = 1:ncZ_ty
%         Rid_Z_aux = Radon_Z__(1,count_aux:count_aux+interp_Z_info(t1,t2));
%         
%         n_levels=floor(log2(length(Rid_Z_aux)));
%         [dc_dwt,levels] = wavedec(real(Rid_Z_aux),n_levels,'db1');
%         D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
%         
%         Rid_Z_aux2 = D;
%         
%         Rid_Z(count_aux:count_aux+interp_Z_info(t1,t2)) = Rid_Z_aux2;
%         count_aux = count_aux + interp_Z_info(t1,t2);
%        
%     end
% end
% 
% Rid_Z(end) = [];
% 


end

