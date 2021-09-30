function [  ] = Ridgelet3d_interp_visual_dev3( Radon_X__,Radon_Y__,Radon_Z__,interp_X_info,interp_Y_info,interp_Z_info)
% 
Radon_X__ = Radon_Z__;
Radon_Y__ = Radon_Z__;
% Radon_Z__ = Radon_Z__;


% Radon_X__ = Radon_X_;
% Radon_Y__ = Radon_Y_;
% Radon_Z__ = Radon_Z_;



% h = volshow(X,'Renderer','MaximumIntensityProjection');

% sample_points_hor = Radon_hor__(2,:);

ncX_r = max(interp_X_info(:));
ncX_tz = length(interp_X_info(:,1));
ncX_ty = length(interp_X_info(1,:));


ncY_r = max(interp_Y_info(:));
ncY_tx = length(interp_Y_info(:,1));
ncY_tz = length(interp_Y_info(1,:));


ncZ_r = max(interp_Z_info(:));
ncZ_tx = length(interp_Z_info(:,1));
ncZ_ty = length(interp_Z_info(1,:));



Radon_X__vis = zeros(ncX_tz,ncX_ty,ncX_r);
Radon_Y__vis = zeros(ncY_tx,ncY_tz,ncY_r);
Radon_Z__vis = zeros(ncZ_tx,ncZ_ty,ncZ_r);


count_aux=1;

for t1 = 1:ncX_tz
    for t2 = 1:ncX_ty
        for r=1:interp_X_info(t1,t2)
            r_loc=(ncX_r-interp_X_info(t1,t2))/2+r;
            Radon_X__vis(t1,t2,r_loc)=Radon_X__(1,count_aux);
            count_aux = count_aux +1;
            
        end
    end
end

volshow(real(Radon_X__vis),'Renderer','MaximumIntensityProjection');

% figure; imagesc(real(Radon_X__vis(:,:,10)));colorbar;


count_aux=1;

for t1 = 1:ncY_tx
    for t2 = 1:ncY_tz
        for r=1:interp_Y_info(t1,t2)
            r_loc=(ncY_r-interp_Y_info(t1,t2))/2+r;
            Radon_Y__vis(t1,t2,r_loc)=Radon_Y__(1,count_aux);
            count_aux = count_aux +1;
            
        end
    end
end

volshow(real(Radon_Y__vis),'Renderer','MaximumIntensityProjection');

% figure; imagesc(real(Radon_Y__vis(:,:,10)));colorbar;



count_aux=1;

for t1 = 1:ncZ_tx
    for t2 = 1:ncZ_ty
        for r=1:interp_Z_info(t1,t2)
            r_loc=(ncZ_r-interp_Z_info(t1,t2))/2+r;
            Radon_Z__vis(t1,t2,r_loc)=Radon_Z__(1,count_aux);
            count_aux = count_aux +1;            
        end
    end
end

volshow(real(Radon_Z__vis),'Renderer','MaximumIntensityProjection');

% figure; imagesc(real(Radon_Z__vis(:,:,10)));colorbar;



% % ncr = length(Radon_vert_aux__);
% t_loc_count=1;
% 
% for r=1:ncr_vert
%     for t = 1:RPinterp_vert_cell_info(r)
% 
%     t_loc=(nct_vert-RPinterp_vert_cell_info(r))/2+t;
%     Radon_vert__vis(t_loc,r)=Radon_vert__(1,t_loc_count);    
%     t_loc_count = t_loc_count +1;
%     
%     end
% end
% % Radon_vert__=[Radon_vert__,Radon_vert_aux__];
% 
% % figure; imagesc(Radon_hor__vis);colorbar;
% figure; imagesc(real(Radon_vert__vis));colorbar;
% % figure; imagesc(imag(Radon_vert__vis));colorbar;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Radon_hor__vis = zeros(nct_hor,ncr_hor);
% 
% t_loc_count=1;
% 
% for r=1:ncr_hor
%     for t = 1:RPinterp_hor_cell_info(r)
% 
%     t_loc=(nct_hor-RPinterp_hor_cell_info(r))/2+t;
%     Radon_hor__vis(t_loc,r)=Radon_hor__(1,t_loc_count);    
%     t_loc_count = t_loc_count +1;
%     
%     end
% end
% 
% figure; imagesc(real(Radon_hor__vis));colorbar;
% % figure; imagesc(imag(Radon_hor__vis));colorbar;
% 
% 
% 
% % figure; imagesc(imag(Radon_hor__vis));colorbar;
% 
% 
% 
% 
% % 
% % 
% % aux_RPinterp_hor_count=1;
% % Radon_hor_aux=[];
% % Radon_hor=[];
% % for i=1:length(sample_points_hor(:,1))
% %     if (sample_points_hor(i,5)==aux_RPinterp_hor_count)
% %         Radon_hor_aux(end+1)=RPinterp_hor(i);
% %     else
% %        
% % %         Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
% % %         display(aux_RPinterp_hor_count);
% %         
% %         Radon_hor_aux=[];  
% %         
% %         aux_RPinterp_hor_count=aux_RPinterp_hor_count+1;
% %         
% %         Radon_hor_aux(end+1)=RPinterp_hor(i);
% % 
% %     end
% % end
% % 
% % Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
% % % display(aux_RPinterp_hor_count);
% % 
% 
% % figure; imagesc(Radon_hor);colourbar;

end

