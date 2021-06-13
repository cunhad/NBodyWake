function [ output_args ] = Ridgelet2d_RP_crude_visual_dev3( Radon_hor__,Radon_vert__,RPinterp_hor_cell_info,RPinterp_vert_cell_info,ncx,ncy )



% sample_points_hor = Radon_hor__(2,:);



nct_hor = max(RPinterp_hor_cell_info);
ncr_hor = length(RPinterp_hor_cell_info);


nct_vert = max(RPinterp_vert_cell_info);
ncr_vert = length(RPinterp_vert_cell_info);






Radon_vert__vis = zeros(nct_vert,ncr_vert);

% ncr = length(Radon_vert_aux__);
t_loc_count=1;

for r=1:ncr_vert
    for t = 1:RPinterp_vert_cell_info(r)

    t_loc=(nct_vert-RPinterp_vert_cell_info(r))/2+t;
    Radon_vert__vis(t_loc,r)=Radon_vert__(1,t_loc_count);    
    t_loc_count = t_loc_count +1;
    
    end
end
% Radon_vert__=[Radon_vert__,Radon_vert_aux__];

% figure; imagesc(Radon_hor__vis);colorbar;
figure; imagesc(real(Radon_vert__vis));colorbar;
% figure; imagesc(imag(Radon_vert__vis));colorbar;





















Radon_hor__vis = zeros(nct_hor,ncr_hor);

t_loc_count=1;

for r=1:ncr_hor
    for t = 1:RPinterp_hor_cell_info(r)

    t_loc=(nct_hor-RPinterp_hor_cell_info(r))/2+t;
    Radon_hor__vis(t_loc,r)=Radon_hor__(1,t_loc_count);    
    t_loc_count = t_loc_count +1;
    
    end
end

figure; imagesc(real(Radon_hor__vis));colorbar;
% figure; imagesc(imag(Radon_hor__vis));colorbar;



% figure; imagesc(imag(Radon_hor__vis));colorbar;




% 
% 
% aux_RPinterp_hor_count=1;
% Radon_hor_aux=[];
% Radon_hor=[];
% for i=1:length(sample_points_hor(:,1))
%     if (sample_points_hor(i,5)==aux_RPinterp_hor_count)
%         Radon_hor_aux(end+1)=RPinterp_hor(i);
%     else
%        
% %         Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
% %         display(aux_RPinterp_hor_count);
%         
%         Radon_hor_aux=[];  
%         
%         aux_RPinterp_hor_count=aux_RPinterp_hor_count+1;
%         
%         Radon_hor_aux(end+1)=RPinterp_hor(i);
% 
%     end
% end
% 
% Radon_hor=[Radon_hor,(ifft(Radon_hor_aux.').')];
% % display(aux_RPinterp_hor_count);
% 

% figure; imagesc(Radon_hor);colourbar;

end

