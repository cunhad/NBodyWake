function [ReducRadon_Z__] = RemoveLowVolRad_dev3(Radon_Z__,interp_Z_info,remov_fact)


% Radon_Z__ = Radon_Z_ ;
% % interp_Z_info = 
% remov_fact = 1/2;


ncZ_tx = length(interp_Z_info(:,1));
ncZ_ty = length(interp_Z_info(1,:));

ncZ = min(interp_Z_info(:));

count_aux=1;
Radon_Z_aux__=[];
ReducRadon_Z__=[];

for t1 = 1:ncZ_tx
    for t2 = 1:ncZ_ty
        
        for r=1:interp_Z_info(t1,t2)
%             r_loc=(ncZ_r-interp_Z_info(t1,t2))/2+r;
            Radon_Z_aux__(end+1)=Radon_Z__(1,count_aux);
            count_aux = count_aux +1;            
        end
        
%         n_levels=floor(log2(length(Radon_Z_aux__)));
%         
%         [dc_dwt,levels] = wavedec(real(Radon_Z_aux__),n_levels,'db1');
%         [Ridgelet_Z_aux_] = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);        
% %         Ridgelet_Z_aux_ = [D,D];
% %         Ridgelet_Z_aux_(end) = [];

        leng = length(Radon_Z_aux__);
        remov_low_ind = ceil(leng/2 - remov_fact*ncZ/2); 
        remov_hig_ind = floor(leng/2 + remov_fact*ncZ/2);
        
        Radon_Z_aux__(remov_hig_ind:end) = [];
        Radon_Z_aux__(1:remov_low_ind) = [];
        
%         Radon_Z_aux2__ = Radon_Z_aux__;
        ReducRadon_Z__=[ReducRadon_Z__,Radon_Z_aux__];
        Radon_Z_aux__=[];
        
    end
end



end

