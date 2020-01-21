function [ anali,signal,equator_phi ] = slices_curv_2a3d_locangreal_all_groups( )

%(example)  [ anali] = slices_curv_2a3d_all('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','/home/asus/Dropbox/extras/storage/graham/small_res/anali_hpx/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,2,2 );

%clearvars

cd('../../preprocessing');

%for this resolution, one in eahc 40 orientations (sqrt(384)) will show the wake, so we
%take 100 aprox = 3840/40

% filename='_2dproj_z3_data_sl';
% nc=1024;
% trsh=20;
% cut=1;
% lev=2;
% sigma = 5;
% slices=32;
% anal_lev=2;
NSIDE=4;
z_glob=str2num('3.000')
z='3.000';

N_angles=12*NSIDE*NSIDE/2;


specs_path_list_nowake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaraCannyl1lr1_anali_hpx_Sa4t1/16Mpc_2048c_1024p_zi63_nowakem')
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

% specs_path_list_wake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaraCannyl1lr1_analiFhm_hpx_Sa2t1/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m')
specs_path_list_wake=strcat('/home/asus/Dropbox/extras/storage/graham/ht/data_cps128_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaraCannyl1lr1_anali_hpx_Sa4t1/16Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m')
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');
% sample_list_wake=sort_nat(sample_list_wake)

for w_nw=1:2
    % for w_nw=2
    
    if w_nw==1
        [~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info('/home/asus/Dropbox/extras/storage/graham/ht/','4Mpc_2048c_1024p_zi63_nowakem','/sample3001/');
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
        ch='_7';
        coul='b';
    else
        [~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info('/home/asus/Dropbox/extras/storage/graham/ht/','4Mpc_2048c_1024p_zi63_wakeGmu2t10m8zi10m','/sample3001/half_lin_cutoff_half_tot_pert_nvpw_v0p6/');
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
        ch='_4';
        coul='r';
    end
    
    if w_nw==1
             signal_nw=[];
    else            
            signal_w=[];
    end
        
    
    
    for sample = 1:length(sample_list)
        
        path_in=strcat(specs_path_list,'/',string(sample_list(sample)),'/data_anali/1lf_0.5rf/NSIDE_',num2str(NSIDE),'/slices_curv_2a3d/hpx/molvd/')
        
        filename=strcat(path_in,'/_',num2str(find(str2num(char(redshift_list))==z_glob)),'_molvp_hpx_slice_cuvr_2a3d_z',z,'_data.txt')
        
        info = dlmread(filename);
        
        if w_nw==1
            signal_nw=[signal_nw info(1:N_angles,3)'];
            signal_sample_nw(sample,:)=info(1:N_angles,3)';
        else
            signal_w=[signal_w info(1:N_angles,3)'];
            signal_sample_w(sample,:)=info(1:N_angles,3)';
        end
        
    end
    
    

    
end

mean_nw=mean(signal_sample_nw(:))
std_nw=std(signal_sample_nw(:))

mean_w=mean(signal_sample_w(:))
std_w=std(signal_sample_w(:))

sorted_signal_sample_nw=sort((signal_sample_nw),2);
sorted_signal_sample_w=sort((signal_sample_w),2);

stn_nw=sort((signal_sample_nw-mean_nw)/std_nw,2);
stn_w=sort((signal_sample_w-mean_nw)/std_nw,2);


sorted_signal_sample_nw=sort((signal_sample_nw),2);
sorted_signal_sample_w=sort((signal_sample_w),2);

thresh=max(signal_sample_nw(:))
% thresh=532;
outlier_w_count=sum(signal_sample_w>thresh,2)
sum(sum(signal_sample_w>thresh,2))



figure;
h1 = histogram(signal_sample_nw(:),'BinWidth',10);
hold on
h2 = histogram(signal_sample_w(:),'BinWidth',10);
xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
set(gca, 'YScale', 'log')

cd('../wake_detection/slices_curv_2d/');



stn_nw_count=sum(stn_nw>4,2);
stn_w_count=sum(stn_w>4,2);



stn_w_self=sort((signal_sample_w-mean_w)/std_w,2);

[h,p,ksstat,cv]  =  kstest(stn_nw(:),'Alpha',0.01)

[h,p] = lillietest((signal_sample_nw(6,:)))

x=signal_sample_nw([4,5],:);
[h,p] = lillietest(x(:))

x=mean(signal_sample_nw,2);
[h,p] = lillietest(x(:))

sorted_signal_sample_nw_tot=sort((signal_sample_nw(:)));
sorted_signal_sample_w_tot=sort((signal_sample_w(:)));


% cd('../wake_detection/slices_curv_2d/');




% a=max(signal_sample_nw')'
% b=max(signal_sample_w')'
% s=std(a)
% m=mean(a)
% (a-m)/s
% (b-m)/s


% thresh=590;
thresh=532;
outlier_w_count=sum(signal_sample_w>thresh,2);
sum_w=sum(outlier_w_count)
sum_w_=floor(sum(outlier_w_count)/24)
prob_w=sum(outlier_w_count)/(3840*24)
outlier_nw_count=sum(signal_sample_nw>thresh,2);
sum_nw=sum(outlier_nw_count)
prob_nw=sum(outlier_nw_count)/3840
sum(outlier_w_count)/(24*sum(outlier_nw_count))


% N=10000;
N=10000;
number_nowake=round(N*prob_nw)
number_wake=round(N*prob_w)

s_value=0;
% for i=0:number_wake
for i=0:20
times=i;
value(i+1)=((1-prob_nw)^(N-times))*((prob_nw)^times)*nchoosek(N,times);
s_value=s_value+value(i+1);
sum_value(i+1)=s_value;
end
sum_value_notrange_nw=1-sum_value;

s_value=0;
% for i=0:number_wake
for i=0:20
times=i;
value(i+1)=((1-prob_w)^(N-times))*((prob_w)^times)*nchoosek(N,times);
s_value=s_value+value(i+1);
sum_value(i+1)=s_value;
end
sum_value_notrange_w=1-sum_value;

sum_value_notrange_nw(14)
sum_value_notrange_w(14)

% 
% thresh=4.5;
% stn_w_count=sum(stn_w>thresh,2);
% sum_w=sum(stn_w_count)
% sum_w_=floor(sum(stn_w_count)/32)
% prob_w=sum(stn_w_count)/(3840*32)
% stn_nw_count=sum(stn_nw>thresh,2);
% sum_nw=sum(stn_nw_count)
% prob_nw=sum(stn_nw_count)/3840
% sum(stn_w_count)/(32*sum(stn_nw_count))
% 



% times=3;
% ((1-2.6*10^-4)^(10^4-times))*((2.6*10^-4)^times)*nchoosek(10^4,times)



% s_value=0;
% for i=0:20
% times=i;
% value(i+1)=((1-2.6*10^-4)^(10^4-times))*((2.6*10^-4)^times)*nchoosek(10^4,times);
% s_value=s_value+value(i+1);
% sum_value(i+1)=s_value;
% end
% sum_value_notrange=1-sum_value;

% 
% all_signal_sample_nw=signal_sample_nw(:);
% all_signal_sample_nw_sorted=sort(all_signal_sample_nw,'descend');
% % rel_frac=ceil(length(all_signal_sample_nw_sorted)/64)
% rel_frac=ceil(length(all_signal_sample_nw_sorted)/8)
% relevand_fraction_signal_sample_nw=all_signal_sample_nw_sorted(1:rel_frac);
% 


% figure;
% h1 = histogram(signal_sample_nw(:),'BinWidth',10);
% hold on
% h2 = histogram(signal_sample_w(:),'BinWidth',10);
% xlabel('$S$ value', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('histogram', 'interpreter', 'latex', 'fontsize', 20);
% legend('G\mu=0','G\mu=1 \times 10^{-7}','location','northeast')
% set(gca, 'YScale', 'log')

% 
% all_signal_sample_w=signal_sample_w(:);
% all_signal_sample_w_sorted=sort(all_signal_sample_w,'descend');
% % rel_frac=ceil(length(all_signal_sample_w_sorted)/64)
% rel_frac=ceil(length(all_signal_sample_w_sorted)/(8*128))
% relevand_fraction_signal_sample_w=all_signal_sample_w_sorted(1:rel_frac);
% 
% max_wake_slices=relevand_fraction_signal_sample_w;
% max_nowake_slices=relevand_fraction_signal_sample_nw;
% mean_wake=mean(max_wake_slices)
% mean_nowake=mean(max_nowake_slices)
% std_wake=std(max_wake_slices,1)
% std_nowake=std(max_nowake_slices,1)
% stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
% stn=mean_stn-std_stn
% significance=abs(mean_wake-mean_nowake)/(std_wake+std_nowake)


% signal_sample_sorted_nw=sort(signal_sample_nw,2);
% signal_sample_sorted_w=sort(signal_sample_w,2);
% 
% % find(signal_sample_w(4,:)==582.62)   %=95,sigma=4.9050499
% find(signal_sample_w(6,:)==696.74)   %=338,sigma=8.0051
% find(signal_sample_w(6,:)==623.64)   %=362,sigma=6.0174
% find(signal_sample_w(4,:)==631.6700)   %=368,sigma=6.2385
% find(signal_sample_w(9,:)==628.4300)   %=305,sigma=6.1477
% % find(signal_sample_w(9,:)==875.8000)   %=338,sigma=12.8738


% find(signal_sample_nw(4,:)==592.39)   %=95,sigma=5.1677
% find(signal_sample_nw(6,:)==557.0900)   %=28,sigma=4.2079
% find(signal_sample_nw(3,:)==546.9300)   %=169,sigma=3.9
% find(signal_sample_nw(5,:)==552.5200)   %=170,sigma=4.0836






% 
% kstest2(signal_sample_nw(1:9),signal_sample_nw(10,:))
%  kstest2(signal_sample_nw([1:8,10]),signal_sample_nw(9,:))
%  kstest2(signal_sample_nw([1:7,9:10]),signal_sample_nw(8,:))
%  kstest2(signal_sample_nw([1:6,8:10]),signal_sample_nw(7,:))
%  kstest2(signal_sample_nw([1:5,7:10]),signal_sample_nw(6,:))
%  kstest2(signal_sample_nw([1:4,6:10]),signal_sample_nw(5,:))
%  kstest2(signal_sample_nw([1:3,5:10]),signal_sample_nw(4,:))
%  kstest2(signal_sample_nw([1:2,4:10]),signal_sample_nw(3,:))
%  kstest2(signal_sample_nw([1,3:10]),signal_sample_nw(2,:))
%  kstest2(signal_sample_nw(2:10),signal_sample_nw(1,:))


% 
% fig=figure;     
%     histogram(signal_nw(:));
%        
%     sorted_signal_nw=sort(signal_nw(:),'descend'); 
%     truncated_sorted_signal_nw=sorted_signal_nw((1:100));
%     fig=figure;     
%     histogram(truncated_sorted_signal_nw);
%     
%  fig=figure;     
%     histogram(signal_w(:));
%        
%     sorted_signal_w=sort(signal_w(:),'descend');    
%     truncated_sorted_signal_w=sorted_signal_w(1:100);
%     fig=figure;     
%     histogram(truncated_sorted_signal_w);
%     
% 
% % mean_wake=mean(truncated_sorted_signal_w)
% % mean_nowake=mean(truncated_sorted_signal_nw)
% % std_nowake=std(truncated_sorted_signal_w,1)
% % std_wake=std(truncated_sorted_signal_nw,1)
% % stn_nowake=(max_nowake_slices-mean_nowake)/std_nowake
% % stn_wake=(max_wake_slices-mean_nowake)/std_nowake
% % mean_stn=mean(stn_wake)
% % std_stn=std(stn_wake,1)
% % mean_stn-std_stn    
% 
% 
% mean_wake=mean(truncated_sorted_signal_w)
% mean_nowake=mean(truncated_sorted_signal_nw)
% std_nowake=std(truncated_sorted_signal_nw,1)
% std_wake=std(truncated_sorted_signal_w,1)
% % max_wake=max(truncated_sorted_signal_w,1)
% % max_nowake=max(truncated_sorted_signal_nw,1)
% % stn_nowake=(max_nowake'-mean_nowake)/std_nowake
% % stn_wake=(max_wake'-mean_nowake)/std_nowake
% stn_nowake=(truncated_sorted_signal_nw-mean_nowake)/std_nowake;
% stn_wake=(truncated_sorted_signal_w-mean_nowake)/std_nowake;
% mean_stn=mean(stn_wake)
% std_stn=std(stn_wake,1)
%     
% kurtosis(sorted_signal_nw)
% kurtosis(sorted_signal_w)


    

% cd('../wake_detection/slices_curv_2d/');

end


function [thetas, phis] = s2let_hpx_sampling_ring(nside)

% s2let_hpx_sampling_ring 
% Compute the Healpix sampling nodes
% Default usage :
%
%   [thetas, phis] = s2let_hpx_sampling_ring(nside)
%
% nside is the input Healpix resolution,
% thetas contains the colatitudes of the nodes,
% phis contains the longitudes of the nodes.
%
% S2LET package to perform Wavelets transform on the Sphere.
% Copyright (C) 2012-2015  Boris Leistedt & Jason McEwen
% See LICENSE.txt for license details

            %     npix = 12 * nside^2;
            %     thetas = zeros(npix,1);
            %     phis = zeros(npix,1);
            %     for i = 1:npix
            %        [t,p] = healpix_pix2ang_ring(i, nside);
            %        thetas(i) = t;
            %        phis(i) = p;
            %     end

                npix = 12 * nside^2;
                ipix = 1:npix;
                nl2 = 2 * nside;
                ncap = 2  * nside * (nside - 1);  % points in each polar cap, =0 for nside =1

                thetas = zeros(npix,1);
                phis = zeros(npix,1);

                ind = (ipix <= ncap); % North polar cap
                    hip   = ipix(ind) .* 0.5;
                    fihip = fix(hip);
        iring = fix( sqrt( hip - sqrt(fihip) ) ) + 1 ;% counted from North pole
        iphi  = ipix(ind) - 2*iring.*(iring - 1) ;

            thetas(ind) = acos( 1.0 - iring.^2.0 ./ (3.0*nside^2.0) );
        phis(ind)   = (iphi - 0.5) .* pi ./ (2.0.*iring);

    ind = (ipix <= nl2*(5*nside+1) & ipix > ncap); % Equatorial region

       ip    = ipix(ind) - ncap - 1;
       nl4   = 4*nside;
       iring = floor( ip / nl4 ) + nside; % counted from North pole
       iphi  = mod(ip,nl4) + 1;

       fodd  = 0.5 * (1 + mod(iring+nside,2)) ; % 1 if iring+nside is odd, 1/2 otherwise
       thetas(ind) = acos( (nl2 - iring) ./ (1.5*nside) );
       phis(ind)   = (iphi - fodd) * pi /(2.0*nside);

    ind = (ipix > nl2*(5*nside+1) & ipix > ncap); % South polar cap

       ip    = npix - ipix(ind) + 1;
       hip   = ip*0.5;
       fihip = round(hip);
       iring = floor( sqrt( hip - sqrt(fihip) ) ) + 1 ; % counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring.*(iring-1));

       thetas(ind) = acos( -1.0 + iring.^2 ./ (3.0*nside^2) );
       phis(ind)   = (iphi - 0.5) * pi ./ (2.0*iring);

end


function [x, y] = ssht_mollweide(thetas, phis,thetas_shift, phis_shift)
% ssht_mollweide - Compute Mollweide projection
%
% Compute Mollweide projection of spherical coordinates.
%
% Usage is given by
%
%   [x,y] = ssht_mollweide(thetas, phis)
%
% where thetas and phis are spherical coordinates and x and y are the
% projected Mollweide coordinates.
%
% Author: Jason McEwen (www.jasonmcewen.org)

MAX_ITERATIONS = 1e5;
TOL = 1e-10;


% Convert theta to longitude.
thetas = pi/2 - thetas;
phis = phis - pi;


% recenter
[rx(1,:),rx(2,:),rx(3,:)] = sph2cart(phis,thetas,1);

theta=pi/2+thetas_shift;
phi=-phis_shift;

% theta=pi/2;
% phi=pi/4;

Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

rx=Rz*rx;
rx=Ry*rx;

[phis,thetas,~] = cart2sph(rx(1,:),rx(2,:),rx(3,:));


t = thetas;
for it = 1:MAX_ITERATIONS

   dt = (t + sin(t) - pi.*sin(thetas)) ./ (1 + cos(t));
   t = t - dt;
   
   if(max(abs(dt)) < TOL)
      break;
   end
   
end
t = t/2;
x = 2 .* sqrt(2) ./ pi .* phis .* cos(t);
y = sqrt(2) .* sin(t);
end


function [f_index_max] = find_peaks(thetas, phis,f,NSIDE,number_of_maxima)

% number_of_maxima=4;

% peak=max(f);
% ave=mean(f);
% f_symmetric=f-ave;
% signal=max(f_symmetric);

f_size=length(f);
f(f_size/2:end)=[];

[f_sort,f_index_max]=sort(f,'descend');

% f_sort=flip(unique(sort(f,'descend')));


% f_symmetric_sort=fliplr(f_symmetric_sort);

% std_f=std(f);
% stn_f=signal/std_f;
% f_index_max=find(f_symmetric==signal);
% f_index_max=zeros(1,number_of_maxima);
% for pk=1:length(f_index_max)
%     f_index_max(1,pk)=find(f==f_sort(pk),1);
% end

phis_max=phis(f_index_max);
thetas_max=thetas(f_index_max);
% display(phis_max(1:10));
% display(thetas_max(1:10));

angle_cutoff=2*sqrt(4*pi/(12*(NSIDE/2)*(NSIDE/2)));

for f_verify =1:number_of_maxima
    
    f_index_max_aux=length(f_index_max);
    
%     display(f_index_max_aux)
        	v1=zeros(1,3);
            [v1(1,1),v1(1,2),v1(1,3)] = sph2cart(phis_max(f_verify),+pi/2-thetas_max(f_verify),1);
            v2=-v1;
    
%     parfor f_compare=f_verify+1:f_index_max_aux
    for f_compare=f_verify+1:f_index_max_aux

            w1=zeros(1,3);
            [w1(1,1),w1(1,2),w1(1,3)] = sph2cart(phis_max(f_compare),+pi/2-thetas_max(f_compare),1);
            w2=-w1;
            
            a=zeros(1,4);
            a(1,1)=acos(dot(v1, w1));
            a(1,2)=acos(dot(v1, w2));
            a(1,3)=acos(dot(v2, w1));
            a(1,4)=acos(dot(v2, w2));
            am=min(a);
            
            if (am<angle_cutoff)
                
                f_index_max(f_compare)=0;
                
            end
    end
    
    f_index_max(f_index_max==0) = [];
    
end

end
