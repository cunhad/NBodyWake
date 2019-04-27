function [ map,anali ] = slices_curv_2a3d_locangreal( root,root_data_2d_in,root_data_2d_out,root_data_2d_anali_out,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,numb_rand,slices,lev,lev_rid,lev_3d,lev_3drig,step_of_degree,wavel_removal_factor,NSIDE,partition2d,partition3rd,sum_depth)

% (example) [ map ] = slices_curv_2d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,8,4,8,2,5 );

%(example) [ map ] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/ht/','/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/','/home/asus/Dropbox/extras/storage/graham/ht/analy_cps32_1024/','4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m','/sample3001/','','3.000xv0.dat',1,1,1,32,8,2,5 );
%(example) for i=1:24; [ count_sum] = proj2d_cic_dm_data_out_part_rand_hpx_slices('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,1,1,i,2,2); end;
%(example)  for i=1:24; [ map ,anali] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,i,2,2,2,5 ); end;
%  
% 
% root='/home/asus/Dropbox/extras/storage/graham/small_res/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/small_res/data_test3/';
% root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/small_res/data3/';
% root_data_2d_anali_out='/home/asus/Dropbox/extras/storage/graham/small_res/anali3/';
% spec='64Mpc_256c_128p_zi63_nowakem';
% aux_path='/sample2001/';
% aux_path_out='';
% filename='10.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% numb_rand=1;
% slices=32;
% lev=2;
% lev_rid=2;
% lev_3d=1;
% lev_3drig=2;
% step_of_degree=(64/180);
% wavel_removal_factor=1/2;
% NSIDE=8;
% partition2d=1;
% partition3rd=2;
% sum_depth=4;
% 
% 
% root='/home/asus/Dropbox/extras/storage/graham/ht/';
% root_data_2d_in='/home/asus/Dropbox/extras/storage/graham/ht/data_hpx_2d/';
% root_data_2d_out='/home/asus/Dropbox/extras/storage/graham/ht/data_tst/';
% root_data_2d_anali_out='/home/asus/Dropbox/extras/storage/graham/ht/data_tst_anali/';
% spec='4Mpc_2048c_1024p_zi63_nowakem';
% aux_path='/sample3004/';
% aux_path_out='';
% filename='3.000xv0.dat';
% lenght_factor=1;
% resol_factor=1;
% numb_rand=95;
% slices=32;
% lev=2;
% lev_rid=1;
% lev_3d=1;
% lev_3drig=1;
% step_of_degree=(180/256);
% wavel_removal_factor=1/2;
% NSIDE=8;
% partition2d=1;
% partition3rd=2;
% sum_depth=4;


%%%%Problems

%%%%    anali3_curv definition

% 

% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
% addpath(genpath('/home/asus/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));


 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_cpp/mex/'));
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct_wrapping_matlab'));
 addpath(genpath('/home/cunhad/projects/rrg-rhb/cunhad/Programs/CurveLab_matlab_3d-0.1-2.1.3/fdct3d'));


cd('../../preprocessing');

[~,redshift_list,nodes_list,size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw] = preprocessing_info(root,spec,aux_path );

% [ size_box,nc,np,zi,wake_or_no_wake,multiplicity_of_files,Gmu,ziw,z,path_file_in,~ ] = preprocessing_part(root,spec,aux_path,filename,1024,1);


z_glob=str2num(filename(1:end-7))


nb=np*resol_factor/lenght_factor;
nc_anal3d=nb/partition2d;


% filename='_2dproj_z3_data_sl';    

F = ones(nb);
X = fftshift(ifft2(F)) * sqrt(prod(size(F)));
%X = F * sqrt(prod(size(F)));
%C = fdct_wrapping(X,0);
C = fdct_wrapping(X,0,2);
%C = fdct_wrapping(F,0);
E = cell(size(C));
for s=1:length(C)
    E{s} = cell(size(C{s}));
    for w=1:length(C{s})
        A = C{s}{w};
        E{s}{w} = sqrt(sum(sum(A.*conj(A))) / prod(size(A)));
    end
end

F=zeros(nb);
C_zero = fdct_wrapping(F,0);

F2 = ones(nb,nb,slices/partition3rd);
% F2 = ones(new_nc,new_nc,slices);
X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
C2 = fdct3d_forward(X2);
E2 = cell(size(C2));
for s=1:length(C2)
    E2{s} = cell(size(C2{s}));
    for w=1:length(C2{s})
        A2 = C2{s}{w};
        E2{s}{w} = sqrt(sum(sum(sum(A2.*conj(A2)))) / prod(size(A2)));
    end
end
%
%
%
F2=zeros(nb,nb,slices/partition3rd);
C_zero2 = fdct3d_forward(F2);


anali=zeros(slices,4,5);
anali_curv=zeros(slices,5,lev);
% anali3_curv=zeros(5,(lev_3d)*partition2d*partition2d*partition3rd);
anali2a3=zeros(slices,4,5);
anali2a3_curv=zeros(slices,5,lev);

if ~ismember(1,sum_depth)
    anali_depth=zeros(slices/sum_depth,4,5);
    anali2a3_depth=zeros(slices/sum_depth,4,5);
    anali2a3_curv_depth=zeros(slices/sum_depth,5,lev);
end

mkdir(root_data_2d_anali_out);
mkdir(root_data_2d_anali_out,strcat(spec,aux_path));


path1=strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/');

path2=dir(char(strcat(path1,"*pv*")));
path2={path2.name};
path2=path2(1);

path3='/2dproj/dm/';


path_data_2d_anali_out=string(strcat(strcat(root_data_2d_anali_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
mkdir(char(strcat(root_data_2d_anali_out,spec,aux_path)),char(strcat('data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));

if ~isempty(root_data_2d_out)
    
    mkdir(root_data_2d_out);
    mkdir(root_data_2d_out,strcat(spec,aux_path));
        
    path_data_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    mkdir(char(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
    
    path_data3_out=string(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt3_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
    mkdir(char(strcat(strcat(root_data_2d_out,spec,aux_path),'data_2d_filt3_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
end


map_3d_slices_filt2d=zeros(nb,nb,slices);
map_3d_slices_filt3d=zeros(nc_anal3d,nc_anal3d,slices);

% for slice_id=1:1
for slice_id=1:slices

    
    filename_2d=string(strcat(path1,path2,path3,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'))
    fid = fopen(filename_2d);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    
%     map2=map;
    map(map<=1)=1;%to remove problem with holes
    map_3d_slices=log(map);
    
    C = fdct_wrapping(map_3d_slices,0);
    Ct = C;
    Ct2=C;    
    for s = 1:length(C)
        for w = 1:length(C{s})
            Ct{s}{w} = C_zero{s}{w};
            Ct2{s}{w} = C_zero{s}{w};
        end
    end
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        
        %right part
                sz=cell(length(C{s}),1);
                for w = 1:length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1}(1)
                    for j=1:sz{1}(2)
                        for w = 1:length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1}(2));
                            a(w)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a);
%                         treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1:length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w)*(C{s}{w}(i_c,j_c))*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1:length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                
                %right2 part
                for w = 1+length(C{s})/4:2*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+length(C{s})/4}(1)
                    for j=1:sz{1+length(C{s})/4}(2)
                        for w = 1+length(C{s})/4:2*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+length(C{s})/4}(2));
                            a(w-length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+length(C{s})/4:2*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-length(C{s})/4),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1+length(C{s})/4:2*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                %right3 part
                for w = 1+2*length(C{s})/4:3*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+2*length(C{s})/4}(1)
                    for j=1:sz{1+2*length(C{s})/4}(2)
                        for w = 1+2*length(C{s})/4:3*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+2*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+2*length(C{s})/4}(2));
                            a(w-2*length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+2*length(C{s})/4:3*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+2*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+2*length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-2*length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-2*length(C{s})/4),Ct{s}{w}(i_c,j_c));
                        end
                    end
                end
                
                for w = 1+2*length(C{s})/4:3*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                %right4 part
                for w = 1+3*length(C{s})/4:4*length(C{s})/4
                    sz{w}=size(C{s}{w});
                end
                for  i=1:sz{1+3*length(C{s})/4}(1)
                    for j=1:sz{1+3*length(C{s})/4}(2)
                        for w = 1+3*length(C{s})/4:4*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+3*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+3*length(C{s})/4}(2));
                            a(w-3*length(C{s})/4)=abs(C{s}{w}(i_c,j_c))/E{s}{w};
                        end
%                         std_a=std(a);
                        avr_a=mean(a);
                        %             treash_logic=double(a>=avr_a+sigma*std_a)   ;
                        %             treash_logic=double(a==max(a));
                        %             normalization=(a-mean(a)/std_a);
                        normalization=(a-avr_a)/avr_a;
                        for w = 1+3*length(C{s})/4:4*length(C{s})/4
                            i_c=ceil(i*sz{w}(1)/sz{1+3*length(C{s})/4}(1));
                            j_c=ceil(j*sz{w}(2)/sz{1+3*length(C{s})/4}(2));
                            %                 Ct{s}{w}(i_c,j_c)=treash_logic(w-3*length(C{s})/4)*C{s}{w}(i_c,j_c)*double(Ct{s}{w}(i_c,j_c)~=0);
                            Ct{s}{w}(i_c,j_c)=max(normalization(w-3*length(C{s})/4),Ct{s}{w}(i_c,j_c));
                            
                        end
                    end
                end
                
                for w = 1+3*length(C{s})/4:4*length(C{s})/4
                    Ct2{s}{w} = Ct{s}{w}.*C{s}{w};
                end
                
                
                a=[];
                treash_logic=[];
                normalization=[];
                
        
        
        for w = 1:length(C{s})
            %                     Ct{s}{w} = C{s}{w};
            %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
            %                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E{s}{w});
%             Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
            %                 Ct{s}{w} = C{s}{w};
            %              curv(slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
            
%             ave=mean(abs(Ct{s}{w}(:)));
%             sig=std(abs(Ct{s}{w}(:)))^2;
%             del=skewness(abs(Ct{s}{w}(:)))*(sig^(3/2));
%             rho=kurtosis(abs(Ct{s}{w}(:)))*(sig^(2));

            ave=mean(abs(Ct2{s}{w}(:))/E{s}{w});
            sig=std(abs(Ct2{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(Ct2{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(Ct2{s}{w}(:))/E{s}{w})*(sig^(2));

            A_=ave;
            B_=sig+A_^2;
            C_=del+3*B_*A_-2*A_^3;
            D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
            
            size_n=prod(size(Ct2{s}{w}(:)));
            size_nt=size_nt+size_n;
            
            As=A_*size_n+As;
            Bs=B_*size_n+Bs;
            Cs=C_*size_n+Cs;
            Ds=D_*size_n+Ds;
            
            
            
        end
        %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
        AT=As/size_nt;
        BT=Bs/size_nt;
        CT=Cs/size_nt;
        DT=Ds/size_nt;
        
        sigma_t=BT-AT^2;
        delta_t=CT-3*BT*AT+2*AT^3;
        rho_t=DT-4*CT*AT+6*BT*AT^2-3*AT^4;
        
        var_t=sigma_t;
        skew_t=delta_t/var_t^(3/2);
        kurt_t=rho_t/var_t^2;
        
        curv_1(aux_count)=AT;
        curv_2(aux_count)=var_t;
        curv_3(aux_count)=skew_t;
        curv_4(aux_count)=kurt_t;
        curv_5(aux_count)=rho_t;
        
        aux_count=aux_count+1;
    end
    
    
    
    
    
    
    
    
    map_3d_slices_filt2d(:,:,slice_id) = real(ifdct_wrapping(Ct2,0));    
    
    map_3d_slices_filt2d(1:16,:,slice_id)=0;
    map_3d_slices_filt2d(end-16:end,:,slice_id)=0;
    map_3d_slices_filt2d(:,1:16,slice_id)=0;
    map_3d_slices_filt2d(:,end-16:end,slice_id)=0;
    
    
    if ~isempty(root_data_2d_out)
    
    fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
    fwrite(fileID,map_3d_slices_filt2d(:,:,slice_id), 'float32','l');
    fclose(fileID);
    end
    
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filt2d(:,:,slice_id),theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
    
    %     boudary_removal_factor=2048/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R_nor(1,:))
        %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev_rid);
        %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
        %                 D(1:floor((448+200)/boudary_removal_factor))=0;
        D(length(xp)-nb*wavel_removal_factor:end)=0;
        D(1:nb*wavel_removal_factor)=0;
        
        R_nor_filt(:,i)=D;
    end
    %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
    %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
    R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
    R_nor_filt(1:nb*wavel_removal_factor,:)=[];
    
    this=map_3d_slices_filt2d(:,:,slice_id);
    anali(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    anali(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
    anali_curv(slice_id,1,:)=curv_1;
    anali_curv(slice_id,2,:)=curv_2;
    anali_curv(slice_id,3,:)=curv_3;
    anali_curv(slice_id,4,:)=curv_4;
    anali_curv(slice_id,5,:)=curv_5;
    
end

dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_curv_z',num2str(z_glob),'_anali_curv.txt'),anali_curv,'delimiter','\t');



% now we do the depth analysis

if ~ismember(1,sum_depth)
    
    slices_depth=slices/sum_depth;
%     map_3d_slices_depth=zeros(nb,nb,slices_depth);
    map_3d_slices_filt2d_depth=zeros(nb,nb,slices_depth);

%         map_3d_slices_depth=map_3d_slices_depth+map_3d_slices;
%         map_3d_slices_filt3d_depth=map_2d_slices_filt2d_depth+map_3d_slices_filt2d;
%         slice_depth_id=floor(slice_id/sum_depth);

    for slice_depth_id=1:slices_depth
        
%         map_3d_slices_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt2d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3);
        map_3d_slices_filt2d_depth(:,:,slice_depth_id)=sum(map_3d_slices_filt2d(:,:,(slice_depth_id-1)*sum_depth+1:(slice_depth_id)*sum_depth),3); 
        
        
        
        theta = 0:step_of_degree:180;
        [R,xp] = radon(map_3d_slices_filt2d_depth(:,:,slice_depth_id),theta);
        
        unit=ones(nb);
        [R_u,xp] = radon(unit,theta);
        
        frac_cut=0.5;
        R_nor=R;
        R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
        R_nor(R_u<=nb*frac_cut)=0;
        
        %     boudary_removal_factor=2048/nb;
        
        n_levels=floor(log2(length(R_nor(:,1))));
        R_nor_filt=zeros(size(R_nor));
        for i=1:length(R_nor(1,:))
            %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
            [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
            D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
            %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
            %                 D(1:floor((448+200)/boudary_removal_factor))=0;
            D(length(xp)-nb*wavel_removal_factor:end)=0;
            D(1:nb*wavel_removal_factor)=0;
            
            R_nor_filt(:,i)=D;
        end
        %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
        %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
        R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
        R_nor_filt(1:nb*wavel_removal_factor,:)=[];
                      
        
         
        this=map_3d_slices_filt2d_depth(:,:,slice_depth_id);
        
        anali_depth(slice_depth_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];      
        anali_depth(slice_depth_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
        anali_depth(slice_depth_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
        anali_depth(slice_depth_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];                            
        
        
    end

    dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurv_z',num2str(z_glob),'_anali_depth',num2str(sum_depth),'.txt'),anali_depth,'delimiter','\t');
    
end



for partition=0:partition3rd-1
    for partition_x=0:partition2d-1
        for partition_y=0:partition2d-1
            
            C_aux=map_3d_slices_filt2d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd);
            
            
            C = fdct3d_forward(C_aux);
            
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero2{s}{w};
                end
            end
            
            aux_count=1;
            for s =length(Ct)-lev_3d:length(Ct)-1
                %                 thresh=0;
%                 thresh = sigma + sigma*(s == length(C));
                
                
                As=0;
                Bs=0;
                Cs=0;
                Ds=0;
                
                size_nt=0;
                
                for w = 1:length(C{s})
                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E2{s}{w})
                    %                                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w};
                    %                   curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                    
                    
                    ave=mean(abs(Ct{s}{w}(:))/E2{s}{w});
                    sig=std(abs(Ct{s}{w}(:))/E2{s}{w})^2;
                    del=skewness(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(3/2));
                    rho=kurtosis(abs(Ct{s}{w}(:))/E2{s}{w})*(sig^(2));
                    
                    A_=ave;
                    B_=sig+A_^2;
                    C_=del+3*B_*A_-2*A_^3;
                    D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
                    
                    size_n=prod(size(Ct{s}{w}(:)));
                    size_nt=size_nt+size_n;
                    
                    As=A_*size_n+As;
                    Bs=B_*size_n+Bs;
                    Cs=C_*size_n+Cs;
                    Ds=D_*size_n+Ds;
                    
                    
                end
                %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                
                AT=As/size_nt;
                BT=Bs/size_nt;
                CT=Cs/size_nt;
                DT=Ds/size_nt;
                
                
                sigma_t=BT-AT^2;
                delta_t=CT-3*BT*AT+2*AT^3;
                rho_t=DT-4*CT*AT+6*BT*AT^2-3*AT^4;
                
                var_t=sigma_t;
                skew_t=delta_t/var_t^(3/2);
                kurt_t=rho_t/var_t^2;
                
                curv_1(aux_count)=AT;
                curv_2(aux_count)=var_t;
                curv_3(aux_count)=skew_t;
                curv_4(aux_count)=kurt_t;
                curv_5(aux_count)=rho_t;
                
                aux_count=aux_count+1;
            end
            
            %             map_3d_slices_filt2a3d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2) = real(fdct3d_inverse(Ct));
            map_3d_slices_filt3d((partition_x)*nb/partition2d+1:(partition_x)*nb/partition2d+nb/partition2d,(partition_y)*nb/partition2d+1:(partition_y)*nb/partition2d+nb/partition2d,(partition)*slices/partition3rd+1:(partition)*slices/partition3rd+slices/partition3rd) = real(fdct3d_inverse(Ct));
            
            anali3_curv_aux(1,partition+1,partition_x+1,partition_y+1,:)=curv_1(:);
            anali3_curv_aux(2,partition+1,partition_x+1,partition_y+1,:)=curv_2(:);
            anali3_curv_aux(3,partition+1,partition_x+1,partition_y+1,:)=curv_3(:);
            anali3_curv_aux(4,partition+1,partition_x+1,partition_y+1,:)=curv_4(:);
            anali3_curv_aux(5,partition+1,partition_x+1,partition_y+1,:)=curv_5(:);
            
        end
    end
end

anali3_curv1_=anali3_curv_aux(1,:,:,:,:);
anali3_curv2_=anali3_curv_aux(2,:,:,:,:);
anali3_curv3_=anali3_curv_aux(3,:,:,:,:);
anali3_curv4_=anali3_curv_aux(4,:,:,:,:);
anali3_curv5_=anali3_curv_aux(5,:,:,:,:);


anali3_curv1=anali3_curv1_(:);
anali3_curv2=anali3_curv2_(:);
anali3_curv3=anali3_curv3_(:);
anali3_curv4=anali3_curv4_(:);
anali3_curv5=anali3_curv5_(:);



anali3_curv(1,:)=anali3_curv1;
anali3_curv(2,:)=anali3_curv2;
anali3_curv(3,:)=anali3_curv3;
anali3_curv(4,:)=anali3_curv4;
anali3_curv(5,:)=anali3_curv5;


dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali3_curv.txt'),anali3_curv,'delimiter','\t');
        

%now the analysis of the 3d filtered slices

for slice_id=1:slices            
    
    C = fdct_wrapping(map_3d_slices_filt3d(:,:,slice_id),0);
    Ct = C;
    Ct2=C;    
   
    
    aux_count=1;
    for s = length(C)-lev:length(C)-1
        
        As=0;
        Bs=0;
        Cs=0;
        Ds=0;
        
        size_nt=0;
        
        
        for w = 1:length(C{s})

            ave=mean(abs(C{s}{w}(:))/E{s}{w});
            sig=std(abs(C{s}{w}(:))/E{s}{w})^2;
            del=skewness(abs(C{s}{w}(:))/E{s}{w})*(sig^(3/2));
            rho=kurtosis(abs(C{s}{w}(:))/E{s}{w})*(sig^(2));

            A_=ave;
            B_=sig+A_^2;
            C_=del+3*B_*A_-2*A_^3;
            D_=rho+4*C_*A_-6*B_*A_^2+3*A_^4;
            
            size_n=prod(size(C{s}{w}(:)));
            size_nt=size_nt+size_n;
            
            As=A_*size_n+As;
            Bs=B_*size_n+Bs;
            Cs=C_*size_n+Cs;
            Ds=D_*size_n+Ds;
            
            
            
        end
        %          curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(slice_id,:,aux_count));
        AT=As/size_nt;
        BT=Bs/size_nt;
        CT=Cs/size_nt;
        DT=Ds/size_nt;
        
        sigma_t=BT-AT^2;
        delta_t=CT-3*BT*AT+2*AT^3;
        rho_t=DT-4*CT*AT+6*BT*AT^2-3*AT^4;
        
        var_t=sigma_t;
        skew_t=delta_t/var_t^(3/2);
        kurt_t=rho_t/var_t^2;
        
        curv3_1(aux_count)=AT;
        curv3_2(aux_count)=var_t;
        curv3_3(aux_count)=skew_t;
        curv3_4(aux_count)=kurt_t;
        curv3_5(aux_count)=rho_t;
        
        aux_count=aux_count+1;
    end
    
    theta = 0:step_of_degree:180;
    [R,xp] = radon(map_3d_slices_filt3d(:,:,slice_id),theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
    
    %     boudary_removal_factor=2048/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R_nor(1,:))
        %                 [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        [dc_dwt,levels] = wavedec(R_nor(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev_3drig);
        %                 D(floor((2465-200)/boudary_removal_factor):end)=0;
        %                 D(1:floor((448+200)/boudary_removal_factor))=0;
        D(length(xp)-nb*wavel_removal_factor:end)=0;
        D(1:nb*wavel_removal_factor)=0;
        
        R_nor_filt(:,i)=D;
    end
    %             R_nor_filt(floor((2465-200)/boudary_removal_factor):end,:)=[];
    %             R_nor_filt(1:floor((448+200)/boudary_removal_factor),:)=[];
    R_nor_filt(length(xp)-nb*wavel_removal_factor:end,:)=[];
    R_nor_filt(1:nb*wavel_removal_factor,:)=[];
    
       
    
        this=map_3d_slices_filt3d(:,:,slice_id);
    
    anali2a3(slice_id,1,:)=[max(this(:)),std(this(:)),max(this(:))/std(this(:)),kurtosis(kurtosis(this)),kurtosis(this(:))];
    anali2a3(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali2a3(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali2a3(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];

    anali2a3_curv(slice_id,1,:)=curv3_1;
    anali2a3_curv(slice_id,2,:)=curv3_2;
    anali2a3_curv(slice_id,3,:)=curv3_3;
    anali2a3_curv(slice_id,4,:)=curv3_4;
    anali2a3_curv(slice_id,5,:)=curv3_5;
    
%     if ~isempty(root_data_2d_out)
%         
%         fileID = fopen(strcat(path_data_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
%         fwrite(fileID,map_3d_slices_filt3d(:,:,slice_id), 'float32','l');
%         fclose(fileID);
%         
%     end
    
    
    if ~isempty(root_data_2d_out)
    
    fileID = fopen(strcat(path_data3_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3a2dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
    fwrite(fileID,map_3d_slices_filt3d(:,:,slice_id), 'float32','l');
    fclose(fileID);
    end


end

dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali2a3.txt'),anali2a3,'delimiter','\t');
dlmwrite(strcat(path_data_2d_anali_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_3dcurv_z',num2str(z_glob),'_anali2a3_curv.txt'),anali2a3_curv,'delimiter','\t');

% 
% 
% path_out=string(strcat(strcat(root_data_2d_anali_out,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
% string(strcat(root_data_2d_anali_out,spec,aux_path))
% string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
% mkdir(char(strcat(root_data_2d_anali_out,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
% % mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));
% 
% % dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
% dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl','.txt'),anali,'delimiter','\t');
% dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_2a3dcurv_sl','.txt'),anali_2a3d,'delimiter','\t');


cd('../wake_detection/slices_curv_2d/');


end

