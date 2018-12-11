function [ map,anali ] = slices_curv_2a3d( root,root_data_2d_in,root_data_2d_anali,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,numb_rand,slice,NSIDE,lev,sigma )

% (example) [ map ] = slices_curv_2d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,8,4,8,2,5 );

%(example) [ map ] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/ht/','/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024/','/home/asus/Dropbox/extras/storage/graham/ht/analy_cps32_1024/','4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m','/sample3001/','','3.000xv0.dat',1,1,1,32,8,2,5 );
%(example) for i=1:24; [ count_sum] = proj2d_cic_dm_data_out_part_rand_hpx_slices('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,1,1,i,2,2); end;
%(example)  for i=1:24; [ map ,anali] = slices_curv_2a3d('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,i,2,2,2,5 ); end;
 
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

F2 = ones(nb,nb,slice/2);
X2 = fftshift(ifft2(F2)) * sqrt(prod(size(F2)));
C2 = fdct3d_forward(X2);
E2 = cell(size(C2));
for s=1:length(C2)
    E{s} = cell(size(C2{s}));
    for w=1:length(C2{s})
        A2 = C2{s}{w};
        E2{s}{w} = sqrt(sum(sum(A2.*conj(A2))) / prod(size(A2)));
    end
end



F2=zeros(nb,nb,slices/2);
C_zero2 = fdct3d_forward(F2);


anali=zeros(slice,4,5);
anali_2a3d=zeros(slice,4,5);


mkdir(root_data_2d_anali);
mkdir(root_data_2d_anali,strcat(spec,aux_path));


    path1=strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/');
    
    path2=dir(char(strcat(path1,"*pv*")));
    path2={path2.name};
    path2=path2(1);
    
    path3='/2dproj/dm/';

    path_out_2d_filt_slices=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data_2d_filt_slices/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));

    
% for slice_id=1:1
for slice_id=1:slice

%     display(strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/'))
    
%     path1=strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/');
%     
%     path2=dir(char(strcat(path1,"*pv*")));
%     path2={path2.name};
%     path2=path2(1);
%     
%     path3='/2dproj/dm/';
    
    
    
    filename_2d=string(strcat(path1,path2,path3,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'))
    fid = fopen(filename_2d);
    map = fread(fid,[nb nb], 'float32','l') ;
    fclose(fid);
    
%     map2=map;
    map(map<=1)=1;%to remove problem with holes
    map2=log(map);
    
    C = fdct_wrapping(map2,0);
    Ct = C;
    for s = 1:length(C)
        for w = 1:length(C{s})
            Ct{s}{w} = C_zero{s}{w};
        end
    end
    
    for s = length(C)-lev:length(C)-1
        %                 thresh=0;
        thresh = sigma + sigma*(s == length(C));
        for w = 1:length(C{s})
            %                     Ct{s}{w} = C{s}{w};
            %                     Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
            Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
        end
    end
    
    map_filt= real(ifdct_wrapping(Ct,0));
    
    map_filt(1:16,:)=0;
    map_filt(end-16:end,:)=0;
    map_filt(:,1:16)=0;
    map_filt(:,end-16:end)=0;
    
            
            map_3d_slices_filt2d(:,:,slice_id) =real(ifdct_wrapping(Ct,0));    
    
    fileID = fopen(strcat(path_out_2d_filt_slices,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
    fwrite(fileID,map_filt, 'float32','l');
    fclose(fileID);
    
    anali(slice_id,1,:)=[max(map_filt(:)),std(map_filt(:)),max(map_filt(:))/std(map_filt(:)),kurtosis(kurtosis(map_filt)),kurtosis(map_filt(:))];
    
    theta = 0:180;
    [R,xp] = radon(map_filt,theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
        
	boudary_removal_factor=1024/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R(1,:))
        [dc_dwt,levels] = wavedec(R(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev);
                        D(floor((1237-256)/boudary_removal_factor):end)=0;
                        D(1:floor((217+256)/boudary_removal_factor))=0;
%         D(1237-256:end)=0;
%         D(1:217+256)=0;
        
        R_nor_filt(:,i)=D;
    end
                R_nor_filt(floor((1237-256)/boudary_removal_factor):end,:)=[];
                R_nor_filt(1:floor((217+256)/boudary_removal_factor),:)=[];
%     R_nor_filt(1237-256:end,:)=[];
%     R_nor_filt(1:217+256,:)=[];
    
    
    
    anali(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
    
end


        for partition=0:1
            
            C_aux=map_3d_slices_filt2d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2);
            
            C = fdct3d_forward(C_aux);
            
            Ct = C;
            for s = 1:length(C)
                for w = 1:length(C{s})
                    Ct{s}{w} = C_zero2{s}{w};
                end
            end
            
            aux_count=1;
            for s = 1:1
                %                 thresh=0;
                thresh = sigma + sigma*(s == length(C));
                for w = 1:length(C{s})
                    Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > 0*E2{s}{w})
                    %                                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w}.* ((C{s}{w}) > -thresh*E{s}{w}&(C{s}{w}) < thresh*E{s}{w});
                    %                   Ct{s}{w} = C{s}{w};
                    %                   curv(w_nw,sample,slice_id,w,aux_count)=kurtosis(abs(C{s}{w}(:)));
                end
                %                 curv2(w_nw,sample,slice_id,aux_count)=kurtosis(curv(w_nw,sample,slice_id,:,aux_count));
                aux_count=aux_count+1;
            end
            
            map_3d_slices_filt2a3d(:,:,(partition)*slices/2+1:(partition)*slices/2+slices/2) = real(fdct3d_inverse(Ct));
            
        end
        
        
for slice_id=1:slice

    map_filt=map_3d_slices_filt2a3d(:,:,slice_id);
        
%     map_filt(1:16,:)=0;
%     map_filt(end-16:end,:)=0;
%     map_filt(:,1:16)=0;
%     map_filt(:,end-16:end)=0;
                
    
    fileID = fopen(strcat(path_out_2d_filt_slices,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_2a3dcurvfilt_z',num2str(z_glob),'_data_sl',num2str(slice_id),'.bin'),'w');
    fwrite(fileID,map_filt, 'float32','l');
    fclose(fileID);
    
    anali_2a3d(slice_id,1,:)=[max(map_filt(:)),std(map_filt(:)),max(map_filt(:))/std(map_filt(:)),kurtosis(kurtosis(map_filt)),kurtosis(map_filt(:))];
    
    theta = 0:180;
    [R,xp] = radon(map_filt,theta);
    
    unit=ones(nb);
    [R_u,xp] = radon(unit,theta);
    
    frac_cut=0.5;
    R_nor=R;
    R_nor(R_u>nb*frac_cut)=R_nor(R_u>nb*frac_cut)./R_u(R_u>nb*frac_cut);
    R_nor(R_u<=nb*frac_cut)=0;
        
	boudary_removal_factor=1024/nb;
    
    n_levels=floor(log2(length(R_nor(:,1))));
    R_nor_filt=zeros(size(R_nor));
    for i=1:length(R(1,:))
        [dc_dwt,levels] = wavedec(R(:,i),n_levels,'db1');
        D = wrcoef('d',dc_dwt,levels,'db1',lev);
                        D(floor((1237-256)/boudary_removal_factor):end)=0;
                        D(1:floor((217+256)/boudary_removal_factor))=0;
%         D(1237-256:end)=0;
%         D(1:217+256)=0;
        
        R_nor_filt(:,i)=D;
    end
                R_nor_filt(floor((1237-256)/boudary_removal_factor):end,:)=[];
                R_nor_filt(1:floor((217+256)/boudary_removal_factor),:)=[];
%     R_nor_filt(1237-256:end,:)=[];
%     R_nor_filt(1:217+256,:)=[];
    
    
    
    anali_2a3d(slice_id,2,:)=[max(R(:)),std(R(:)),max(R(:))/std(R(:)),kurtosis(kurtosis(R)),kurtosis(R(:))];
    anali_2a3d(slice_id,3,:)=[max(R_nor(:)),std(R_nor(:)),max(R_nor(:))/std(R_nor(:)),kurtosis(kurtosis(R_nor)),kurtosis(R_nor(:))];
    anali_2a3d(slice_id,4,:)=[max(R_nor_filt(:)),std(R_nor_filt(:)),max(R_nor_filt(:))/std(R_nor_filt(:)),kurtosis(kurtosis(R_nor_filt)),kurtosis(R_nor_filt(:))];
    
    
end
        


% mkdir(root_data_2d_anali);
% mkdir(root_data_2d_anali,strcat(spec,aux_path));

% strcat(root_data_2d_in,spec,aux_path,'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/')

path_out=string(strcat(strcat(root_data_2d_anali,spec,aux_path),'data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'));
string(strcat(root_data_2d_anali,spec,aux_path))
string(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/'))
mkdir(char(strcat(root_data_2d_anali,spec,aux_path)),char(strcat('data/',aux_path_out,num2str(lenght_factor),'lf_',num2str(resol_factor),'rf','/NSIDE_',num2str(NSIDE),'/anglid_',num2str(numb_rand),'/',path2,'/2dproj/dm/')));
% mkdir(char(strcat(root_data_2d_anali,spec,aux_path)));

% dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl',num2str(count_slice),'.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_curv_sl','.txt'),anali,'delimiter','\t');
dlmwrite(strcat(path_out,'_',num2str(find(str2num(char(redshift_list))==z_glob)),'_2dproj_z',num2str(z_glob),'_data_2a3dcurv_sl','.txt'),anali_2a3d,'delimiter','\t');


cd('../wake_detection/slices_curv_2d/');


end

