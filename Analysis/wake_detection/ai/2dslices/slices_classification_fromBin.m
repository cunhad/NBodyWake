function [] = slices_classification_fromBin()



clearvars

filename='_2dproj_z3_data_sl32All';
nc=512;
NSIDE=4;
depth=4
slices=32;

specs_path_list_nowake='~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');
% sample_list_wake=sort_nat(sample_list_wake)

sample_id_range=[1 : length(sample_list_nowake)];

count=1;

% for w_nw=1:2
for w_nw=2
    
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
%         ch='_7';
        ch='_1';
        coul='b';
    else
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
%         ch='_4';
        ch='_1';
        coul='r';
    end
    
    
    for sample = 1:length(sample_id_range)
        %     for sample = 1
        
           list_of_angle_paths=dir(char(strcat(specs_path_list,'/',string(sample_list(sample)),'/data/1lf_0.5rf/NSIDE_',num2str(NSIDE),'/anglid_*')));
        list_of_angle_paths={list_of_angle_paths.name};
        angle_id_list=1:length(list_of_angle_paths);
        
        for angle_id=angle_id_list

            %             angle_id=1
            
            for slice_id=1:slices
%             for slice_id=1:1
%             for slice_id=1      
                
                
                path1=strcat(specs_path_list,'/',string(sample_list(sample)),'/data/1lf_0.5rf/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/');
                
                path2=dir(char(strcat(path1,"*pv*")));
                path2={path2.name};
                path2=path2(1);
                
                path3='/2dproj/dm/';

                file=strcat(specs_path_list,'/',string(sample_list(sample)),'/data/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/',string(path2),'/2dproj/dm/',ch,filename,num2str(slice_id),'.bin');                
%                 file=strcat(specs_path_list,'/',string(sample_list(sample)),'/data/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/',string(path2),'/2dproj/dm/',ch,filename,num2str(slice_id),'.bin');
                
                
                file=char(file);
                  
                list{count}=file;
                count=count+1;
                         
            end     
            
        end       
        
    end            
    
end


label= categorical(abs(double(contains(string(list),'nowake'))-1));

% imds = imageDatastore(list,'ReadFcn',@read_slices_binAll,'FileExtensions','.bin','Labels',label);
% 
% 
% imds = imageDatastore({list{1}},'ReadFcn',@read_slices_binAll,'FileExtensions','.bin','Labels',label(1));
% 
% 
% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% imds = imageDatastore({filename},'ReadFcn',@read_slices_binAll,'FileExtensions','.bin','Labels',label(1));
% 
% imds = imageDatastore({filename},'ReadFcn',@(filename) read_slices_binAll(filename,1),'FileExtensions','.bin','Labels',label(1));
% 
% imds1 = imageDatastore({filename},'ReadFcn',@(filename) read_slices_binAll(filename,1),'FileExtensions','.bin','Labels',label(1));
% imds2 = imageDatastore({filename},'ReadFcn',@(filename) read_slices_binAll(filename,2),'FileExtensions','.bin','Labels',label(1));
% 
% 
% imds = combine(imds1,imds2)
% 
% imds = imageDatastore({filename},'ReadFcn',@(filename) read_slices_binAll(filename,1),'FileExtensions','.bin','Labels',label(1));
% 
% 
% dsseq = combine(imds1,imds2,ReadOrder="sequential")
% 
% a = readimage(dsseq,1);
% 
% 
% 
% 
% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% 
% imds1 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl1,'FileExtensions','.bin','Labels',label(1));
% imds2 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl2,'FileExtensions','.bin','Labels',label(1));
% 
% imds1 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl1,'FileExtensions','.bin','LabelSource' ,'foldernames');
% imds2 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl2,'FileExtensions','.bin','LabelSource' ,'foldernames');
% 
% imds = combine(imds1,imds2)
% 
% 
% imds = transform(imds1,imds2, @(x1,x2) [x1,x2])
% 
% tile = read(imds);
% 
% 
% transform(imds1,imds2, @(x) imresize(x,32))
% 
% @(x) imresize(x,targetSize)
% 
% 
% 
% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% 
% imds1 = imageDatastore({filename},'ReadFcn',@read_slices_bin_slices,'FileExtensions','.bin','Labels',label(1:2),'ReadSize',2);
% 
% 
% list_ = list;
% list_{end+1} = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% label= categorical(abs(double(contains(string(list),'nowake'))-1));
% 
% imds = imageDatastore(list_,'ReadFcn',@read_slices_binAll,'FileExtensions','.bin','Labels',label);
% 
% 
% matlab.io.datastore.ImageDatastore
% 
% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% imds = MySequenceDatastore(filename);
% 
% 
% filename1 = matlab.io.datastore.DsFileReader(filename);
% imds1 = imageDatastore({filename1},'ReadFcn',@read_slices_bin_sl1,'FileExtensions','.bin','LabelSource' ,'foldernames');
% 
% 
% 
% 
% filename = '~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin'
% 
% imds1 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl1,'FileExtensions','.bin','Labels',label(1));
% imds2 = imageDatastore({filename},'ReadFcn',@read_slices_bin_sl2,'FileExtensions','.bin','Labels',label(1));
% imds = combine(imds1,imds2);





imds = imageDatastore(list,'ReadFcn',@read_slices_bin_slices,'FileExtensions','.bin','Labels',label);






% imdsnew = transform(imds1,imds2, @(x) x)

% imds = MySequenceDatastore({filename},'ReadFcn',@read_slices_binAll,'FileExtensions','.bin','Labels',label(1));

% imdsnew = transform(imds1,imds2, @(x) x)


% imds = MySequenceDatastore({filename});


% imds = MySequenceDatastore({list{1});


% bimds = blockedImageDatastore({filename})

% a = readimage(imds,1);
% figure; imshow(a);


figure; image(readimage(imds,1)); colorbar;

% 
% dataOut = read(imds);
% 
% imshow(imtile(dataOut));


% categorical([1,1]);
% imds.UnderlyingDatastores = 

% Define the convolutional neural network architecture


layers = [
% imageInputLayer([1024 1024 3])
imageInputLayer([512 512 1])

convolution2dLayer(2,1,'Padding','same')
batchNormalizationLayer
reluLayer
convolution2dLayer(2,2,'Padding','same')
batchNormalizationLayer
reluLayer
maxPooling2dLayer(2,'Stride',2)
convolution2dLayer(4,4,'Padding','same')
batchNormalizationLayer
reluLayer
fullyConnectedLayer(1)
softmaxLayer
classificationLayer];



%minibatch size

miniBatchSize=2;

%Specify Training Options

options = trainingOptions('sgdm', ...
    'ExecutionEnvironment','multi-gpu', ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',2, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imds, ...
    'ValidationFrequency',2);


net = trainNetwork(imds,layers,options);

net = trainNetwork(imds.UnderlyingDatastores,layers,options);


imds.UnderlyingDatastores
% 
% 
% fid = fopen(list{1});
% slice_2d = fread(fid,[512 512], 'float32','l') ;
% % data=transpose(data);
% fclose(fid);
% 
% figure; image(slice_2d);
% 
% 
% 
% 
% % fid = fopen(list{1});
% fid = fopen('~/Dropbox/extras/storage/graham/ht/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample5051/half_lin_cutoff_half_tot_pert_nvpw_v0p6/data/1lf_0.5rf/NSIDE_4/anglid_1/-43-113--256pv_0.20448--0.62099-0.7854ra/2dproj/dm/_1_2dproj_z3_data_slAll.bin');
% slice_3d = fread(fid,512*512*32, 'float32','l') ;  %Optimize this
% slice_3d = reshape(slice_3d,512,512,32);
% % data=transpose(data);
% fclose(fid);
% 
% figure; image(slice_3d(:,:,1));

% slice = 1;
% pos1 = 
% fid = fopen(list{1});
% slice_2d = fread(fid,[512 512], 'float32','l') ;
% % data=transpose(data);
% fclose(fid);

end


