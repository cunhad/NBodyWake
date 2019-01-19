function [  ] = slices2d_ds_classification_curv()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% slices2d_ds( root,root_data_2d_in,root_data_2d_anali,spec,aux_path,aux_path_out,filename,lenght_factor,resol_factor,numb_rand,slice,NSIDE ,num_cores)

%(example)  [ map ,anali] = slices2d_ds('/home/asus/Dropbox/extras/storage/graham/small_res/','/home/asus/Dropbox/extras/storage/graham/small_res/data_test2/','/home/asus/Dropbox/extras/storage/graham/small_res/anali/','64Mpc_256c_128p_zi63_nowakem','/sample2001/','','10.000xv0.dat',1,1,1,2,2,1);

 

% myCluster = parcluster('local');
% myCluster.NumWorkers=num_cores;
% saveProfile(myCluster);
% 
% p = parpool(num_cores);




filename='_2dproj_z3_data_sl';
nc=1024;
% trsh=20;
% cut=1;
% lev=2;
% sigma = 5;
slices=32;
% anal_lev=2;

specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2a3dcurvfilt/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2a3dcurvfilt/4Mpc_2048c_1024p_zi63_wakeGmu1t10m7zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw');
% sample_list_wake=sort_nat(sample_list_wake)

sample_id_range=[1 : length(sample_list_nowake)];

count=1;

for w_nw=1:2
% for w_nw=1:1
    
    if w_nw==1
        specs_path_list=specs_path_list_nowake;
        sample_list=sample_list_nowake;
        ch='_7';
        coul='b';
    else
        specs_path_list=specs_path_list_wake;
        sample_list=sample_list_wake;
        ch='_4';
        coul='r';
    end
    
    
        for sample = 1:length(sample_id_range)
%     for sample = 1:1
        
%         map_3d_slices=zeros(nc,nc,slices);
%         map_3d_slices_filt2d=zeros(nc,nc,slices);
        
                for slice_id=1:slices
%         for slice_id=1:1
            
            sample_id=(slices*(sample-1))+slice_id;
            
            filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/2d_curvfilt/',ch,filename,num2str(slice_id),'_curvfilt2a3d','.bin');
            
            filename_nowake=char(filename_nowake);
%             %             fid = fopen(filename_nowake);
%             %         scalefactor = fread(fid, [1 1], 'float32','l') ;
% %             display(filename_nowake)
%             slice_2d_ds = fileDatastore(filename_nowake,'ReadFcn',@read_slices_bin,'FileExtensions','.bin');
%             slice_2d=cell2mat(tall(slice_2d_ds));
%             
%             figure; imagesc([2/1024:4/1024:4],[2/1024:4/1024:4],gather(log(slice_2d))); colorbar; axis('image');
%             xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%             ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%             set(gca,'FontName','FixedWidth');
%             set(gca,'FontSize',16);
%             set(gca,'linewidth',2);
            
        list{count}=filename_nowake;
        count=count+1;
            %             fclose(fid);
            
        end
        
    end
    
    
    
end


%labels

label_numbers_path='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_1024_2a3dcurvfilt_all/';
label_number=dlmread(strcat(label_numbers_path,'curvfilt2a3d_label.txt'));

sum_label_number=sum(label_number);

label= categorical(abs(double(label_number)));

% label= categorical(abs(double(contains(list,'nowake'))-1));
%  imds = imageDatastore(list,'ReadFcn',@read_slices_bin,'FileExtensions','.bin','Labels',label);

%load datastore

imds = imageDatastore(list,'ReadFcn',@read_slices_bin,'FileExtensions','.bin','Labels',label);
labelCount = countEachLabel(imds);
img = readimage(imds,1);
size(img);


numTrainFiles = ceil(sum_label_number/2);
% numTrainFiles = 320;
% [imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.9,'randomize');


% Define the convolutional neural network architecture

layers = [
    imageInputLayer([1024 1024 1])
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

% see the image

% figure; imshow(readimage(imds,1));





% %or
% 
% inputSize = imdsTrain.SequenceDimension;
% numClasses = imdsTrain.NumClasses;
% numHiddenUnits = 100;
% 
% layers = [ ...
%     sequenceInputLayer(inputSize)
%     lstmLayer(numHiddenUnits,'OutputMode','last')
%     fullyConnectedLayer(numClasses)
%     softmaxLayer
%     classificationLayer];

%minibatch size

miniBatchSize=2;

%Specify Training Options

options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');
%Train Network Using Training Data


net = trainNetwork(imdsTrain,layers,options);

YPred = classify(net,imdsValidation,'MiniBatchSize',miniBatchSize);
YValidation = imdsValidation.Labels;
accuracy = sum(YPred == YValidation)/numel(YValidation)

end




