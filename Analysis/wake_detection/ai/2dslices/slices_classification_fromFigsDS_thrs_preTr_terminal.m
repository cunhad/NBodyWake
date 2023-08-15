function [] = slices_classification_fromFigsDS_thrs_preTr_terminal()



% myCluster = parcluster('local');
% myCluster.NumWorkers=32;
% saveProfile(myCluster);
% 
% % parpool('Processes',32);
% 
% parpool('local',32);


% tic;

% input vars
path_soft_links = '/scratch/cunhad/data_cps32_512_hpx_2d_NSIDE4_figs_thr40';
path_soft_linksTot = '/scratch/cunhad/soft_links_Figs/data_cps32_512_hpx_2d_NSIDE4';
% filename='_1_2dproj_z3_data_sl32All';
% filename='_1_2dproj_z3_data_slAll'; 
% slices_sz=32;
% angle_sz=96;

%analysis vars
percentage_to_train = 80;
miniBatchSize=8; %defout is 128, reduced due to memory



% sample_list_nowake=dir(strcat(path_soft_links,'/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_nowakem/sample*'));
% sample_list_nowake_folder=sample_list_nowake.folder;
% sample_list_nowake={sample_list_nowake.name};
% 
% sample_list_wake=dir(strcat(path_soft_links,'/data_cps32_512_hpx_2d_NSIDE4/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m/sample*'));
% sample_list_wake_folder = sample_list_wake.folder;
% sample_list_wake={sample_list_wake.name};
% 
% sample_list_com = intersect(sample_list_nowake,sample_list_wake, 'stable');
% idx_test = randperm(length(sample_list_com),floor(length(sample_list_com)*percentage_to_test/100)) ;
% 
% 
% 
% list = strings(2,length(sample_list_nowake),angle_sz,slices_sz);
% 
% 
% for w_nw=1:2
%     if w_nw==1
%         sample_list=sample_list_nowake;
%         sample_list_folder = sample_list_nowake_folder;
%     else
% %         sample_list=sample_list_wake;
%         sample_list=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');
%         sample_list_folder = sample_list_wake_folder;
%     end
%     sample_list_sz = length(sample_list); %required for parfor
%     parfor sample_id = 1:sample_list_sz
%         sample = sample_list(sample_id);
%         for angle_id = 1:angle_sz
%             angle = char(strcat('anglid_',num2str(angle_id)));
%             for slice_id=1:slices_sz
%                 file=strcat(sample_list_folder,'/',sample,'/',angle,'/',filename,num2str(slice_id),'.png');
%                 file=char(file);
%                 list(w_nw,sample_id,angle_id,slice_id)=file;
%             end
%         end
%     end
% end
% 
% list = list(:);
% 
% toc;
% 
% list_train = list(contains(string(list'),sample_list_com(idx_test)));
% list_validate = list(~contains(string(list'),sample_list_com(idx_test)));
% 
% label_train= categorical(abs(double(contains(string(list_train),'nowake'))-1));
% label_validate= categorical(abs(double(contains(string(list_validate),'nowake'))-1));
% 
% % imds_train = imageDatastore(list_train,'ReadFcn',@read_slices_bin_slices,'FileExtensions','.bin','Labels',label_train);
% % imds_validate = imageDatastore(list_validate,'ReadFcn',@read_slices_bin_slices,'FileExtensions','.bin','Labels',label_validate);
% 
% imds_train = imageDatastore(list_train,'Labels',label_train);
% imds_validate = imageDatastore(list_validate,'Labels',label_validate);

imds = imageDatastore(path_soft_links,"IncludeSubfolders",true,'LabelSource','foldernames');
%imds = imageDatastore(path_soft_links,"IncludeSubfolders",true);
%imds.Labels = categorical(abs(double(contains(string(imds.Files),'nowake'))-1));
[imds_train, imds_validate] = splitEachLabel(imds, percentage_to_train/100, 'randomize');

% figure; image(readimage(imds_train,1)); colorbar;

% layers = [
% % imageInputLayer([1024 1024 3])
% imageInputLayer([1065 1065 3])
% % imageInputLayer([512 512 1])
% 
% convolution2dLayer(2,1,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(2,2,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% convolution2dLayer(4,4,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% fullyConnectedLayer(2)
% softmaxLayer
% classificationLayer];




net = inceptionresnetv2;
inputSize = net.Layers(1).InputSize;


% deepNetworkDesigner(net)


% Replace Final Layers


lgraph = layerGraph(net);

[learnableLayer,classLayer] = findLayersToReplace(lgraph);

numClasses = numel(categories(imds_train.Labels));

if isa(learnableLayer,'nnet.cnn.layer.FullyConnectedLayer')
    newLearnableLayer = fullyConnectedLayer(numClasses, ...
        'Name','new_fc', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
    
elseif isa(learnableLayer,'nnet.cnn.layer.Convolution2DLayer')
    newLearnableLayer = convolution2dLayer(1,numClasses, ...
        'Name','new_conv', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
end

lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);

newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);


% Freeze Initial Layers


layers = lgraph.Layers;
connections = lgraph.Connections;

layers(1:50) = freezeWeights(layers(1:50));
lgraph = createLgraphUsingConnections(layers,connections);


% Train Network

augimdsTrain = augmentedImageDatastore(inputSize(1:2),imds_train);
augimdsValidation = augmentedImageDatastore(inputSize(1:2),imds_validate);




% layers = [
%     % imageInputLayer([1024 1024 3])
%     imageInputLayer([1065 1065 3])
%     % imageInputLayer([512 512 1])
% 
%     convolution2dLayer(5,4,'Stride',1)
%     batchNormalizationLayer
%     reluLayer
% 
%     convolution2dLayer(5,4,'Stride',1)
%     batchNormalizationLayer
%     reluLayer  
%     
%     convolution2dLayer(5,4,'Stride',1)
%     batchNormalizationLayer
%     reluLayer
% 
%     convolution2dLayer(5,4,'Stride',1)
%     batchNormalizationLayer
%     reluLayer  
%     
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer    
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%     
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer 
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%    
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer 
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%     
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer 
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%     
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer 
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer 
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
% 
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer  
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
% 
%     convolution2dLayer(5,10,'Stride',1)
%     batchNormalizationLayer
%     reluLayer     
% 
%     convolution2dLayer(5,20,'Stride',1)
%     batchNormalizationLayer
%     reluLayer   
% 
%     maxPooling2dLayer(2,'Stride',1)
% 
%     convolution2dLayer(5,4,'Stride',1)
%     batchNormalizationLayer
%     reluLayer
% 
% 
% 
% %     flattenLayer
% 
% %     dropoutLayer
% % 
% %     fullyConnectedLayer(10)
% 
% %     dropoutLayer
% 
%     fullyConnectedLayer(4)
% 
%     dropoutLayer
% 
%     fullyConnectedLayer(2)
%     
%     softmaxLayer
% 
%     classificationLayer];



% layers = [
% % imageInputLayer([1024 1024 3])
% imageInputLayer([1065 1065 3])
% 
% convolution2dLayer(2,1,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(2,2,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% convolution2dLayer(4,4,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(4,8,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% convolution2dLayer(8,8,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(8,16,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% convolution2dLayer(16,16,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(16,32,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% convolution2dLayer(32,32,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% convolution2dLayer(32,64,'Padding','same')
% batchNormalizationLayer
% reluLayer
% 
% maxPooling2dLayer(2,'Stride',2)
% 
% % convolution2dLayer(64,64,'Padding','same')
% % batchNormalizationLayer
% % reluLayer
% % convolution2dLayer(64,128,'Padding','same')
% % batchNormalizationLayer
% % reluLayer
% % maxPooling2dLayer(2,'Stride',2)
% 
% fullyConnectedLayer(256)
% fullyConnectedLayer(256)
% fullyConnectedLayer(64)
% fullyConnectedLayer(2)
% softmaxLayer
% classificationLayer];

%Specify Training Options

options = trainingOptions('sgdm', ...
    'ExecutionEnvironment','multi-gpu', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',2, ...
    'Shuffle','every-epoch', ...
    'ValidationData',augimdsValidation, ...
    'MiniBatchSize',miniBatchSize, ...
    'ValidationFrequency',5);
%     'Plots','training-progress',...

delete(gcp('nocreate'));

% numGPUs = gpuDeviceCount("available");
% parpool(numGPUs);

% net = trainNetwork(imds_train,layers,options)

net = trainNetwork(augimdsTrain,lgraph,options);  

YPred = classify(net,augimdsValidation,'MiniBatchSize',miniBatchSize);
% YValidation = augimdsValidation.Labels;
[~,info]=read(augimdsValidation)
YValidation = info.Label;
accuracy = sum(YPred == YValidation)/numel(YValidation)


% FractionPredic_with_wake=sum(YPred_wake=='4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m')/numel(YPred)
% sum(YPred_wake=='4Mpc_2048c_1024p_zi63_nowakem')
YPred_categ = (YPred=='4Mpc_2048c_1024p_zi63_nowakem');
% YPred_categ(end)=0;
YValidation_categ = (YValidation=='4Mpc_2048c_1024p_zi63_nowakem');

YPred_wake_where_theis=categorical(double(~YPred_categ).*double(~YValidation_categ));
FractionYPred_wake_where_thereis=sum(YPred_wake_where_theis==categorical(1))/sum(YValidation_categ==logical(1))

YPred_wake_where_theisnot=categorical(double(~YPred_categ).*double(YValidation_categ));
FractionYPred_wake_where_thereisnot=sum(YPred_wake_where_theisnot==categorical(1))/sum(YValidation_categ==logical(0))





imdsTot = imageDatastore(path_soft_linksTot,"IncludeSubfolders",true,'LabelSource','foldernames');
augimdsTot = augmentedImageDatastore(inputSize(1:2),imdsTot);


YPredTot = classify(net,augimdsTot,'MiniBatchSize',miniBatchSize);                                                                                                                                             
% YValidationTot = augimdsTot.Labels;
[~,info]=read(augimdsTot)
YValidationTot = info.Label;
accuracyTot = sum(YPredTot == YValidationTot)/numel(YValidationTot)


% FractionPredic_with_wake=sum(YPred_wake=='4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m')/numel(YPred)
% sum(YPred_wake=='4Mpc_2048c_1024p_zi63_nowakem')
YPredTot_categ = (YPredTot=='4Mpc_2048c_1024p_zi63_nowakem');
% YPred_categ(end)=0;
YValidationTot_categ = (YValidationTot=='4Mpc_2048c_1024p_zi63_nowakem');

YPred_wake_where_theisTot=categorical(double(~YPredTot_categ).*double(~YValidationTot_categ));
FractionYPred_wake_where_thereisTot=sum(YPred_wake_where_theisTot==categorical(1))/sum(YValidationTot_categ==logical(1))

YPred_wake_where_theisnotTot=categorical(double(~YPredTot_categ).*double(YValidationTot_categ));
FractionYPred_wake_where_thereisnotTot=sum(YPred_wake_where_theisnotTot==categorical(1))/sum(YValidationTot_categ==logical(0))

end
