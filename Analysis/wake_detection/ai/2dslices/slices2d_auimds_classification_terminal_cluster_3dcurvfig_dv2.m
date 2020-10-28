function [  ] = slices2d_auimds_classification_terminal_cluster_3dcurvfig_dv2()
%UNTITLED Summary of this function goes here
%   we will also classify the non-parallel to wake data



% filename='_col_2dproj_z3_data_sl';
% filename='_col_2dproj_3dcurv_z3_visual_sl';
filename='_2dproj_3dcurv_z3_visual_sl';

% nc=512;
NSIDE=4;
% trsh=20;
% cut=1;
% lev=2;
% sigma = 5;
slices=32;
% anal_lev=2;


N_angles=12*NSIDE*NSIDE/2;
N_angles_t=12*NSIDE*NSIDE;

specs_path_list_nowake='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_nowakem'
% specs_path_list_nowake='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpx_2d_NSIDE4_ai_fig/4Mpc_2048c_1024p_zi63_nowakem'
% specs_path_list_nowake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_nowakem'
sample_list_nowake=dir(strcat(specs_path_list_nowake,'/sample*'));
sample_list_nowake={sample_list_nowake.name};
% sample_list_nowake=sort_nat(sample_list_nowake)

specs_path_list_wake='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpxNSIDE4_wake_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
% specs_path_list_wake='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
% specs_path_list_wake='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_hpx_2d_NSIDE4_wake_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
% specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_wake_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
% specs_path_list_wake='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_hpx_2d_NSIDE4_wake_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
sample_list_wake=dir(strcat(specs_path_list_wake,'/sample*'));
sample_list_wake={sample_list_wake.name};
% sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw');
sample_list_wake=strcat(sample_list_wake,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');
% sample_list_wake=sort_nat(sample_list_wake)


specs_path_list_wake_nonparal='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_cps32_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
% specs_path_list_wake_nonparal='/home/asus/Dropbox/extras/storage/graham/ht/data_cps32_512_hpxNSIDE4_2dclaral1lr1na1024_and_3dparclaral1lr1_visual_fig/4Mpc_2048c_1024p_zi63_wakeGmu4t10m8zi10m'
sample_list_wake_nonparal=dir(strcat(specs_path_list_wake_nonparal,'/sample*'));
sample_list_wake_nonparal={sample_list_wake_nonparal.name};
sample_list_wake_nonparal=strcat(sample_list_wake_nonparal,'/half_lin_cutoff_half_tot_pert_nvpw_v0p6');


sample_id_range=[1 : length(sample_list_nowake)];

% label_numbers_path='/home/cunhad/projects/rrg-rhb/cunhad/simulations/cubep3m/ht/data_hpx_2a3d_anali_all_group/';
% label_number=dlmread(strcat(label_numbers_path,'curvfilt2a3d_label_large.txt'));
% label_numbers_statistics=dlmread(strcat(label_numbers_path,'curvfilt2a3d_label_large_numbers_statistics.txt'));

% 
% label_number=dlmread(strcat(label_numbers_path,'curvfilt2a3d_label_notover.txt'));
% label_numbers_statistics=dlmread(strcat(label_numbers_path,'curvfilt2a3d_label_notover_numbers_statistics.txt'));

count=1;

for w_nw=1:3
% for w_nw=1:1
% for w_nw=1
    
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
    
    if w_nw==3
        count=1;
        specs_path_list=specs_path_list_wake_nonparal;
        sample_list=sample_list_wake_nonparal;
    end
    
    
    for sample = 1:length(sample_id_range)
%     for sample = 1:1
%     sample = 1        
%         map_3d_slices=zeros(nc,nc,slices);
%         map_3d_slices_filt2d=zeros(nc,nc,slices);

%         list_of_angle_paths=dir(char(strcat(specs_path_list,'/',string(sample_list(sample)),'/data/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/anglid_*')));
        list_of_angle_paths=dir(char(strcat(specs_path_list,'/',string(sample_list(sample)),'/visual_3dfilt/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/anglid_*')));
        list_of_angle_paths={list_of_angle_paths.name};
        
        for angle_id=1:length(list_of_angle_paths)
%             for angle_id=1:100
%             angle_id=1
            
            for slice_id=1:slices
                %         for slice_id=1:1
%                 slice_id=1


%                 sample_id=(slices*(sample-1))+slice_id;
%                 
%                 filename_nowake=strcat('',specs_path_list,'/',string(sample_list(sample)),'/data/1lf_1rf_0-0-0pv_1.5708-0-0ra/2dproj/dm/2d_curvfilt/',ch,filename,num2str(slice_id),'_curvfilt2a3d','.bin');                
%                 filename_nowake=char(filename_nowake);
%                 %             %             fid = fopen(filename_nowake);
%                 %             %         scalefactor = fread(fid, [1 1], 'float32','l') ;
%                 % %             display(filename_nowake)
%                 %             slice_2d_ds = fileDatastore(filename_nowake,'ReadFcn',@read_slices_bin,'FileExtensions','.bin');
%                 %             slice_2d=cell2mat(tall(slice_2d_ds));
%                 %
%                 %             figure; imagesc([2/1024:4/1024:4],[2/1024:4/1024:4],gather(log(slice_2d))); colorbar; axis('image');
%                 %             xlabel('$Z(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 %             ylabel('$Y(Mpc/h)$', 'interpreter', 'latex', 'fontsize', 20);
%                 %             set(gca,'FontName','FixedWidth');
%                 %             set(gca,'FontSize',16);
%                 %             set(gca,'linewidth',2);
%                 
%                 list{count}=filename_nowake;
%                 count=count+1;
%                 %             fclose(fid);
                
                
                path1=strcat(specs_path_list,'/',string(sample_list(sample)),'/visual_3dfilt/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/');
                
                path2=dir(char(strcat(path1,"*pv*")));
                path2={path2.name};
                path2=path2(1);
                
                path3='/2dproj/dm/';
                
%                 path_in=strcat(specs_path_list,'/',string(sample_list(sample)),'/visual_2dfilt/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/',string(path2),'/2dproj/dm/',ch,filename,num2str(slice_id),'_log_fig.png');
                path_in=strcat(specs_path_list,'/',string(sample_list(sample)),'/visual_3dfilt/1lf_0.5rf','/NSIDE_',num2str(NSIDE),'/',list_of_angle_paths(angle_id),'/',string(path2),'/2dproj/dm/',ch,filename,num2str(slice_id),'.png');
                
                if w_nw==3
                    list_nonparal{count}=char(path_in);
                    count=count+1;
                else
                    list{count}=char(path_in);
                    count=count+1;
                end
                
            end
            
        end
        
    end
    
    
    
end

%labels

label= categorical(abs(double(contains(string(list),'nowake'))-1));
label_logical= logical(abs(double(contains(string(list),'nowake'))-1));
label_nonparal= categorical(abs(double(contains(string(list_nonparal),'nowake'))-1));
% label_num
% list_nowake=list(~label_logical);
% list_wake=list(label_logical);


imds = imageDatastore(list,'Labels',label);
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.6,'randomized');
imdsValidation_nonparal = imageDatastore(list_nonparal,'Labels',label_nonparal);

labelCount = countEachLabel(imds);
img = readimage(imds,1);
size(img);

% figure; imshow(readimage(imds,1));

% Augmenter
% 
% imageAugmenter = imageDataAugmenter( ...
%     'RandRotation',[0,360], ...
%     'RandXReflection',1,...
%     'RandYReflection',1);
% %     'RandScale',[2 2],...
% %     'RandXTranslation',[]);
% 
% % auimdsTrain = augmentedImageDatastore([1024 1024],imdsTrain,'DataAugmentation',imageAugmenter);
% % auimdsValidation = augmentedImageDatastore([1024 1024],imdsValidation,'DataAugmentation',imageAugmenter);
% 
% auimdsTrain = augmentedImageDatastore([224 224],imdsTrain,'ColorPreprocessing','gray2rgb');
% auimdsValidation = augmentedImageDatastore([224 224],imdsValidation,'ColorPreprocessing','gray2rgb');


% %display data
% 
% data = readByIndex(auimdsTrain,2);
% figure; imshow(cell2mat(data{1:1,{'input'}}));

% Define the convolutional neural network architecture
% 
% layers = [
% %     imageInputLayer([1024 1024 3])
%     imageInputLayer([1065 1065 3])
%     
%     convolution2dLayer(3,8,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
%     maxPooling2dLayer(2,'Stride',2)
%     
%     convolution2dLayer(3,16,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
%     maxPooling2dLayer(2,'Stride',2)
%     
%     convolution2dLayer(3,32,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
%     convolution2dLayer(3,32,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
%     fullyConnectedLayer(2)
%     softmaxLayer
%     classificationLayer];



layers = [
% imageInputLayer([1024 1024 3])
imageInputLayer([1065 1065 3])

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
convolution2dLayer(4,8,'Padding','same')
batchNormalizationLayer
maxPooling2dLayer(2,'Stride',2)
reluLayer
convolution2dLayer(8,8,'Padding','same')
batchNormalizationLayer
reluLayer
convolution2dLayer(8,16,'Padding','same')
batchNormalizationLayer
reluLayer
maxPooling2dLayer(2,'Stride',2)
convolution2dLayer(16,16,'Padding','same')
batchNormalizationLayer
reluLayer
convolution2dLayer(16,32,'Padding','same')
batchNormalizationLayer
reluLayer
maxPooling2dLayer(2,'Stride',2)
convolution2dLayer(32,32,'Padding','same')
batchNormalizationLayer
reluLayer
convolution2dLayer(32,64,'Padding','same')
batchNormalizationLayer
reluLayer
maxPooling2dLayer(2,'Stride',2)
% convolution2dLayer(64,64,'Padding','same')
% batchNormalizationLayer
% reluLayer
% convolution2dLayer(64,128,'Padding','same')
% batchNormalizationLayer
% reluLayer
% maxPooling2dLayer(2,'Stride',2)
fullyConnectedLayer(256)
fullyConnectedLayer(256)
fullyConnectedLayer(64)
fullyConnectedLayer(2)
softmaxLayer
classificationLayer];


%minibatch size

miniBatchSize=32;

%Specify Training Options

options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',16, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',300);
%     'Verbose',false, ...
%     'Plots','training-progress');
%Train Network Using Training Data


net = trainNetwork(imdsTrain,layers,options);

YPred = classify(net,imdsValidation,'MiniBatchSize',miniBatchSize);
YValidation = imdsValidation.Labels;
accuracy = sum(YPred == YValidation)/numel(YValidation)

YPred_wake=YPred;
FractionPredic_with_wake=sum(YPred_wake==categorical(1))/numel(YPred_wake)

YPred_wake_where_theis=categorical(str2num(char(YPred_wake)).*str2num(char(YValidation)));
FractionYPred_wake_where_thereis=sum(YPred_wake_where_theis==categorical(1))/sum(YPred_wake==categorical(1))

YPred_wake_where_theisnot=categorical(str2num(char(YPred_wake)).*abs(1-str2num(char(YValidation))));
FractionYPred_wake_where_thereisnot=sum(YPred_wake_where_theisnot==categorical(1))/sum(YPred_wake==categorical(1))


YValidation_wake=YValidation(YValidation==categorical(1));



YPred_nonparal = classify(net,imdsValidation_nonparal,'MiniBatchSize',miniBatchSize);
YValidation_nonparal = imdsValidation_nonparal.Labels;
accuracy_nonparal = sum(YPred_nonparal == YValidation_nonparal)/numel(YValidation_nonparal)
YPred_wake_nonparal=YPred_nonparal;
FractionPredic_with_wake_nonparal=sum(YPred_wake_nonparal==categorical(1))/numel(YPred_wake_nonparal)
YPred_wake_where_theis_nonparal=categorical(str2num(char(YPred_wake_nonparal)).*str2num(char(YValidation_nonparal)));
FractionYPred_wake_where_thereis_nonparal=sum(YPred_wake_where_theis_nonparal==categorical(1))/sum(YPred_wake_nonparal==categorical(1))



% 
% 
% net = googlenet;
% % analyzeNetwork(net);
% % net.Layers(1)
% inputSize = net.Layers(1).InputSize;
% 
% if isa(net,'SeriesNetwork') 
%   lgraph = layerGraph(net.Layers); 
% else
%   lgraph = layerGraph(net);
% end 
% 
% edit(fullfile(matlabroot,'examples','nnet','main','findLayersToReplace.m'))
% 
% [learnableLayer,classLayer] = findLayersToReplace(lgraph);
% [learnableLayer,classLayer] 
% 
% 
% numClasses = numel(categories(imdsTrain.Labels));
% 
% if isa(learnableLayer,'nnet.cnn.layer.FullyConnectedLayer')
%     newLearnableLayer = fullyConnectedLayer(numClasses, ...
%         'Name','new_fc', ...
%         'WeightLearnRateFactor',10, ...
%         'BiasLearnRateFactor',10);
%     
% elseif isa(learnableLayer,'nnet.cnn.layer.Convolution2DLayer')
%     newLearnableLayer = convolution2dLayer(1,numClasses, ...
%         'Name','new_conv', ...
%         'WeightLearnRateFactor',10, ...
%         'BiasLearnRateFactor',10);
% end
% 
% lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);
% 
% newClassLayer = classificationLayer('Name','new_classoutput');
% lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);
% 
% layers = lgraph.Layers;
% connections = lgraph.Connections;
% 
% layers(1:10) = freezeWeights(layers(1:10));
% lgraph = createLgraphUsingConnections(layers,connections);
% 
% options = trainingOptions('sgdm', ...
%     'MiniBatchSize',2, ...
%     'MaxEpochs',6, ...
%     'InitialLearnRate',0.01, ...
%     'Shuffle','every-epoch', ...
%     'ValidationData',auimdsValidation, ...
%     'ValidationFrequency',30, ...
%     'Verbose',false, ...
%     'Plots','training-progress');
% %     'InitialLearnRate',3e-4, ...
% net = trainNetwork(auimdsTrain,lgraph,options);
% 
% 
%  YPred = classify(net,auimdsValidation,'MiniBatchSize',miniBatchSize);
%  YValidation = imdsValidation.Labels;
%  accuracy = sum(YPred == YValidation)/numel(YValidation)


end

