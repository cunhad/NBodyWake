%vdsrImagePatchDatastore   Super-resolution image datastore for training VDSR
%
%   A vdsrImagePatchDatastore object encapsulates a datastore which
%   creates batches of upsampled low-resolution patches and corresponding
%   residual image patches to be fed to a super-resolution deep neural
%   network for training.
%
%   vdsrImagePatchDatastore properties:
%       MiniBatchSize           - Number of patches in a mini-batch
%       BatchesPerImage         - Number of batches of patches to be extracted per image
%       PatchSize               - Size of the image patches
%       ScaleFactor             - Super-resolution scale factors used for training
%
%   vdsrImagePatchDatastore methods:
%       vdsrImagePatchDatastore - Construct a vdsrImagePatchDatastore
%
%
% Example
%
% -------
%
% folders = fullfile(matlabroot,'toolbox','matlab',{'demos','imagesci'});
%
% exts = {'.jpg','.png','.tif'};
%
% imds = imageDatastore(folders,'FileExtensions',exts);
%
% ds = vdsrImagePatchDatastore(imds,...
%     'MiniBatchSize',64,...
%     'PatchSize', 41,...
%     'BatchesPerImage',1,...
%     'ScaleFactor',4);
%
% miniBatch = read(ds);
%
% figure
% montage(miniBatch.upsampledPatches)
% title('Upsampled Low-Resolution Patches');
% figure
% montage(miniBatch.residualPatches,'DisplayRange',[0 0.05])
% title('Residual Patches (Amplified)');

%   Copyright 2017-2018 The MathWorks, Inc.

% Reference
% ---------
%
% Kim, J., J. K. Lee, and K. M. Lee. "Accurate Image
% Super-Resolution Using Very Deep Convolutional Networks." Proceedings
% of the IEEE Conference on Computer Vision and Pattern Recognition. 2016,
% pp. 1646-1654.

classdef vdsrImagePatchDatastore < matlab.io.Datastore &...
        matlab.io.datastore.MiniBatchable &...
        matlab.io.datastore.Shuffleable
    
    % Required property definition from MiniBatchable mixin
    properties
        MiniBatchSize
    end
    
    % Required property definition from MiniBatchable mixin
    properties (SetAccess = protected)
        NumObservations
    end
    
    properties (Access = private, Hidden, Dependent)
        TotalNumberOfMiniBatches
    end
    
    properties (Access = private)
        InputImageDatastore
        BatchesPerImage
        PatchSize
        CurrentFullImage
        CurrentImageIndex
        CurrentMiniBatchIndex
        NumBatchesReadFromCurrentImage
        ScaleFactor
    end
    
    methods
        
        function ds = vdsrImagePatchDatastore(imds,varargin)
            %
            %   store = vdsrImagePatchDatastore(imds) creates a randomly
            %            cropped upsampled low resolution image and
            %            residual images image patch pair Datastore using 
            %            images from ImageDatastore imds.
            %
            %   store = vdsrImagePatchDatastore(__, Name, Value,__) creates a
            %            randomly cropped pristine and noisy patch pair image
            %            Datastore with additional parameters controlling the data
            %            generation process.
            %
            %   Parameters are:
            %
            %   MiniBatchSize             : Integer specifying the size of
            %                               the mini-batch.
            %                               Default is 64.
            %
            %   BatchesPerImage           : Integer specifying the number
            %                               of batches generated from an
            %                               image.
            %                               Default is 1.
            %
            %   PatchSize                 : Size of the random crops. It
            %                               can be an integer scalar
            %                               specifying same row and column
            %                               sizes or a two element integer
            %                               vector specifying different row
            %                               and column sizes.
            %                               Default is 41.
            %
            %   ScaleFactor               : Vector specifying the scale
            %                               factors used for training.
            %                               Default is 4.
            %
            %   NOTE: This function requires the Deep Learning Toolbox.
            %
            %   Class Support
            %   -------------
            %
            %   imds is an ImageDatastore.
            %
            %   Example
            %   -------
            %
            %   imds = imageDatastore(pathToGrayscaleNaturalImageData);
            %
            %   store = vdsrImagePatchDatastore(imds,...
            %       'MiniBatchSize',64,...
            %       'PatchSize', 41,...
            %       'BatchesPerImage',1,...
            %       'ScaleFactor',3);
            %
            %   layers = VDSRLayers();
            %
            %   opts = trainingOptions('sgdm');
            %
            %   net = trainNetwork(store,layers,opts);
            
              
            options = parseInputs(varargin{:});
            ds.BatchesPerImage = options.BatchesPerImage;
            ds.MiniBatchSize = options.MiniBatchSize;
            ds.NumObservations = length(imds.Files) * ds.MiniBatchSize * ds.BatchesPerImage;
            
            if isscalar(options.PatchSize)
                ds.PatchSize = [options.PatchSize options.PatchSize];
            else
                ds.PatchSize = options.PatchSize;
            end
            
            ds.ScaleFactor = options.ScaleFactor;
            
            % Copy datastore to preserve state of imds input.
            ds.InputImageDatastore = imds.copy();
            
            ds.reset();     
        end
        
    end
    
    % Required method definitions from matlab.io.Datastore
    methods
        
        function [miniBatchTable,info] = read(ds)
            
            if (ds.NumBatchesReadFromCurrentImage == ds.BatchesPerImage)
                ds.NumBatchesReadFromCurrentImage = 0;
                ds.CurrentImageIndex = ds.CurrentImageIndex + 1;
                ds.CurrentFullImage = ds.InputImageDatastore.readimage(ds.CurrentImageIndex);
                
                if numel(size(ds.CurrentFullImage)) ~= 3 || size(ds.CurrentFullImage,3) ~=3
                    error('Expected input image to be in RGB colorspace.');
                end
                
            end
            
            % Use only the luminance component for training
            I = rgb2ycbcr(ds.CurrentFullImage);
            I = im2double(I(:,:,1));            
            
            % Randomly apply one value from ScaleFactor
            if isvector(ds.ScaleFactor)
                scaleFactor = ds.ScaleFactor(randi([1 numel(ds.ScaleFactor)],1));
            else
                scaleFactor = ds.ScaleFactor;
            end
            
            lowresImage = imresize(I,1/scaleFactor,'bicubic');
            upsampledImage = imresize(lowresImage,[size(I,1) size(I,2)],'bicubic');
            
            upsampledPatches = cell(ds.MiniBatchSize,1);
            residualPatches = cell(ds.MiniBatchSize,1);

            % Define image augmentation options
            augOptions = imageDataAugmenter( ...
                'RandXReflection',true, ...
                'RandYReflection',true, ...
                'RandRotation',@()90*randi([0 1],1));
                
            for i = 1:ds.MiniBatchSize

                % Randomly crop the image. The vector roi specifies the size
                % and position of the crop rectangle as [xmin ymin width height]
                cropWidth = ds.PatchSize(1);
                cropHeight = ds.PatchSize(2);
                patchStartCol = randi([1,size(I,2)-cropWidth+1],1);
                patchStartRow = randi([1,size(I,1)-cropHeight+1],1);
                roi = [patchStartCol patchStartRow cropWidth cropHeight];
                
                referencePatch = imcrop(I,roi);
                upsampledPatch = imcrop(upsampledImage,roi);
                residualPatch = referencePatch-upsampledPatch;
                
                % Apply random image augmentation        
                aug = augment(augOptions,{upsampledPatch,residualPatch});      
                upsampledPatches{i} = aug{1};
                residualPatches{i} = aug{2};
                
            end
            
            ds.NumBatchesReadFromCurrentImage = ds.NumBatchesReadFromCurrentImage + 1;
            ds.CurrentMiniBatchIndex = ds.CurrentMiniBatchIndex + 1;
            miniBatchTable = [table(upsampledPatches) table(residualPatches)];
            info.CurrentImageIndexFromDatastore = ds.CurrentImageIndex;
            info.BatchIndexFromCurrentImage = ds.NumBatchesReadFromCurrentImage;
            
        end
        
        function TF = hasdata(ds)
            
            outOfData = (ds.CurrentImageIndex >= length(ds.InputImageDatastore.Files)) &&...
                (ds.NumBatchesReadFromCurrentImage >= ds.BatchesPerImage);
            
            TF = ~outOfData;
            
        end
        
        function reset(ds)
            
            ds.CurrentFullImage = ds.InputImageDatastore.readimage(1);
            ds.CurrentImageIndex = 1;
            ds.NumBatchesReadFromCurrentImage = 0;
            ds.CurrentMiniBatchIndex = 0;
            
        end
        
        function tnmb = get.TotalNumberOfMiniBatches(ds)
            
            tnmb = floor(ds.NumObservations/ds.MiniBatchSize) + ...
                (mod(ds.NumObservations, ds.MiniBatchSize) > 0)*1;
            
        end
        
    end
    
    methods (Hidden)
        function frac = progress(ds)
            frac = ds.CurrentMiniBatchIndex / ds.TotalNumberOfMiniBatches;
        end
    end
    
    % Required method definitions from matlab.io.datastore.Shuffleable
    % mixin.
    methods
        
        function dsrand = shuffle(ds)
            
            % To shuffle, shuffle underlying ImageDatastore
            dsrand = copy(ds);
            dsrand.InputImageDatastore = shuffle(dsrand.InputImageDatastore);
            
        end
    
    % End of required method definitions from matlab.io.datastore.Shuffleable
    % mixin.
    end
    
    
    methods(Static, Hidden = true)
        function self = loadobj(S)
            self = vdsrImagePatchDatastore(S.imds, ...
                'MiniBatchSize', S.MiniBatchSize,...
                'BatchesPerImage', S.BatchesPerImage,...
                'PatchSize', [S.PatchSize(1) S.PatchSize(2)], ...
                'ScaleFactor', S.ScaleFactor);
        end
    end
    
    methods (Hidden)
        function S = saveobj(self)
            
            % Serialize vdsrImagePatchDatastore object
            S = struct('imds',self.imds,...
                'MiniBatchSize',self.MiniBatchSize,...
                'BatchesPerImage',self.BatchesPerImage,...
                'PatchSize',self.PatchSize,...
                'ScaleFactor', self.ScaleFactor);
        end
        
    end
    
end

function B = validateImagedatastore(ds)

validateattributes(ds, {'matlab.io.datastore.ImageDatastore'}, ...
    {'nonempty','vector'}, mfilename, 'IMDS');
validateattributes(ds.Files, {'cell'}, {'nonempty'}, mfilename, 'IMDS');

B = true;

end

function options = parseInputs(varargin)

parser = inputParser();
parser.addParameter('BatchesPerImage',1,@validateBatchesPerImage);
parser.addParameter('PatchSize',41,@validatePatchSize);
parser.addParameter('MiniBatchSize',64,@validateMiniBatchSize);
parser.addParameter('ScaleFactor',4,@validateResizeFactors);

parser.parse(varargin{:});
options = parser.Results;

end


function B = validateResizeFactors(resizeFactors)

attributes = {'nonempty','real','vector', ...
    'positive','integer','finite','nonsparse','nonnan','nonzero'};

validateattributes(resizeFactors,images.internal.iptnumerictypes, attributes,...
    mfilename,'ScaleFactor');

B = true;

end

function B = validateBatchesPerImage(BatchesPerImage)

attributes = {'nonempty','real','scalar', ...
    'positive','integer','finite','nonsparse','nonnan','nonzero'};

validateattributes(BatchesPerImage,images.internal.iptnumerictypes, attributes,...
    mfilename,'BatchesPerImage');

B = true;

end


function B = validateMiniBatchSize(miniBatchSize)

attributes = {'nonempty','real','scalar', ...
    'positive','integer','finite','nonsparse','nonnan','nonzero'};

validateattributes(miniBatchSize,images.internal.iptnumerictypes, attributes,...
    mfilename,'MiniBatchSize');

B = true;

end

function B = validatePatchSize(PatchSize)

attributes = {'nonempty','real','vector', ...
    'positive','integer','finite','nonsparse','nonnan','nonzero'};

validateattributes(PatchSize,images.internal.iptnumerictypes, attributes,...
    mfilename,'PatchSize');

if numel(PatchSize) > 2
    error('Invalid PatchSize');
end

B = true;

end

