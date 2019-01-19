classdef associatedImageDatastore < matlab.io.Datastore & ...
                       matlab.io.datastore.MiniBatchable & ...
                       matlab.io.datastore.Shuffleable  & ...
                       matlab.io.datastore.PartitionableByIndex
    
    properties
        MiniBatchSize
        OutputImds
    end
    
    properties(SetAccess = protected)
        NumObservations
    end

    properties(Access = private)
        % This property is inherited from Datastore
        CurrentFileIndex
        % These custom properties store copies of the two ImageDatastores
        InputImds
%         OutputImds
    end


    methods
        
        function ds = associatedImageDatastore(inputImds,outputImds,miniBatchSize)
            % Construct an associatedImageDatastore object
            ds.InputImds = copy(inputImds);
            ds.OutputImds = outputImds;
            ds.InputImds.ReadSize = miniBatchSize;
            ds.NumObservations = length(inputImds.Files);
            ds.MiniBatchSize = miniBatchSize;
            ds.CurrentFileIndex = 1;
        end

        function tf = hasdata(ds)
            % Return true if more data is available
            tf = hasdata(ds.InputImds);
        end

%         function [data,info] = read(ds)            
%             % Read one batch of data
%             inputImageData = read(ds.InputImds)
% %             ds.CurrentFileIndex
% %             ds.MiniBatchSize
%             outputImageData = ds.OutputImds(1:16)
% %             outputImageData = ds.OutputImds;
%             size(inputImageData)
%             size(outputImageData)
%             data = table(inputImageData,outputImageData);            
%             info.batchSize = size(data,1);
%             ds.CurrentFileIndex = ds.CurrentFileIndex + info.batchSize;
%             info.currentFileIndex = ds.CurrentFileIndex;  
%         end

function [data,info] = read(ds)
            % [data,info] = read(ds) read one mini-batch of data.
            
            miniBatchSize = ds.MiniBatchSize;
            info = struct;
            
            for i = 1:miniBatchSize
                predictors{i,1} = read(ds.InputImds);
                responses(i,1) = ds.OutputImds(ds.CurrentFileIndex);
                ds.CurrentFileIndex = ds.CurrentFileIndex + 1;
            end
            
            data = preprocessData(ds,predictors,responses);
        end
        
        function data = preprocessData(ds,predictors,responses)
            % data = preprocessData(ds,predictors,responses) preprocesses
            % the data in predictors and responses and returns the table
            % data
            
            miniBatchSize = ds.MiniBatchSize;
            
            % Pad data to length of longest sequence.
            sequenceLengths = cellfun(@(X) size(X,2),predictors);
            maxSequenceLength = max(sequenceLengths);
            for i = 1:miniBatchSize
                X = predictors{i};
                
                % Pad sequence with zeros.
                if size(X,2) < maxSequenceLength
                    X(:,maxSequenceLength) = 0;
                end
                
                predictors{i} = X;
            end
            
            % Return data as a table.
            data = table(predictors,responses);
        end

        function reset(ds)
            % Reset to the start of the data
            reset(ds.InputImds);
%             ds.OutputImds=[];
            ds.CurrentFileIndex = 1;
        end
        
        function dsnew = shuffle(ds)
            dsnew = copy(ds);
            shuffledIndexOrder = randperm(ds.NumObservations);
            dsnew.InputImds.Files = dsnew.InputImds.Files(shuffledIndexOrder);
            dsnew.OutputImds.Files = dsnew.OutputImds.Files(shuffledIndexOrder);
        end
        
        function dsnew = partitionByIndex(ds,indices)  
           dsnew = copy(ds);
           dsnew.InputImds.Files = dsnew.InputImds.Files(indices);
           dsnew.OutputImds.Files = dsnew.OutputImds.Files(indices);
        end
        
    end 

    methods (Hidden = true)

        function frac = progress(ds)
            % Determine percentage of data read from datastore
            frac = (ds.CurrentFileIndex-1)/ds.NumObservations;
        end

    end

end % end class definition