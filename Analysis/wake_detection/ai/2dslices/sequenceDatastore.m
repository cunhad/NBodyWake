classdef sequenceDatastore < matlab.io.Datastore & ...
        matlab.io.datastore.MiniBatchable & ...
        matlab.io.datastore.Shuffleable  & ...
        matlab.io.datastore.PartitionableByIndex
    
    properties
        Datastore
        Labels
        NumClasses
        SequenceDimension
        MiniBatchSize
        list
    end
    
    properties(SetAccess = protected)
        NumObservations
    end
    
    properties(Access = private)
        CurrentFileIndex
        InputImds
    end
    
    methods
        function ds = sequenceDatastore(folder)
            % ds = sequenceDatastore(folder) creates a sequence datastore
            % from the data in folder.
            
            ds.list=folder;
            
            % Create file datastore.
            fds = fileDatastore(ds.list, ...
                'ReadFcn',@readSequence, ...
                'IncludeSubfolders',true);
            ds.Datastore = fds;
            
            ds.InputImds = copy(ds);
            
            % Read labels from folder names.
            numObservations = numel(fds.Files);
            for i = 1:numObservations
                file = fds.Files{i};
%                 filepath = fileparts(file);
%                 [~,label] = fileparts(filepath);
                if (contains(file,'nowake'))                    
                labels{i,1} = 0;
                else
                labels{i,1} = 1;
                end
            end
            ds.Labels = labels;
            ds.NumClasses = numel(unique(char(string(labels))));
            
            % Determine sequence dimension.
            X = preview(fds);
%             ds.SequenceDimension = size(X,1);
            ds.SequenceDimension = size(X);
            % Initialize datastore properties.
            ds.MiniBatchSize = 128;
            ds.NumObservations = numObservations;
            ds.CurrentFileIndex = 1;
        end
        
        function tf = hasdata(ds)
            % tf = hasdata(ds) returns true if more data is available.
            
            tf = ds.CurrentFileIndex + ds.MiniBatchSize - 1 ...
                <= ds.NumObservations;
        end
        
        function [data,info] = read(ds)
            % [data,info] = read(ds) read one mini-batch of data.
            
            miniBatchSize = ds.MiniBatchSize;
            info = struct;
            
            for i = 1:miniBatchSize
                predictors{i,1} = read(ds.Datastore);
                responses(i,1) = ds.Labels(ds.CurrentFileIndex);
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
            % reset(ds) resets the datastore to the start of the data.
            
            reset(ds.Datastore);
            ds.CurrentFileIndex = 1;
        end
        
        function dsNew = shuffle(ds)
            % dsNew = shuffle(ds) shuffles the files and the corresponding
            % labels in the datastore.
            
            % Create copy of datastore.
            dsNew = copy(ds);
            dsNew.Datastore = copy(ds.Datastore);
            fds = dsNew.Datastore;
            
            % Shuffle files and corresponding labels.
            numObservations = dsNew.NumObservations;
            idx = randperm(numObservations);
            fds.Files = fds.Files(idx);
            dsNew.Labels = dsNew.Labels(idx);
        end
        
        function dsnew = partitionByIndex(ds,indices)  
%            dsnew = copy(ds);
%            dsnew.InputImds.Files = dsnew.InputImds.Files(indices);
%            dsnew.OutputImds.Files = dsnew.OutputImds.Files(indices);
% 
            dsnew = copy(ds);
            dsnew.Datastore = copy(ds.Datastore);
            fds = dsnew.Datastore;
            
            % Shuffle files and corresponding labels.
%             numObservations = dsnew.NumObservations;
%             idx = randperm(numObservations);
indices
            (size(fds.Files))            
            fds.Files = fds.Files(indices);
            dsnew.Labels = dsnew.Labels(indices);
            dsnew.Labels = 1;
            (size(fds.Files))
            
            
        end
    end
    
    methods (Hidden = true)
        function frac = progress(ds)
            % frac = progress(ds) returns the percentage of observations
            % read in the datastore.
            
            frac = (ds.CurrentFileIndex - 1) / ds.NumObservations;
        end
    end
end

function data = readSequence(filename)
% data = readSequence(filename) reads the sequence X from the MAT file
% filename
fid = fopen(filename);
data=zeros(1024,1024,1);
data(:,:,1) = fread(fid,[1024 1024], 'float32','l') ;
% data=transpose(data);
fclose(fid);
end