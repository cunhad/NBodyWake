classdef MySequenceDatastore < matlab.io.Datastore & ...
                       matlab.io.datastore.MiniBatchable
    
    properties
        Datastore
        Labels
        NumClasses
        SequenceDimension
        MiniBatchSize
    end
    
    properties(SetAccess = protected)
        NumObservations
    end

    properties(Access = private)
        % This property is inherited from Datastore
        CurrentFileIndex
        FileSet matlab.io.datastore.DsFileSet
    end


    methods
        
        function ds = MySequenceDatastore(location)
            % Construct a MySequenceDatastore object

            % Create a file datastore. The readSequence function is
            % defined following the class definition.
%             ds.FileSet = matlab.io.datastore.ImageDatastore(location,'FileExtensions','.bin',...
%                 'ReadFcn',@readSequence);
%             ds.CurrentFileIndex = 1;
            fds = fileDatastore(location, ...
                'ReadFcn',@readSequence, ...
                'IncludeSubfolders',true);
            ds.Datastore = fds;

            % Read labels from folder names
            numObservations = numel(fds.Files);
            for i = 1:numObservations
                file = fds.Files{i};
                filepath = fileparts(file);
                [~,label] = fileparts(filepath);
                labels{i,1} = label;
            end
            ds.Labels = categorical(labels);
            ds.NumClasses = numel(unique(labels));
            
            % Determine sequence dimension. When you define the LSTM
            % network architecture, you can use this property to
            % specify the input size of the sequenceInputLayer.
            X = preview(fds);
            ds.SequenceDimension = size(X,1);
            
            % Initialize datastore properties.
            ds.MiniBatchSize = 128;
            ds.NumObservations = numObservations;
            ds.CurrentFileIndex = 1;
        end

        function tf = hasdata(ds)
            % Return true if more data is available
            tf = ds.CurrentFileIndex + ds.MiniBatchSize - 1 ...
                <= ds.NumObservations;
        end

        function [data,info] = read(ds)            
            % Read one mini-batch batch of data
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
            % Reset to the start of the data
            reset(ds.Datastore);
            ds.CurrentFileIndex = 1;
        end
        
    end 

    methods (Hidden = true)

        function frac = progress(ds)
            % Determine percentage of data read from datastore
            frac = (ds.CurrentFileIndex - 1) / ds.NumObservations;
        end

    end

end % end class definition

function data = readSequence(filename)
% data = readSequence(filename) reads the sequence X from the MAT-file
% filename

fid = fopen(filename);
slice_3d = fread(fid,512*512*32, 'float32','l') ;  %Optimize this, so only the needed data is loaded
% data=transpose(data);
fclose(fid);
slice_3d = reshape(slice_3d,512,512,32);
data = slice_3d(:,:,1);

end