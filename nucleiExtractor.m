classdef  nucleiExtractor < handle
    
    properties 
        listing;
        numStacks;
        hasNextStack;
        currentStack;
        sourceDir;
        segChannel;
        minVol;
        minFilledVol;
        sepWavelengthsFlag;
        waveLengthsStrings;
    end
    
    methods 
        
        function extractor = nucleiExtractor(...
                sourceDir,segChannel,minVol,minFilledVol,varargin)
            
            % Store instantiation inputs
            extractor.sourceDir = sourceDir;
            extractor.segChannel = segChannel;
            extractor.minVol = minVol;
            extractor.minFilledVol = minFilledVol;
            
			% Deal with separate wavelength input option
            if isempty(varargin)
                extractor.sepWavelengthsFlag = false;
            else
                extractor.sepWavelengthsFlag = true;
                extractor.waveLengthsStrings = varargin{1};
				if isempty(extractor.waveLengthsStrings)
					extractor.sepWavelengthsFlag = false;
				end
            end
                
            % --- Get files to read
            
            % Recursively lists all matching files in the source directory
            extractor.listing = rdir(sourceDir);
            if extractor.sepWavelengthsFlag
                % only contain strings pointing to first wavelength
                pathList = {extractor.listing.name};
                pathList = cellfun(@(xx)sprintf('%s',xx),pathList,...
                    'UniformOutput',false);
                pathList = [pathList(:)];
                extractor.listing = extractor.listing(...
                    cellfun(@(xx)~isempty(xx),...
                    strfind(pathList,varargin{1}{1})));
            end
            
            extractor.numStacks = numel(extractor.listing); % number of source files
            extractor.listing = extractor.listing(...
                randperm(extractor.numStacks)); % randomize stack order

            extractor.currentStack = 1;
            extractor.hasNextStack = ...
                extractor.currentStack<=extractor.numStacks;
            
        end
        
        function [nuclei,voxelSizes] = retrieveNuclei(extractor)
            
            if extractor.hasNextStack
                % If the directory still has stacks left
                
                filepath = extractor.listing(extractor.currentStack).name;
                reader_1 = bfGetReader(filepath);
                omeMeta = reader_1.getMetadataStore();
                
                % --- get the voxel edge sizes
                voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0);
                voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
                rawVoxelSizeX = voxelSizeX.doubleValue();
                voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0);
                voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
                rawVoxelSizeY = voxelSizeY.doubleValue();
                
                voxelSizes = [rawVoxelSizeY,rawVoxelSizeX];
                
                % --- get the spatial stack dimensions
                rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, voxels
                rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, voxels
                rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % image height, voxels
                
                stackDim = [rawStackSizeY,rawStackSizeX,rawStackSizeZ];
                
                numImages = omeMeta.getImageCount();
                
                % get number of color channels
                if ~extractor.sepWavelengthsFlag
                    numChannels = omeMeta.getChannelCount(0);
                else
                    numChannels = numel(extractor.waveLengthsStrings);
                end
                
                if extractor.sepWavelengthsFlag
                    % open additional image reader instances for the other
                    % color channels, in case color channels are stored in
                    % separate files
                    
                    extraReaders = cell(1,numChannels-1);
                    
                    for cc = 2:numChannels
                        
                        sourceString = ...
                            extractor.listing(extractor.currentStack).name;
                        sourceString = strrep(sourceString,...
                            extractor.waveLengthsStrings{1},...
                            extractor.waveLengthsStrings{cc});
                        extraReaders{cc-1} = bfGetReader(sourceString);
                        
                    end
                    
                end
                
                nuclei = cell(1,numImages);
                
                central_z_ind = round(stackDim(3)./2);

                for bb = 1:numImages
                    
%                     disp(bb./numImages)
                    
                    % --- read in the stack
                    
                    rawStack_cell = cell(1,numChannels);
                    
                    reader_1.setSeries(bb-1);
                    
                    if extractor.sepWavelengthsFlag
                        for ee = 1:numChannels-1
                            extraReaders{ee}.setSeries(bb-1);
                        end
                    end
                    
                    for cc = 1:numChannels
                        
                        rawStack_cell{cc} = ...
                            zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ);
                        
                        for zz = 1:rawStackSizeZ
                            
                            if cc == 1 || ~extractor.sepWavelengthsFlag
                                planeInd = ...
                                    reader_1.getIndex(zz-1,cc-1,0)+1;
                                rawStack_cell{cc}(:,:,zz) = ...
                                    bfGetPlane(reader_1,planeInd);
                            else
                                planeInd = ...
                                    reader_1.getIndex(zz-1,0,0)+1;
                                rawStack_cell{cc}(:,:,zz) = ...
                                    bfGetPlane(extraReaders{cc-1},planeInd);
                            end
                        end
                        
                    end
                    
                    % --- segmentation

                    central_DNA_image = squeeze(...
                        rawStack_cell{extractor.segChannel}(:,:,central_z_ind));

                    maxInt = max(central_DNA_image(:));
                    binEdges = unique(round(linspace(0,maxInt,500)));
                    
                    counts = histc(central_DNA_image(:),binEdges);
                    
                    counts = counts(1:end-1);
                    binEdges = binEdges(1:end-1);
                    
                    binEdges = (double(binEdges(counts>0))).';
                    counts = double(counts(counts>0));
                    
                    counts = counts./max(counts);
                    
                    [threshVal,~] = otsuLimit(binEdges,counts,[0,Inf]);
                    
                    segStack = rawStack_cell{extractor.segChannel}>threshVal;
                    
                    connectedComponents = bwconncomp(segStack,26);
                    regionProperties = regionprops(connectedComponents,...
                        'Area','FilledArea','BoundingBox');
                    
                    volumes = [regionProperties.Area];
                    filledVolumes = [regionProperties.Area];
                    
                    regionProperties = regionProperties(...
                        volumes>(extractor. minVol./prod(voxelSizes)) ...
                        & (filledVolumes./volumes)>extractor.minFilledVol);
                    
                    numRegions = numel(regionProperties);
                    
                    nuclei{bb} = cell(1,numRegions);
                    
                    for rr = 1:numRegions
                        
                        nuclei{bb}{rr} = cell(1,numChannels);
                        
                        extractLimits = regionProperties(rr).BoundingBox;
                        lowLims = extractLimits(1:3)+0.5;
                        upLims = lowLims + extractLimits(4:6)-1;
                        
                        for cc = 1:numChannels
                            
                            nuclei{bb}{rr}{cc} = ...
                                rawStack_cell{cc}(lowLims(2):upLims(2),...
                                lowLims(1):upLims(1),lowLims(3):upLims(3));
                            
                        end
                        
                    end
                    
                end
                
                reader_1.close()
                   
                nuclei = [nuclei{:}];
                
                % Increment to next stack
                extractor.currentStack = ...
                    extractor.currentStack + 1;
                
            else
                
                nuclei = [];
                voxelSizes = [];
                
            end
            extractor.hasNextStack = ...
                    extractor.currentStack<=extractor.numStacks;
        end
        
    end
    
end