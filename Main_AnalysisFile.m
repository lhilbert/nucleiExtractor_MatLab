clear all

addpath('./bfmatlab') % make sure bfmatlab functions are available

% Should the code use parallel processing for speed-up?
parfor_flag = false; % Change this value to true to activate parallel
ps = parallel.Settings;
ps.Pool.AutoCreate = false;

% Where are your data located?
mainDirectory = ...
    './ExampleData';

% How to find all files for a given condition?
% The search mechanism is recursive, and has a directory/path wild card,
% which is **, and also a file name wild card, which is *
% So, for example '**/*CTRL*.tif' will find all files within all
% subdirectories that have CTRL anywhere in their name. No need to keep
% very close track of where which files are stored.
conditionDirs = {...
    '/Ctrl_1/*.tif',...
	'/Ctrl_2/*.tif',...
	'/ActA_1/*.tif',...
	'/ActA_2/*.tif'};


% What labels should be used to name the different conditions?
condLabels = {...
    'Control (1)','Control (2)',...
	'Activin A (1)','Activin A (2)',...
	};

% flag for separate wavelengths
%useStrings = {'w0000','w0001','w0002','w0003'};
wavelengthsStrings = {[],[],[],[]};
% This option is retained for the odd cases where acquisition software
% decides to save image data into separate files for different wavelengths.
% If that is not the case, simply use an empty square bracket array, as
% done above

% Maximal number of image stacks to use per condition, Inf to use all
% This option is useful for large datasets, when you first need to analyze
% just a subpart fo the dataset to get the analysis procedure to work
maximumStackCap = 20;

% Choose which channel to use for segmentation, typically DNA stain
seg_channel = 1; % Find out channel number by opening example data in FIJI

% Inclusion criteria for the 3D segmentation step
minVol = 4.^2; % in cubic micrometers
minFilledVol = 0.2; % Fraction of the bounding box that must be filled

%% Extract structural features

% Which of the conditions to actually use?
useConds = [1:4];

% Do not change from here ...
conditionDirs = conditionDirs(useConds);
condLabels = condLabels(useConds);
wavelengthsStrings = wavelengthsStrings(useConds);
numConds = numel(conditionDirs);
% ... to here.

% --- Preparation of arrays to store the information obtained from the
% individual nuclei, sorted by the different conditions

nucStack_cell = cell(1,numConds); % Substacks containing nuclei
voxelDim_cell = cell(1,numConds); % Physical dimensions of xyz voxels
nucImage_cell = cell(1,numConds); % Mid-section images of the nuclei

% Storage of geometric features
area_cell = cell(1,numConds);
solidity_cell = cell(1,numConds);
eccentricity_cell = cell(1,numConds);

% Storage of color intensity, coefficient of variation, and covariance
int_cell = cell(1,numConds);
cov_cell = cell(1,numConds);
covar_cell = cell(1,numConds);

for nn = 1:numConds

    fprintf('Processing condition %d of %d\n\n',nn,numConds)

	% Construct the search string for this condition
    thisSearchString = fullfile(mainDirectory,conditionDirs{nn});
    
    % Now the "magic" happens, we will instantiate an extractor object,
    % which will subsequently serve up the data associated with the nuclei
    % detected in all the source data fitting the search string	
	extractor = nucleiExtractor( ...
		thisSearchString,seg_channel,minVol,minFilledVol,...
		wavelengthsStrings{nn});
    % Note: you can also leave out the last argument, in case all
    % wavelengths are provided within the same file. The function call
    % works in both cases.
	
	% The extractor object provides the number of substacks it has
	% obtained. Each substack contains a nucleus (assuming that the raw
	% data and the subsequent analysis are appropriate for the
	% sub-segmentation of nuclei, that is).
    numStacks = extractor.numStacks;
    
	% Here, we apply a cut-off in terms of stacks to work with, so we
	% can run fast analyses on only a part of the data. This is useful for
	% troubleshooting the image analysis without always processing the
	% entire dataset.
    numStacks = min([numStacks,maximumStackCap]);
    
    % Allocation of storage variables
	voxelDim_cell{nn} = cell(1,numStacks);
	nucImage_cell{nn} = cell(1,numStacks);
	
	area_cell{nn} = cell(1,numStacks);
    solidity_cell{nn} = cell(1,numStacks);
    eccentricity_cell{nn} = cell(1,numStacks);

    int_cell{nn} = cell(1,numStacks);
    cov_cell{nn} = cell(1,numStacks);
	covar_cell{nn} = cell(1,numStacks);
	
	% From here ...
	if ~parfor_flag
		parfor_progress(numStacks);
	end
	% ... to here only to give a progress bar. Never mind that it's called
	% parfor_progress, it is not used in parallel fashion in this code!
	
    for mm = 1:numStacks
        
		% Get all nuclei from one of the image stacks inside the extractor
		% object. currentNuclei is a cell array containing the substacks of
		% the individual nuclei inside of this particular stack.
        [currentNuclei,voxelSizes] = extractor.retrieveNuclei;
        
		voxelDim_cell{nn}{mm} = voxelSizes;

		pixelArea = voxelSizes(1).*voxelSizes(2);
                
        numNuclei = numel(currentNuclei);        
        nucImage_cell{nn}{mm} = cell(1,numNuclei);
        
        int_cell{nn}{mm} = cell(1,numNuclei);
        cov_cell{nn}{mm} = cell(1,numNuclei);
        covar_cell{nn}{mm} = cell(1,numNuclei);
        
        area_cell{nn}{mm} = zeros(1,numNuclei);
        solidity_cell{nn}{mm} = zeros(1,numNuclei);
        eccentricity_cell{nn}{mm} = zeros(1,numNuclei);
		
        for ll = 1:numNuclei
			
			% Quantify different properties form the substack of a given
			% nucleus
            [nuc_int,cyto_int,nuc_area,nuc_solidity,nuc_eccentricity,...
                CoV_vec,Covar_matrix,nucImage] = ...
                analyzeSingleStack(currentNuclei{ll},...
                seg_channel,pixelArea,0,Inf);

			if ~isnan(nuc_area)
				
				nucImage_cell{nn}{mm}{ll} = nucImage;
				
				% Here, we use for the intensity the measure of
				% nuclear-over-cytoplasmic ratio. Different ones are
				% possible, such as the difference or the direct value.
				int_cell{nn}{mm}{ll} = nuc_int./cyto_int;
				
				cov_cell{nn}{mm}{ll} = CoV_vec;
				covar_cell{nn}{mm}{ll} = Covar_matrix;
				
				area_cell{nn}{mm}(ll) = nuc_area;
				solidity_cell{nn}{mm}(ll) = nuc_solidity;
				eccentricity_cell{nn}{mm}(ll) = nuc_eccentricity;
				
            end
            
		end
		
		% Do not change from here ...
		if ~parfor_flag
			parfor_progress;
		end
		% ... to here.
		
	end
	
	% Do not change from here ...
	if ~parfor_flag
		parfor_progress(0);
	end
    % ... to here.
	
	% Now we will pool the data from the different stacks into bigger,
	% joint cell arrays. up until this point, they are still cell arrays
	% within a cell array, as they came from the separate main stacks.
	
	nucImage_cell{nn} = [nucImage_cell{nn}{:}];
    
	int_cell{nn} = [int_cell{nn}{:}];
	cov_cell{nn} = [cov_cell{nn}{:}];
	covar_cell{nn} = [covar_cell{nn}{:}];

	area_cell{nn} = [area_cell{nn}{:}];
	solidity_cell{nn} = [solidity_cell{nn}{:}];
	eccentricity_cell{nn} = [eccentricity_cell{nn}{:}];
    
	% Now we still have to kick out those cases where no actual nucleus
	% could be successfully segmented. We do this by checking the first
	% channels intensity value, whether it is a real/integer number.
    
	validInds = cellfun(@(elmt)~isempty(elmt),int_cell{nn});
    
	nucImage_cell{nn} = nucImage_cell{nn}(validInds);
	
	int_cell{nn} = int_cell{nn}(validInds);
	cov_cell{nn} = cov_cell{nn}(validInds);
	covar_cell{nn} = covar_cell{nn}(validInds);
	
	area_cell{nn} = area_cell{nn}(validInds);
	solidity_cell{nn} = solidity_cell{nn}(validInds);
	eccentricity_cell{nn} = eccentricity_cell{nn}(validInds);
	
end


%% Distributions of population parameters

% --- Normalization
% Indicate for each condition which other condition should be used to
% normalize. Put 0 for no normalization.
normRefs = [0,0,0,0];

numChannels = numel(int_cell{1}{1}); % Get the number of color channels
% Only works if the number of color channels is the same across the whole
% dataset

channel_ints = cell(1,numChannels);
channel_covs = cell(1,numChannels);
crossChannel_covars = cell(numChannels,numChannels);

for channel = 1:numChannels
	
	channel_ints{channel} = cellfun(...
		@(cond_elmt)cellfun(@(elmt)elmt(channel),cond_elmt),...
		int_cell,'UniformOutput',false);

	channel_scaling_vec = ...
		cellfun(@(xx)prctile(xx,80),channel_ints{channel});
	
	for cond = 1:numConds
		if normRefs(cond)>0
			channel_ints{channel}{cond} = channel_ints{channel}{cond} ...
				./channel_scaling_vec(normRefs(cond));
		end
	end

	channel_covs{channel} = cellfun(...
		@(cond_elmt)cellfun(@(elmt)elmt(channel),cond_elmt),...
		cov_cell,'UniformOutput',false);

	for crossChannel = 1:numChannels
		crossChannel_covars{channel,crossChannel} = cellfun(...
			@(cond_elmt)cellfun(...
			@(elmt)elmt(channel,crossChannel),cond_elmt),...
			covar_cell,'UniformOutput',false);
	end
	
end

% --- Selection criteria to remove not nucleus-shaped objects

minArea = 10;
minRNA = -Inf;
maxRNA = +Inf;
minSolidity = 0.7;
maxEccentricity = 0.9;

markerSize = 2;

include_inds = cell(1,numConds);

for cc = 1:numConds
    
    include_inds{cc} = ...
        area_cell{cc}>=minArea&...
        solidity_cell{cc}>=minSolidity&...
        eccentricity_cell{cc}<=maxEccentricity;    
end

DNA_vals_gated = cellfun(@(xx,inds)xx(inds),...
	channel_ints{1},include_inds,...
    'UniformOutput',false);
GFP_vals_gated = cellfun(@(xx,inds)xx(inds),...
	channel_ints{2},include_inds,...
    'UniformOutput',false);

DNA_cov_gated = cellfun(@(xx,inds)xx(inds),...
	channel_covs{1},include_inds,...
    'UniformOutput',false);
GFP_cov_gated = cellfun(@(xx,inds)xx(inds),...
	channel_covs{2},include_inds,...
    'UniformOutput',false);

DNA_GFP_covar_gated = cellfun(@(xx,inds)xx(inds),...
	crossChannel_covars{1,2},include_inds,...
    'UniformOutput',false);


% Pick which features should be shown as distributions
valCell_pool = {DNA_vals_gated,GFP_vals_gated,...
	DNA_cov_gated,GFP_cov_gated,DNA_GFP_covar_gated};
ylabel_pool = {'I_{DNA} gated','I_{GFP} gated',...
	'C_{DNA}','C_{GFP}','CoVar_{DNA-GFP}'};

numVals = numel(valCell_pool);

figure(1)
clf

for ll = 1:numVals
    
    subplot(1,numVals,ll)
    
    valCell = valCell_pool{ll};
    
    hold on
    
    plotData = [valCell{:}]';
    groupingVec = [];
    for cc = 1:numConds
        
        groupingVec = [groupingVec;cc.*ones(numel(valCell{cc}),1)];
        
    end
    
    plotHandles = distributionPlot(plotData,'groups',groupingVec,...
        'addSpread',1,'distWidth',0.0,'showMM',0,'histOpt',0,'divFactor',8);
    
    hold on
    
    numHandles = numel(unique(groupingVec));
    
    for kk = 1:numHandles
                    
        set(plotHandles{4}{1}(kk),'MarkerSize',markerSize,...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor',[0.5,0.5,0.5],...
            'Marker','o')
        
    end
    
    plotHandles = zeros(1,2);    
    
    meanVector = cellfun(@mean,valCell);
    SEMVector = cellfun(@std,valCell)...
        ./sqrt(cellfun(@(xx)numel(xx),DNA_vals_gated));
    
    for kk = 1:numConds
        plot(kk.*[1,1],meanVector(kk)+[SEMVector(kk),-SEMVector(kk)],'k-',...
            'LineWidth',1)
        plot(kk+[-0.2,+0.2],meanVector(kk).*[1,1],'k-',...
            'LineWidth',1)
        plot(kk+[-0.1,+0.1],meanVector(kk).*[1,1]+SEMVector(kk),'k-',...
            'LineWidth',1)
        plot(kk+[-0.1,+0.1],meanVector(kk).*[1,1]-SEMVector(kk),'k-',...
            'LineWidth',1)
    end
    
    set(gca,'XLim',[0.5,numConds+0.5],'Box','on')
    set(gca,'XTick',1:numConds,'XTickLabel',condLabels)
    
    ylabel(ylabel_pool{ll})
    
end