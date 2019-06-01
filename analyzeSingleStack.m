function [nuc_int,cyto_int,nuc_area,nuc_solidity,nuc_eccentricity,...
    CoV_vec,Covar_matrix,nucImage] = ...
    analyzeSingleStack(rawStack,segChannel,pixelArea,minArea,maxArea)

plotSegmentation = false;

% --- find center slice by maximum contrast

[rawStackSizeY,rawStackSizeX,rawStackSizeZ] = ...
    size(rawStack{segChannel});

numChannels = numel(rawStack);

contrastVec = zeros(1,rawStackSizeZ);

for ll = 1:rawStackSizeZ
    
    % Use maximal intensity to pick center section
    section = squeeze(rawStack{segChannel}(:,:,ll));
    contrastVec(ll) = mean(section(:));
    
end

[~,maxContrastInd] = max(contrastVec);

maxContrastImage = ...
    cellfun(@(elmt)squeeze(elmt(:,:,maxContrastInd)),...
    rawStack,'UniformOutput',false);

% --- segmentation

medianFiltImage = ...
    medfilt2(maxContrastImage{segChannel},[5 5]);

maxInt = max(medianFiltImage(:));
binEdges = unique(round(linspace(0,maxInt,500)));

counts = histc(medianFiltImage(:),binEdges);

counts = counts(1:end-1);
binEdges = binEdges(1:end-1);

binEdges = (double(binEdges(counts>0))).';
counts = double(counts(counts>0));

counts = counts./max(counts);

[threshVal,otsuMetric] = otsuLimit(binEdges,counts,[0,Inf]);

segImage = medianFiltImage>=1.2.*threshVal;

% --- erode/dilate treatment with hole-filling

dilateCycles = 6;
dilateMeasure = 6;
erodeCycles = 12;

% Make neighborhood mask for operations
neighborhood = ...
    [0,1,0;...
    1,1,1;...
    0,1,0];

for kk = 1:dilateCycles
    segImage = imdilate(segImage,neighborhood);
end

% Fill holes in segmented image
segImage = imfill(segImage,'holes');

for kk = 1:erodeCycles
    segImage = imerode(segImage,neighborhood);
end

% Fill holes in segmented image
segImage = imfill(segImage,'holes');

%--- erode/dilate treatment over



if plotSegmentation
    
    figure(1)
    
    subplot(3,2,1)
    
    plot(binEdges,counts,'k-')
    hold on
    plot([1,1].*threshVal,[0,1],'k--')
    hold off
    
    xlabel('Intensity')
    ylabel('Frequency')
    
    subplot(3,2,3)
    
    plot(1:rawStackSizeZ,contrastVec,'k-')
    
    hold on
    
    hold off
    
    xlabel('z section')
    ylabel('Contrast')
    
    
    
    subplot(3,2,2)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel})
    
    axis equal
    
    title('Segmentation')
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    colormap(gray)
    
    
    
    subplot(3,2,4)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel})
    
    axis equal
    
    title('Quantification')
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    colormap(gray)
    
    
    
    
    subplot(3,2,6)
    
    imagesc([0,rawStackSizeX],[0,rawStackSizeY],...
        maxContrastImage{segChannel}.*segImage)
    
    axis equal
    
    set(gca,'XLim',[0,rawStackSizeX],'YLim',[0,rawStackSizeY],...
        'YDir','normal')
    
    title(sprintf('%3.3f um',sum(segImage(:))))
    
    colormap(jet)
    
end


% --- Segmentation and quantification

segIntImg = maxContrastImage;
segBinImg = segImage;

% --- segment out nucleus closest to center of slice

regions = bwconncomp(segBinImg);
regionProps = regionprops(regions,'Area');
regionAreas = [regionProps(:).Area];


% Restrict to objects in area range
validVolInds = regionAreas.*pixelArea>=minArea ...
    & regionAreas.*pixelArea<=maxArea;
nucleiRegions = struct;
nucleiRegions.Connectivity = regions.Connectivity;
nucleiRegions.ImageSize = regions.ImageSize;
nucleiRegions.NumObjects = sum(validVolInds);
nucleiRegions.PixelIdxList = regions.PixelIdxList(validVolInds);

% Extract further properties only for the valid regions

nucleiProps = regionprops(nucleiRegions,segIntImg{segChannel},...
    'Area','WeightedCentroid','MeanIntensity','BoundingBox','Image',...
    'Solidity','Eccentricity');

% Get object closest to center of the image
nucleiCentroid = {nucleiProps(:).WeightedCentroid};
[~,minInd] = min(cellfun(@(elmt)...
    sum((elmt-[rawStackSizeX./2,rawStackSizeY./2]).^2),nucleiCentroid));

if numel(nucleiCentroid)<1
    
    nuc_area = NaN;
    nuc_solidity = NaN;
    nuc_eccentricity = NaN;
    nuc_int = zeros(1,numChannels);
    nuc_int(:) = NaN;
    cyto_int = zeros(1,numChannels);
    cyto_int(:) = NaN;
    CoV_vec = zeros(1,numChannels);
    CoV_vec(:) = NaN;
    Covar_matrix = zeros(numChannels,numChannels);
    Covar_matrix(:) = NaN;
    
    nucImage = {};
    
    return;
    
end

maxArea = nucleiProps(minInd).Area.*pixelArea;
maxSolidity = nucleiProps(minInd).Solidity;
maxEccentricity = nucleiProps(minInd).Eccentricity;

maxBBox = nucleiProps(minInd).BoundingBox;
maxImage = nucleiProps(minInd).Image;

nuc_area = maxArea;
nuc_solidity = maxSolidity;
nuc_eccentricity = maxEccentricity;

segInds = nucleiRegions.PixelIdxList(minInd);
segInds = segInds{1};

segMask = false(size(segBinImg));
segMask(segInds) = true;

% --- quantification of nuclear and cytoplasmic intensity

nucBBox = maxBBox;

% Make mask to determine cytoplasmic intensity

totalDil = dilateCycles+erodeCycles;

% Determine small bounding box (without dilate)
smallMinY = nucBBox(2)+0.5;
smallMinX = nucBBox(1)+0.5;
smallMaxY = nucBBox(2)+nucBBox(4)-0.5;
smallMaxX = nucBBox(1)+nucBBox(3)-0.5;


% Determine extended bounding box (after dilate)
fullExtMinY = smallMinY-totalDil;
fullExtMinX = smallMinX-totalDil;
fullExtMaxY = smallMaxY+totalDil;
fullExtMaxX = smallMaxX+totalDil;

% Limit extended bounding box to within image limits
extMinY = max(1,fullExtMinY);
yLoDiff = extMinY - fullExtMinY;
extMinX = max(1,fullExtMinX);
xLoDiff = extMinX - fullExtMinX;
extMaxY = min(rawStackSizeY,fullExtMaxY);
yHiDiff = fullExtMaxY - extMaxY;
extMaxX = min(rawStackSizeX,fullExtMaxX);
xHiDiff = fullExtMaxX - extMaxX;

% Extended bounding box size
extSizeY = extMaxY - extMinY + 1;
extSizeX = extMaxX - extMinX + 1;

% Inclusion mask
inclMask = zeros(extSizeY,extSizeX);
inclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
    (1+totalDil-xLoDiff):(end-totalDil+xHiDiff))...
    = maxImage;

% Exclusion mask
exclMask = segMask(extMinY:extMaxY,extMinX:extMaxX);

dilateNeighborhood = zeros(3,3);
dilateNeighborhood(2,2) = 1;
dilateNeighborhood(1,2) = 1;
dilateNeighborhood(3,2) = 1;
dilateNeighborhood(2,1) = 1;
dilateNeighborhood(2,3) = 1;

for dd = 1:dilateCycles
    inclMask = imdilate(inclMask,dilateNeighborhood);
    exclMask = imdilate(exclMask,dilateNeighborhood);
end

for dd = 1:dilateMeasure
    inclMask = imdilate(inclMask,dilateNeighborhood);
end

measureMask = inclMask & ~exclMask;

cyto_int = zeros(1,numChannels);
nuc_int = zeros(1,numChannels);
CoV_vec = zeros(1,numChannels);
Covar_matrix = zeros(numChannels,numChannels);
nucImage = cell(1,numChannels);
nucInts = cell(1,numChannels);

for cc = 1:numChannels
    nucImage{cc} = maxContrastImage{cc}(extMinY:extMaxY,extMinX:extMaxX);
    nucInts{cc} = maxContrastImage{cc}(segMask);
    nucInts{cc} = nucInts{cc}(:);
    nuc_int(cc) = mean(nucInts{cc});
    cyto_int(cc) = mean(nucImage{cc}(measureMask));
    CoV_vec(cc) = std(nucInts{cc})./nuc_int(cc);
end

for cc = 1:numChannels
    
    for dd = cc:numChannels
        Covar_matrix(cc,dd) = mean(...
            (nucInts{cc}./nuc_int(cc)-1).*(nucInts{dd}./nuc_int(dd)-1))...
            ./sqrt(mean((nucInts{cc}./nuc_int(cc)-1).^2))...
            ./sqrt(mean((nucInts{dd}./nuc_int(dd)-1).^2));
        Covar_matrix(dd,cc) = Covar_matrix(cc,dd);
    end
        
end

end