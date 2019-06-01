This repository provides the nucleiExtractor object, which can be used to extract substacks from 3D microscopy data. These substacks each contain one nucleus. This analysis code is best used on micorscopy data where the nuclie have been labelled with a DNA stain. It handles only single time points, but multiple positions and also multiple color channels are fine.

The repository furhter contains an example script (MainAnalysisFile.m) and all additional script/code required to execute the analysis of an example data set from our lab. The data are contained in the ExampleData folder.

As the code is using the bfmatlab importer from BioFormats, it should work out of the box with most manufacturer's microscopes.

Credits go to
> The Open Microscopy Environment (OME) for supplying BioFormats, in this case for MatLab via the bfmatlab tool
> parfor, distributionPlot from the MatLab File Exchange
