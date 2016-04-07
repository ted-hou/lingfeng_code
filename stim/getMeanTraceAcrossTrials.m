function nMean = getMeanTraceAcrossTrials(nReshaped)
nStimLocs = length(nReshaped);
nTraces = size(nReshaped(1).nPreStim, 1);

for iStimLoc = 1:nStimLocs
	nPreStim 	= nReshaped(iStimLoc).nPreStim;

	nMean(iStimLoc).nPreStim 	= squeeze(mean(mean(nReshaped(iStimLoc).nPreStim(:, :, :, 2:end), 4), 2));
	nMean(iStimLoc).nPostStim 	= squeeze(mean(mean(nReshaped(iStimLoc).nPostStim(:, :, :, 1:end-1), 4), 2));
	nMean(iStimLoc).nStim 		= mean(mean(nReshaped(iStimLoc).nStim, 3), 2);
	nMean(iStimLoc).nAll 		= horzcat(nMean(iStimLoc).nPreStim, nMean(iStimLoc).nStim, nMean(iStimLoc).nPostStim);

	arr = nReshaped(iStimLoc).nPreStim; 	nMean(iStimLoc).stdPreStim 	= squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshaped(iStimLoc).nPostStim; 	nMean(iStimLoc).stdPostStim = squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshaped(iStimLoc).nStim; 		nMean(iStimLoc).stdStim		= std(reshape(arr, [nTraces, size(arr, 2)*size(arr, 3)]), 0, 2);
	nMean(iStimLoc).stdAll = horzcat(nMean(iStimLoc).stdPreStim, nMean(iStimLoc).stdStim, nMean(iStimLoc).stdPostStim);
end