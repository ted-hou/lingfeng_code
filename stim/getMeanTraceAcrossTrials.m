function nMean = getMeanTraceAcrossTrials(nReshaped)
nStimLocs = length(nReshaped);
nTraces = size(nReshaped(1).preStim, 1);

for iStimLoc = 1:nStimLocs
	nPreStim 	= nReshaped(iStimLoc).preStim;

	nMean(iStimLoc).preStim 	= squeeze(mean(mean(nReshaped(iStimLoc).preStim(:, :, :, 2:end), 4), 2));
	nMean(iStimLoc).postStim 	= squeeze(mean(mean(nReshaped(iStimLoc).postStim(:, :, :, 1:end-1), 4), 2));
	nMean(iStimLoc).stim 		= mean(mean(nReshaped(iStimLoc).stim, 3), 2);
	nMean(iStimLoc).all 		= horzcat(nMean(iStimLoc).preStim, nMean(iStimLoc).stim, nMean(iStimLoc).postStim);

	arr = nReshaped(iStimLoc).preStim; 	nMean(iStimLoc).stdPreStim 	= squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshaped(iStimLoc).postStim; nMean(iStimLoc).stdPostStim = squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshaped(iStimLoc).stim; 	nMean(iStimLoc).stdStim		= std(reshape(arr, [nTraces, size(arr, 2)*size(arr, 3)]), 0, 2);
	nMean(iStimLoc).stdAll = horzcat(nMean(iStimLoc).stdPreStim, nMean(iStimLoc).stdStim, nMean(iStimLoc).stdPostStim);
end