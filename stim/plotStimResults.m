function [n, P] = plotStimResults(dF, spikeThreshold, n, P, stimParams, expName, expId, forceStackHeight)
% If there are more than 4 stimulation positions, might need more colors. But 4 should be enough.
if isempty(expName)
	expName = '';
end

if isempty(expId)
	expId = 0;
end

colors = 'rgbc';

scopeStimParams 	= getScopeStimParams(stimParams);
scopeStimArtefact 	= stimParams.scopeStimArtefact;
stimLog 			= stimParams.stimLog;
stimsPerTrain		= stimParams.stimsPerTrain;
mirrorPosList		= stimParams.mirrorPosList;
scopeFramePeriod	= 1/stimParams.scanFrameRate;
nTraces 			= size1(dF);
nFramesPerTrace		= size2(dF);
nStimLocs 			= length(mirrorPosList);
%% 1. Remove stim artefact, interpolate for stim frames
dF = removeStimArtefact(dF, scopeStimArtefact);

dFStackHeight = 0.25*maxel(dF);

dFStacked = dF;
dFStacked = horzcat(dFStacked, max(dF, [], 2));
dFStacked = sortrows(dFStacked, size(dFStacked, 2));
dFStacked = dFStacked(:,1:(size(dFStacked,2) - 1));

for iTrace = 1:nTraces
	dFStacked(iTrace, :) = dFStacked(iTrace, :) + (iTrace - 1)*dFStackHeight;
end


% time in seconds
t = scopeFramePeriod * (1:nFramesPerTrace);

% Plot all traces + stim times
% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% hold on;
% for iMirrorPos = 1:nStimLocs;
% 	stimTrace = zeros(1, nFramesPerTrace);
% 	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
% 	stimTrace(frames') = (1 + nTraces)*dFStackHeight;
% 	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
% end

% plot(t, dFStacked', 'k'), hold off;

% ylim([0, (1 + nTraces)*dFStackHeight])
% xlim([0, maxel(t)])
% xlabel('Time (s)')
% ylabel('Neurons/GCaMP Signal (\DeltaF/F)')
% title('\DeltaF/F')
% ax = gca;
% ax.YTick = dFStackHeight*(0.5:5:(nTraces - 0.5));	
% ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
% ax.TickLength = [0 0];


%% 2. Deconvolve the traces using oopsi (fnnd) Volgostein (2010)
if isempty(n) || isempty (P) % Skipped if deconvolution already done since this is quite slow.
	for iTrace = 1:nTraces
		[n(iTrace,1:nFramesPerTrace), P{iTrace}] = fast_oopsi(dF(iTrace,:));
	end
end

% Produce a stacked version of spike probabilities for plotting
nStackHeight = 1;

nStacked = n;
nStacked = horzcat(nStacked, max(dF, [], 2));
nStacked = sortrows(nStacked, size(nStacked, 2));
nStacked = nStacked(:,1:(end - 1));

for iTrace = 1:nTraces
	nStacked(iTrace, :) = nStacked(iTrace, :) + (iTrace - 1)*nStackHeight;
end

% Plot all spike probabilities + stim times
% figure('units', 'normalized', 'outerposition', [0 0 1 1])
% hold on;
% for iMirrorPos = 1:nStimLocs;
% 	stimTrace = zeros(1, nFramesPerTrace);
% 	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
% 	stimTrace(frames') = (1 + nTraces)*nStackHeight;
% 	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
% end
% plot(t, nStacked', 'k'), hold off;

% ylim([0, (1 + nTraces)*nStackHeight])
% xlim([0, maxel(t)])
% xlabel('Time (s)')
% ylabel('Neurons/Spike Probability')
% title('Deconvoluted Spike Probability')
% ax = gca;
% ax.YTick = nStackHeight*(0.5:5:(nTraces - 0.5));	
% ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
% ax.TickLength = [0 0];

% Plot all spikes + stim times
if isempty(spikeThreshold)
	spikeThreshold = 0.7;
end
nThresholded = zeros(size(n));
nThresholded(n>spikeThreshold) = 1;
nStackedAndThresholded = horzcat(nThresholded, max(dF, [], 2));
nStackedAndThresholded = sortrows(nStackedAndThresholded, size(nStackedAndThresholded, 2));
nStackedAndThresholded = nStackedAndThresholded(:,1:(end - 1));

for iTrace = 1:nTraces
	nStackedAndThresholded(iTrace, :) = nStackedAndThresholded(iTrace, :) + (iTrace - 1)*nStackHeight;
end

% figure('units', 'normalized', 'outerposition', [0 0 1 1])
% hold on;
% for iMirrorPos = 1:nStimLocs;
% 	stimTrace = zeros(1, nFramesPerTrace);
% 	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
% 	stimTrace(frames') = (1 + nTraces)*nStackHeight;
% 	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
% end
% plot(t, nStackedAndThresholded', 'k'), hold off;

% ylim([0, (1 + nTraces)*nStackHeight])
% xlim([0, maxel(t)])
% xlabel('Time (s)')
% ylabel('Neurons/Spikes')
% title(['Deconvoluted Spikes (Threshold = ', num2str(spikeThreshold), ')'])
% ax = gca;
% ax.YTick = nStackHeight*(0.5:5:(nTraces - 0.5));	
% ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
% ax.TickLength = [0 0];

%% 4. Plot spike timing (relative to stimulation offset) using the deconvoluted spike probability traces
% i) 	Extract and reshape traces (split into bins/trials/...)
% ii) 	Get mean binned traces for plotting
% iii) 	Run some statistical tests
% iv)	Now plot dat ting

% First sort the raw traces according to max spike amplitude
nSorted = n;
% nSorted = horzcat(nSorted, max(dF, [], 2));
% nSorted = sortrows(nSorted, size(nSorted, 2));
% nSorted = nSorted(:,1:(end - 1));

% i) Reshape the deconvolute spike probability traces n(iTrace, iFrame) into 
% 		nReshaped(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex) and 
% 		nReshapedBaseline(iTrace, iFrame, iBin, iTrial)
% [nReshaped, nReshapedBaseline] = reshapeTraces(nSorted, nBins, stimParams);

% framesPerBin = size2(nReshapedBaseline);

% % ii) Get binned mean traces for plotting
% for iStimLoc = 1:nStimLocs
% 	nMean(iStimLoc).n = squeeze(mean(mean(nReshaped(iStimLoc).n, 4), 2));
% end
% nMeanBaseline = squeeze(mean(mean(nReshapedBaseline, 4), 2));

% % iii) Perform two-sample t test between stim and baseline traces, assuming unknown and unequal variances.
% for iStimLoc = 1:nStimLocs
% 	for iTrace = 1:nTraces
% 		for iBin = 1:nBins
% 			xStim = nReshaped(iStimLoc).n(iTrace, :, iBin, :);
% 			xBase = nReshapedBaseline(iTrace, :, iBin, :);
% 			h(iStimLoc).h(iTrace, iBin) = ttest2(xStim(:), xBase(:), 'Alpha', 0.01, 'Tail', 'right', 'Vartype', 'unequal');
% 		end
% 	end
% end

% % iv) Now plot dat ting
% nMeanStackHeight = maxel([nMean.n]);
% nMeanStacked = nMean;
% nMeanStackedBaseline = nMeanBaseline;

% t = scopeFramePeriod*(0.5*framesPerBin:framesPerBin:(nBins - 0.5)*framesPerBin);

% for iStimLoc = 1:nStimLocs
% 	nMeanStacked(iStimLoc).n = nMeanStacked(iStimLoc).n + repmat((0:(nTraces - 1))', [1, nBins]).*nMeanStackHeight;
% 	nMeanStackedBaseline = nMeanStackedBaseline + repmat((0:(nTraces - 1))', [1, nBins]).*nMeanStackHeight;
% end

% figWidth = 0.1;

% for iStimLoc = 1:nStimLocs
% 	figure('units', 'normalized', 'outerposition', [min(figWidth*(iStimLoc - 1), 1 - figWidth), 0, figWidth, 1]), hold on
% 	for iTrace = 1:nTraces
% 		xmarkers = t(find(h(iStimLoc).h(iTrace, :)==1));
% 		ymarkers = nMeanStacked(iStimLoc).n(iTrace, find(h(iStimLoc).h(iTrace, :)==1));
% 		plot(t, nMeanStacked(iStimLoc).n(iTrace, :), 'k', xmarkers, ymarkers, 'k*', 'MarkerSize', 3);
% 		plot(t, nMeanStackedBaseline(iTrace, :), 'k--');
% 		plot(t, nMeanStackHeight*(iTrace - 1) + zeros(size(t)), 'k:');
% 	end
% 	plot(t, nMeanStackHeight*iTrace + zeros(size(t)), 'k:');
% 	hold off
% 	xlabel('Time After Stimulation (s)')
% 	ylabel('Neurons/Mean Spike Probability Across Trials')
% 	title(['Mean Traces Across Trials, Stim Loc ', num2str(iStimLoc)])
% 	ylim([0, (1 + nTraces)*nMeanStackHeight])
% 	xlim([0, maxel(t)])
% 	ax = gca;
% 	ax.YTick = nMeanStackHeight*(0.5:5:(nTraces - 0.5));	
% 	ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
% 	ax.TickLength = [0 0];
% end


%% 5. Plot aligned stim traces
% i) 	Extract deconvoluted traces from t-5 till t+5, (t=stim time), reshaped into nReshapedStim(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex)
% ii) 	Calculate the mean

% i)	Extract traces 5 seconds before and after stimulation
timeRadius = 15;
secondsPerBin = 0.5;
% secondsPerBin = 0.1;
[nReshapedStim] = reshapeStimTraces(nSorted, timeRadius, secondsPerBin, stimParams);

% ii)	Calculate mean traces for each ROI
for iStimLoc = 1:nStimLocs
	nPreStim 	= nReshapedStim(iStimLoc).nPreStim;

	% nMeanStim(iStimLoc).nPreStim 	= squeeze(mean(mean(nReshapedStim(iStimLoc).nPreStim, 4), 2));
	% nMeanStim(iStimLoc).nPostStim 	= squeeze(mean(mean(nReshapedStim(iStimLoc).nPostStim, 4), 2));
	nMeanStim(iStimLoc).nPreStim 	= squeeze(mean(mean(nReshapedStim(iStimLoc).nPreStim(:, :, :, 2:end), 4), 2));
	nMeanStim(iStimLoc).nPostStim 	= squeeze(mean(mean(nReshapedStim(iStimLoc).nPostStim(:, :, :, 1:end-1), 4), 2));
	nMeanStim(iStimLoc).nStim 		= mean(mean(nReshapedStim(iStimLoc).nStim, 3), 2);
	nMeanStim(iStimLoc).nAll 		= horzcat(nMeanStim(iStimLoc).nPreStim, nMeanStim(iStimLoc).nStim, nMeanStim(iStimLoc).nPostStim);

	arr = nReshapedStim(iStimLoc).nPreStim; 	nMeanStim(iStimLoc).stdPreStim 	= squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshapedStim(iStimLoc).nPostStim; 	nMeanStim(iStimLoc).stdPostStim = squeeze(std(reshape(permute(arr, [1, 2, 4, 3]), [nTraces, size(arr, 2)*size(arr, 4), size(arr, 3)]), 0, 2));
	arr = nReshapedStim(iStimLoc).nStim; 		nMeanStim(iStimLoc).stdStim		= std(reshape(arr, [nTraces, size(arr, 2)*size(arr, 3)]), 0, 2);
	nMeanStim(iStimLoc).stdAll = horzcat(nMeanStim(iStimLoc).stdPreStim, nMeanStim(iStimLoc).stdStim, nMeanStim(iStimLoc).stdPostStim);
end

% iii)	Plot dat
if isempty(forceStackHeight)
	nMeanStimStackHeight = maxel([nMeanStim.nAll]+[nMeanStim.stdAll]);
else
	nMeanStimStackHeight = forceStackHeight;
end
nMeanStimStacked = nMeanStim;

t = secondsPerBin*horzcat((0.5 - size2(nMeanStim(1).nPreStim)):1:-0.5, [0], 0.5:1:(size2(nMeanStim(1).nPostStim) - 0.5));

for iStimLoc = 1:nStimLocs
	nMeanStimStacked(iStimLoc).nAll = nMeanStimStacked(iStimLoc).nAll + repmat((0:(nTraces - 1))', [1, size2(nMeanStimStacked(1).nAll)]).*nMeanStimStackHeight;
end

figWidth = 0.1;

for iStimLoc = 1:nStimLocs
	figure('units', 'normalized', 'outerposition', [min(figWidth*(nStimLocs*expId + iStimLoc - 1), 1 - figWidth), 0, figWidth, 1]), hold on
	for iTrace = 1:nTraces
		% shadedErrorBar(t, nMeanStimStacked(iStimLoc).nAll(iTrace, :), nMeanStim(iStimLoc).stdAll(iTrace, :));			% Mean trace
		plot(t, nMeanStimStacked(iStimLoc).nAll(iTrace, :), 'k');			% Mean trace
		% plot(t, min(nMeanStimStackHeight*(iTrace - 1) + 1, nMeanStimStacked(iStimLoc).nAll(iTrace, :) + nMeanStim(iStimLoc).stdAll(iTrace, :)), 'k:');
		% plot(t, max(nMeanStimStackHeight*(iTrace - 1), nMeanStimStacked(iStimLoc).nAll(iTrace, :) - nMeanStim(iStimLoc).stdAll(iTrace, :)), 'k:');
		plot(t, nMeanStimStackHeight*(iTrace - 1) + zeros(size(t)), 'k:');
	end
	plot(t, nMeanStimStackHeight*iTrace + zeros(size(t)), 'k:');
	line([0, 0],[0, (1 + nTraces)*nMeanStimStackHeight], 'Color', 'k', 'LineStyle', ':');
	hold off
	xlabel('Time (s)')
	ylabel(['Neurons/Mean Spike Probability (Grid height = ', num2str(nMeanStimStackHeight), ')'])
	title([expName, num2str(iStimLoc)])
	ylim([0, (1 + nTraces)*nMeanStimStackHeight])
	xlim([t(1), t(end)])
	ax = gca;
	ax.YTick = nMeanStimStackHeight*(0.5:5:(nTraces - 0.5));	
	ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
	ax.TickLength = [0 0];
end

