function [n, P, nReshaped, nReshapedBaseline, h] = plotStimResults(obj, dF, spikeThreshold, nBins, n, P, stimParams)
% If there are more than 4 stimulation positions, might need more colors. But 4 should be enough.
colors = 'rgbc';

scopeStimParams 	= getScopeStimParams(stimParams);
scopeStimArtefact 	= stimParams.scopeStimArtefact;
stimLog 			= stimParams.stimLog;
stimsPerTrain		= stimParams.stimsPerTrain;
mirrorPosList		= stimParams.stimsPerTrain;
scopeFramePeriod	= obj.metaDataSI.SI4.scanFramePeriod;

%% 1. Remove stim artefact, interpolate for stim frames
dF = removeStimArtefact(dF, scopeStimArtefact);

dFStackHeight = 0.25*maxel(dF);

dFStacked = dF;
dFStacked = horzcat(dFStacked, max(dF, [], 2));
dFStacked = sortrows(dFStacked, size(dFStacked, 2));
dFStacked = dFStacked(:,1:(size(dFStacked,2) - 1));

for iTrace = 1:size(dF, 1)
	dFStacked(iTrace, :) = dFStacked(iTrace, :) + (iTrace - 1)*dFStackHeight;
end


% time in seconds
t = scopeFramePeriod * (1:size(dF, 2));

% Plot all traces + stim times
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on;
for iMirrorPos = 1:length(mirrorPosList);
	stimTrace = zeros(1, size(dF, 2));
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + size(dF, 1))*dFStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end

plot(t, dFStacked', 'k'), hold off;

ylim([0, (1 + size(dF, 1))*dFStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/GCaMP Signal (\DeltaF/F)')
title('\DeltaF/F')


%% 2. Deconvolve the traces using oopsi (fnnd) Volgostein (2010)
if isempty(n) || isempty (P) % Skipped if deconvolution already done since this is quite slow.
	for iTrace = 1:size(dF, 1)
		[n(iTrace,1:size(dF, 2)), P{iTrace}] = fast_oopsi(dF(iTrace,:));
	end
end

% Produce a stacked version of spike probabilities for plotting
nStackHeight = 1;

nStacked = n;
nStacked = horzcat(nStacked, max(dF, [], 2));
nStacked = sortrows(nStacked, size(nStacked, 2));
nStacked = nStacked(:,1:(end - 1));

for iTrace = 1:size(dF, 1)
	nStacked(iTrace, :) = nStacked(iTrace, :) + (iTrace - 1)*nStackHeight;
end

% Plot all spike probabilities + stim times
figure('units', 'normalized', 'outerposition', [0 0 1 1])
hold on;
for iMirrorPos = 1:length(mirrorPosList);
	stimTrace = zeros(1, size(n, 2));
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + size(n, 1))*nStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end
plot(t, nStacked', 'k'), hold off;

ylim([0, (1 + size(n, 1))*nStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/Spike Probability')
title('Deconvoluted Spike Probability')

% Plot all spikes + stim times
if isempty(spikeThreshold)
	spikeThreshold = 0.7;
end
nThresholded = zeros(size(n));
nThresholded(n>spikeThreshold) = 1;
nStackedAndThresholded = horzcat(nThresholded, max(dF, [], 2));
nStackedAndThresholded = sortrows(nStackedAndThresholded, size(nStackedAndThresholded, 2));
nStackedAndThresholded = nStackedAndThresholded(:,1:(end - 1));

for iTrace = 1:size(dF, 1)
	nStackedAndThresholded(iTrace, :) = nStackedAndThresholded(iTrace, :) + (iTrace - 1)*nStackHeight;
end

figure('units', 'normalized', 'outerposition', [0 0 1 1])
hold on;
for iMirrorPos = 1:length(mirrorPosList);
	stimTrace = zeros(1, size(n, 2));
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + size(n, 1))*nStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end
plot(t, nStackedAndThresholded', 'k'), hold off;

ylim([0, (1 + size(n, 1))*nStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/Spikes')
title(['Deconvoluted Spikes (Threshold = ', num2str(spikeThreshold), ')'])


%% 4. Plot spike timing (relative to stimulation offset) using the deconvoluted spike probability traces
% i) 	Extract and reshape traces (split into bins/trials/...)
% ii) 	Get mean binned traces for plotting
% iii) 	Run some statistical tests
% iv)	Now plot dat ting

% First sort the raw traces according to max spike amplitude
nSorted = n;
nSorted = horzcat(nSorted, max(dF, [], 2));
nSorted = sortrows(nSorted, size(nSorted, 2));
nSorted = nSorted(:,1:(end - 1));

% i) Reshape the deconvolute spike probability traces n(iTrace, iFrame) into 
% 		nReshaped(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex) and 
% 		nReshapedBaseline(iTrace, iFrame, iBin, iTrial)
[nReshaped, nReshapedBaseline] = reshapeTraces(nSorted, nBins, stimParams);

nStimLocs 		= length(nReshaped);
nTraces 		= size1(nReshapedBaseline);
framesPerBin 	= size2(nReshapedBaseline);

% ii) Get binned mean traces for plotting
for iStimLoc = 1:nStimLocs
	nMean(iStimLoc).n = squeeze(mean(mean(nReshaped(iStimLoc).n, 4), 2));
end
nMeanBaseline = squeeze(mean(mean(nReshapedBaseline, 4), 2));

% iii) Perform two-sample t test between stim and baseline traces, assuming unknown and unequal variances.
for iStimLoc = 1:nStimLocs
	for iTrace = 1:nTraces
		for iBin = 1:nBins
			xStim = nReshaped(iStimLoc).n(iTrace, :, iBin, :);
			xBase = nReshapedBaseline(iTrace, :, iBin, :);
			h(iStimLoc).h(iTrace, iBin) = ttest2(xStim(:), xBase(:), 'Alpha', 0.01, 'Tail', 'right', 'Vartype', 'unequal');
		end
	end
end

% iv) Now plot dat ting
nMeanStackHeight = maxel([nMean.n]);
nMeanStacked = nMean;
nMeanStackedBaseline = nMeanBaseline;

t = scopeFramePeriod*(1:framesPerBin:nBins*framesPerBin);

for iStimLoc = 1:nStimLocs
	nMeanStacked(iStimLoc).n = nMeanStacked(iStimLoc).n + repmat((0:(nTraces - 1))', [1, nBins]).*nMeanStackHeight;
	nMeanStackedBaseline = nMeanStackedBaseline + repmat((0:(nTraces - 1))', [1, nBins]).*nMeanStackHeight;
end

figWidth = 0.1;

for iStimLoc = 1:nStimLocs
	figure('units', 'normalized', 'outerposition', [min(figWidth*(iStimLoc - 1), 1 - figWidth), 0, figWidth, 1]), hold on
	for iTrace = 1:nTraces
		xmarkers = t(find(h(iStimLoc).h(iTrace, :)==1));
		ymarkers = nMeanStacked(iStimLoc).n(iTrace, find(h(iStimLoc).h(iTrace, :)==1));
		plot(t, nMeanStacked(iStimLoc).n(iTrace, :), 'k', xmarkers, ymarkers,'k*');
		plot(t, nMeanStackedBaseline(iTrace, :), 'k--');
		plot(t, nMeanStackHeight*(iTrace - 1) + zeros(size(t)), 'k:');
	end
	plot(t, nMeanStackHeight*iTrace + zeros(size(t)), 'k:');
	hold off
	xlabel('Time After Stimulation (s)')
	ylabel('Neurons/Mean Spike Probability Across Trials')
	title(['Mean Traces Across Trials, Stim Loc ', num2str(iStimLoc)])
	ylim([0, (1 + nTraces)*nMeanStackHeight])
	xlim([0, maxel(t)])
end
