function [n, P, nMean, nMeanBinned, nMeanBinnedBaseline, nVarBinnedBaseline] = plotStimResults(obj, dF, scopeStimArtefact, stimLog, n, P, spikeThreshold, nBins, stimsPerTrain, mirrorPosList)
% function [n P] = plotStimResults(obj, dF, scopeStimArtefact, stimLog, n, P, plotMethod, spikeThreshold, nBins)


% Separate stimulations based on mirror position
if isempty(stimsPerTrain)
    stimsPerTrain = 7;
end
if isempty(mirrorPosList)
    mirrorPosList = {...
        [-1.2,-0.8],...
        [2.8,3.2]...
    };
end
colors = 'rg';

scopeStimParams = scopeStimArtefact;

for iStim = 1:length(stimLog)
	mirrorPos = stimLog(iStim).mirrorPos;
	iPulses = ((iStim - 1)*stimsPerTrain + 1):iStim*stimsPerTrain;
	scopeStimParams(iPulses, 2:3) = repmat(mirrorPos, length(iPulses), 1);
	scopeStimParams(iPulses, 4) = 0;
	for iMirrorPos = 1:length(mirrorPosList)
		if isequal(single(mirrorPos), single(mirrorPosList{iMirrorPos}))
			scopeStimParams(iPulses, 4) = iMirrorPos;
		end
	end		
end

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
t = obj.metaDataSI.SI4.scanFramePeriod * (1:size(dF, 2));

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
if isempty(n) || isempty (P)
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


%% 4. Plot spike timing (relative to stimulation offset) histogram
% We can use the deconvoluted spike probability traces for this
nSorted = n;
nSorted = horzcat(nSorted, max(dF, [], 2));
nSorted = sortrows(nSorted, size(nSorted, 2));
nSorted = nSorted(:,1:(end - 1));

% Group into bins and take the average

% Binning Method: Mean spiking probability in bin (spike probability for each bin)
for iTrace = 1:size1(nSorted)
	[nMean(:, :, iTrace), nVar(:, :, iTrace)] = averageStimTrials(nSorted(iTrace, :), stimLog, scopeStimParams);
	[nMeanBaseline(iTrace, :), nVarBaseline(iTrace, :)] = averageBaselineTrials(nSorted(iTrace, :), stimLog, scopeStimParams);
end

if isempty(nBins)
	nBins = 20;
end
binSize = floor(size2(nMean)/nBins);

for iBin = 1:nBins
	p = nMean(:, ((iBin - 1)*binSize + 1):iBin*binSize, :);
	nMeanBinned(:, iBin, :) = mean(p, 2);

	vars = nVar(:, ((iBin - 1)*binSize + 1):iBin*binSize);
	nVarBinned = sum(vars, 2)./(binSize^2); % Variance is additive assuming independence between spike probabilities at different timepoints.

	pBaseline = nMeanBaseline(:, ((iBin - 1)*binSize + 1):iBin*binSize);
	nMeanBinnedBaseline(:, iBin) = mean(pBaseline, 2);

	varBaseline = nVarBaseline(:, ((iBin - 1)*binSize + 1):iBin*binSize);
	nVarBinnedBaseline = sum(varBaseline, 2)./(binSize^2); % Variance is additive assuming independence between spike probabilities at different timepoints.
end

% Stack the binned mean spiking probabilities for each cell
% nMeanBinnedStackHeight = 0.5;
nMeanBinnedStackHeight = 1*max(nMeanBinned(:));
nMeanBinnedStacked = nMeanBinned;
nMeanBinnedStackedBaseline = nMeanBinnedBaseline;
tBinned = obj.metaDataSI.SI4.scanFramePeriod*(1:binSize:nBins*binSize);

for iTrace = 1:size(nMeanBinned, 3)
	nMeanBinnedStacked(:, :, iTrace) = nMeanBinnedStacked(:, :, iTrace) + (iTrace - 1)*nMeanBinnedStackHeight;
	nMeanBinnedStackedBaseline(iTrace, :) = nMeanBinnedStackedBaseline(iTrace, :) + (iTrace - 1)*nMeanBinnedStackHeight;
end

figWidth = 0.1;

for iStimLoc = 1:size(nMeanBinned, 1)
	figure('units', 'normalized', 'outerposition', [min(figWidth*(iStimLoc - 1), 1 - figWidth), 0, figWidth, 1]), hold on
	for iTrace = 1:size(nMeanBinned, 3)
		plot(tBinned, nMeanBinnedStacked(iStimLoc, :, iTrace), 'k');
		plot(tBinned, nMeanBinnedStackedBaseline(iTrace, :), 'k--');
		plot(tBinned, nMeanBinnedStackHeight*(iTrace - 1) + zeros(size(tBinned)), 'k:');
	end
	plot(tBinned, iTrace + zeros(size(tBinned)), 'k:');
	hold off
	xlabel('Time After Stimulation (s)')
	ylabel('Neurons/Mean Spike Probability Across Trials')
	title(['Mean Traces Across Trials, Stim Loc ', num2str(iStimLoc)])
	ylim([0, (1 + size(nMeanBinned, 3))*nMeanBinnedStackHeight])
	xlim([0, maxel(tBinned)])
end
