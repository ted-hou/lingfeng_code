function [constrainedFoopsi, spikesReshaped, spikesMean] = plotStimResults(dF, constrainedFoopsi, stimParams, timeRadius, expName, expId, forceStackHeight)
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
mirrorPosList		= stimParams.mirrorPosList;
scopeFramePeriod	= 1/stimParams.scanFrameRate;
nTraces 			= size1(dF);
nFramesPerTrace		= size2(dF);
nStimLocs 			= length(mirrorPosList);
%% 1. Remove stim artefact, interpolate for stim frames
dF = removeStimArtefact(dF, scopeStimArtefact);

dFStackHeight = 1*maxel(dF);

dFStacked = dF;
for iTrace = 1:nTraces
	dFStacked(iTrace, :) = dFStacked(iTrace, :) + (iTrace - 1)*dFStackHeight;
end


% time in seconds
t = scopeFramePeriod * (1:nFramesPerTrace);

% Plot all traces + stim times
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
hold on;
for iMirrorPos = 1:nStimLocs;
	stimTrace = zeros(1, nFramesPerTrace);
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + nTraces)*dFStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end

plot(t, dFStacked', 'k'), hold off;

ylim([0, (1 + nTraces)*dFStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/GCaMP Signal (\DeltaF/F)')
title('\DeltaF/F')
ax = gca;
ax.YTick = dFStackHeight*(0.5:5:(nTraces - 0.5));	
ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
ax.TickLength = [0 0];


%% 2. Deconvolve the traces using constrained foopsi from Volgostein
if isempty(constrainedFoopsi) % Skipped if deconvolution already done
    constrainedFoopsi = deconvolve(dF);
end

calcium = constrainedFoopsi.c;
spikes = constrainedFoopsi.deconvInSpikes;

% Plot decovoluted calcium soncentrations
calciumStackHeight = max(calcium(:));
calciumStacked = calcium;
for iTrace = 1:nTraces
	calciumStacked(iTrace, :) = calciumStacked(iTrace, :) + (iTrace - 1)*calciumStackHeight;
end

figure('units', 'normalized', 'outerposition', [0 0 1 1])
hold on;
for iMirrorPos = 1:nStimLocs;
	stimTrace = zeros(1, nFramesPerTrace);
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + nTraces)*calciumStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end
plot(t, calciumStacked', 'k'), hold off;

ylim([0, (1 + nTraces)*calciumStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/Deconvoluted Calcium Concentrations')
title('Deconvoluted Calcium Concentrations')
ax = gca;
ax.YTick = calciumStackHeight*(0.5:5:(nTraces - 0.5));	
ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
ax.TickLength = [0 0];

% Plot spike rate
spikesStackHeight = max(spikes(:));
spikesStacked = horzcat(spikes, max(dF, [], 2));
spikesStacked = sortrows(spikesStacked, size(spikesStacked, 2));
spikesStacked = spikesStacked(:,1:(end - 1));

for iTrace = 1:nTraces
	spikesStacked(iTrace, :) = spikesStacked(iTrace, :) + (iTrace - 1)*spikesStackHeight;
end

figure('units', 'normalized', 'outerposition', [0 0 1 1])
hold on;
for iMirrorPos = 1:nStimLocs;
	stimTrace = zeros(1, nFramesPerTrace);
	frames = scopeStimArtefact(scopeStimParams(:, 4)==iMirrorPos);
	stimTrace(frames') = (1 + nTraces)*spikesStackHeight;
	plot(t, stimTrace', [colors(iMirrorPos),'-.'])
end
plot(t, spikesStacked', 'k'), hold off;

ylim([0, (1 + nTraces)*spikesStackHeight])
xlim([0, maxel(t)])
xlabel('Time (s)')
ylabel('Neurons/Spike Rate')
title('Deconvoluted Spike Rate')
ax = gca;
ax.YTick = spikesStackHeight*(0.5:5:(nTraces - 0.5));	
ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
ax.TickLength = [0 0];

%% 4. Plot aligned stim traces
% i) 	Extract deconvoluted traces from t-timeRadius till t+timeRadius, (t=stim time), reshaped into spikesReshaped(iStimLoc).calcium(iTrace, iFrame, iBin, iTrialNewIndex)
% ii) 	Calculate the mean

% i)	Extract traces 5 seconds before and after stimulation
secondsPerBin = scopeFramePeriod;
% secondsPerBin = 0.1;
spikesReshaped = reshapeTraces(spikes, timeRadius, secondsPerBin, stimParams);

% ii)	Calculate mean traces for each ROI
spikesMean = getMeanTraceAcrossTrials(spikesReshaped);

% iii)	Plot dat
if isempty(forceStackHeight)
	spikesMeanStackHeight = maxel([spikesMean.all]+[spikesMean.stdAll]);
else
	spikesMeanStackHeight = forceStackHeight;
end
spikesMeanStacked = spikesMean;

t = secondsPerBin*horzcat((0.5 - size2(spikesMean(1).preStim)):1:-0.5, [0], 0.5:1:(size2(spikesMean(1).postStim) - 0.5));

for iStimLoc = 1:nStimLocs
	spikesMeanStacked(iStimLoc).all = spikesMeanStacked(iStimLoc).all + repmat((0:(nTraces - 1))', [1, size2(spikesMeanStacked(1).all)]).*spikesMeanStackHeight;
end

figWidth = 0.1;

for iStimLoc = 1:nStimLocs
	figure('units', 'normalized', 'outerposition', [min(figWidth*(nStimLocs*expId + iStimLoc - 1), 1 - figWidth), 0, figWidth, 1]), hold on
	for iTrace = 1:nTraces
		% shadedErrorBar(t, spikesMeanStacked(iStimLoc).all(iTrace, :), spikesMean(iStimLoc).stdAll(iTrace, :));	% Mean trace with errorbar (std)
		plot(t, spikesMeanStacked(iStimLoc).all(iTrace, :), 'k');			% Mean trace
		plot(t, spikesMeanStackHeight*(iTrace - 1) + zeros(size(t)), 'k:');
	end
	plot(t, spikesMeanStackHeight*iTrace + zeros(size(t)), 'k:');
	line([0, 0],[0, (1 + nTraces)*spikesMeanStackHeight], 'Color', 'k', 'LineStyle', ':');
	hold off
	xlabel('Time (s)')
	ylabel(['Cell ID/Mean Spike Rate (Grid height = ', num2str(spikesMeanStackHeight), ')'])
	title([expName, num2str(iStimLoc)])
	ylim([0, (1 + nTraces)*spikesMeanStackHeight])
	xlim([t(1), t(end)])
	ax = gca;
	ax.YTick = spikesMeanStackHeight*(0.5:5:(nTraces - 0.5));	
	ax.YTickLabel = cellfun(@num2str, num2cell(1:5:nTraces), 'UniformOutput', false);
	ax.TickLength = [0 0];
end

