%% reshapeTraces: Reshape the traces array n(iTrace, iFrame) into nReshaped(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex)
function [nReshaped] = reshapeStimTraces(n, secondsPerBin, stimParams)

% Params for extracting trials
% scopeStimParams 	= getScopeStimParams(stimParams);
% stimLog 			= stimParams.stimLog;
% stimsPerTrain		= stimParams.stimsPerTrain;
% stimOnsetFrames 	= scopeStimParams(1:stimsPerTrain:end, 1);							% Trial end
% stimOffsetFrames	= scopeStimParams(stimsPerTrain:stimsPerTrain:end, 1);				% Trial start
% frameRadius		= floor(timeRadius*stimParams.scanFrameRate);					% This many frames before the stim train and this many frames after the stim train will be extracted
% framesPerBin 		= floor(secondsPerBin*stimParams.scanFrameRate);
% nBinsHalf			= floor(frameRadius/framesPerBin);
% framesStim			= stimOffsetFrames(1) - stimOnsetFrames(1) + 1;
% framesPerTrial		= 2*frameRadius + framesStim;
% nTrials				= length(stimLog);

scopeStimParams 	= getScopeStimParams(stimParams);
stimLog 			= stimParams.stimLog;
stimsPerTrain		= stimParams.stimsPerTrain;
stimOnsetFrames 	= scopeStimParams(1:stimsPerTrain:end, 1);							% Trial end
stimOffsetFrames	= scopeStimParams(stimsPerTrain:stimsPerTrain:end, 1);				% Trial start
nTrials				= length(stimLog);
framesPerBin 		= floor(secondsPerBin*stimParams.scanFrameRate);

stimLocationIds = unique(nonzeros(scopeStimParams(:, 4)));
trialIndices = zeros(numel(stimLocationIds), 1);

% Separate into trials
for iTrial = 1:nTrials
	iStimLoc = scopeStimParams(iTrial*stimsPerTrain, 4);
	trialIndices(iStimLoc) = trialIndices(iStimLoc) + 1;
	iTrialForIndexing = trialIndices(iStimLoc);

	% A trial starts after stim ends, and ends before the next stim begins
	trialStartFrame = stimOffsetFrames(iTrial) + 1;
	if iTrial < nTrials
		trialEndFrame	= stimOnsetFrames(iTrial + 1) - 1;
		framesThisTrial = trialEndFrame - trialStartFrame + 1;
	else
		trialEndFrame	= trialStartFrame + framesThisTrial - 1;						% For the last trial, use the length of the previous trial
	end

	% Perform standardized binning for current trial: all frames in this trial is divided into many bins of equal length (length defined by input param). If this is not possible, the bin with extra/less frames should be in the center of this ITI.
	nBins = round(framesThisTrial/framesPerBin);
	extraFrames = framesThisTrial - nBins*framesPerBin; % This value is negative for a shorter-than-standard bin. e.g., binSize = 3; nFrames = 100; then nBins = 33; extraFrames = 100 - 99 = 1; So the center(ish) bin should have size 4 instead of 3. And yes I had to write this out.
	iCenterBin = ceil((nBins + 1)/2);

	nThisTrial = n(:, trialStartFrame:trialEndFrame);

	% Group frames during stimulation into one bin
	nReshaped(iStimLoc).n(:, :, 1, iTrialForIndexing) = n(:, stimOnsetFrames(iTrial):stimOffsetFrames(iTrial));

	lastFrame = 0;
	for iBin = 1:nBins
		if iBin == iCenterBin
			binTrace = nThisTrial(:, (lastFrame + 1):(lastFrame + framesPerBin + extraFrames));
			if extraFrames < 0
				binTrace = horzcat(binTrace, nan(size(binTrace, 1), -extraFrames));
			end
			nReshaped(iStimLoc).n(:, 1:size(binTrace ,2), iBin + 1, iTrialForIndexing) = binTrace;
			lastFrame = lastFrame + framesPerBin + extraFrames;
		else
			binTrace = nThisTrial(:, (lastFrame + 1):(lastFrame + framesPerBin));
			if extraFrames > 0
				binTrace = horzcat(binTrace, nan(size(binTrace, 1), extraFrames));
			end
			nReshaped(iStimLoc).n(:, 1:size(binTrace ,2), iBin + 1, iTrialForIndexing) = binTrace;
			lastFrame = lastFrame + framesPerBin;
		end
	end
end

