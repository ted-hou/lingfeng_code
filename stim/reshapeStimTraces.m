%% reshapeTraces: Reshape the traces array n(iTrace, iFrame) into nReshaped(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex)
function [nReshaped] = reshapeStimTraces(n, timeRadius, secondsPerBin, stimParams)

% Params for extracting trials
scopeStimParams 	= getScopeStimParams(stimParams);
stimLog 			= stimParams.stimLog;
stimsPerTrain		= stimParams.stimsPerTrain;
stimOnsetFrames 	= scopeStimParams(1:stimsPerTrain:end, 1);							% Trial end
stimOffsetFrames	= scopeStimParams(stimsPerTrain:stimsPerTrain:end, 1);				% Trial start
frameRadius			= floor(timeRadius*stimParams.scanFrameRate);					% This many frames before the stim train and this many frames after the stim train will be extracted
framesPerBin 		= floor(secondsPerBin*stimParams.scanFrameRate);
nBinsHalf			= floor(frameRadius/framesPerBin);
framesStim			= stimOffsetFrames(1) - stimOnsetFrames(1) + 1;
framesPerTrial		= 2*frameRadius + framesStim;
nTrials				= length(stimLog);


% Separate according to stimulated location
stimLocationIds = unique(nonzeros(scopeStimParams(:, 4)));
trialIndices 	= zeros(numel(stimLocationIds), 1);

% Separate into trials
for iTrial = 1:nTrials
	iStimLoc 				= scopeStimParams(iTrial*stimsPerTrain, 4);
	trialIndices(iStimLoc) 	= trialIndices(iStimLoc) + 1; 
	iTrialForIndexing 		= trialIndices(iStimLoc);

	nReshaped(iStimLoc).nStim(:, :, iTrialForIndexing) = ...
		n(:, stimOnsetFrames(iTrial):stimOffsetFrames(iTrial));

	for iBin = 1:nBinsHalf
		nReshaped(iStimLoc).nPreStim(:, :, iBin, iTrialForIndexing) = ...
			n(:, (stimOnsetFrames(iTrial) - (nBinsHalf - iBin + 1)*framesPerBin):(stimOnsetFrames(iTrial) - (nBinsHalf - iBin)*framesPerBin - 1));
		nReshaped(iStimLoc).nPostStim(:, :, iBin, iTrialForIndexing) = ...
			n(:, (stimOffsetFrames(iTrial) + (iBin - 1)*framesPerBin + 1):(stimOffsetFrames(iTrial) + iBin*framesPerBin));
	end

end

