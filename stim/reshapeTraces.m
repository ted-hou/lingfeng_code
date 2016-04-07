%% reshapeTraces: Reshape the traces array n(iTrace, iFrame) into nReshaped(iStimLoc).n(iTrace, iFrame, iBin, iTrialNewIndex)
function [nReshaped, nReshapedBaseline] = reshapeTraces(n, nBins, stimParams)

% Params for extracting trials
scopeStimParams 	= getScopeStimParams(stimParams);
stimLog 			= stimParams.stimLog;
stimsPerTrain		= stimParams.stimsPerTrain;
stimProtocolStart 	= scopeStimParams(1, 1);
stimProtocolEnd 	= scopeStimParams(end, 1);
stimOnsetFrames 	= scopeStimParams(1:stimsPerTrain:end, 1);							% Trial end
stimOffsetFrames	= scopeStimParams(stimsPerTrain:stimsPerTrain:end, 1);				% Trial start
framesPerTrial 		= minel(stimOnsetFrames(2:end) - stimOffsetFrames(1:end - 1)) - 1;	% Use the minimal ISI as standard trial length, for longer trials, events after framesPerTrial is discarded
framesPerBin 		= floor(framesPerTrial/nBins);
nTrials				= length(stimLog);
nTrialsPreStim		= floor((stimProtocolStart - 1)/framesPerTrial);
nTrialsPostStim		= floor((size2(n) - stimProtocolEnd)/framesPerTrial) - 1;			% -1 to remove the first trial immediately after the final stimulation train



% Separate according to stimulated location
stimLocationIds = unique(nonzeros(scopeStimParams(:, 4)));
for iStimLoc = stimLocationIds'
	nReshaped(iStimLoc).n = [];
end

% Separate into trials
for iTrial = 1:nTrials
	iStimLoc 			= scopeStimParams(iTrial*stimsPerTrain, 4);
	iTrialForIndexing 	= size(nReshaped(iStimLoc).n, 4) + 1;

	% Separate into timebins
	for iBin = 1:nBins
		nReshaped(iStimLoc).n(:, 1:framesPerBin, iBin, iTrialForIndexing) = ...
			n(:, (stimOffsetFrames(iTrial) + 1 + (iBin - 1)*framesPerBin):(stimOffsetFrames(iTrial) + iBin*framesPerBin));
	end
end


% Now let's reshape the baseline traces n(iTrace, iFrame) into nReshapedBaseline(iTrace, iFrame, iBin, iTrial)
% First segment pre-stim baseline traces into trials
for iTrial = 1:nTrialsPreStim
	for iBin = 1:nBins
		nReshapedPreStim(:, 1:framesPerBin, iBin, iTrial) = ...
			n(:, (1 + (iTrial - 1)*framesPerTrial + (iBin - 1)*framesPerBin):((iTrial - 1)*framesPerTrial + iBin*framesPerBin));
	end
end

% Then segment post-stim baseline traces into trials
for iTrial = 1:nTrialsPostStim
	for iBin = 1:nBins
		nReshapedPostStim(:, 1:framesPerBin, iBin, iTrial) = ...
			n(:, (1 + stimProtocolEnd + iTrial*framesPerTrial + (iBin - 1)*framesPerBin):(stimProtocolEnd + iTrial*framesPerTrial + iBin*framesPerBin));
	end
end

% Concatenate reshaped pre/post-stim baseline traces
nReshapedBaseline = cat(4, nReshapedPreStim, nReshapedPostStim);
