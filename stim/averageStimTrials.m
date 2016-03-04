%% averageStimTrials: Get the mean response of a trace to repeated stimulation trials
function [dFMean dFVar dFTrial] = averageStimTrials(dF, stimLog, scopeStimParams)

% Get the minimal ITI among all trials, discard events occuring after minInterval for each trial.
nTrials			= size2(stimLog);
nPulsesPerTrain = size1(scopeStimParams)./nTrials;

if isinteger(nPulsesPerTrain)
	error('Number of pulses per train must be an integer!');
end

stimOnsetFrames 	= scopeStimParams(1:nPulsesPerTrain:end, 1);
stimOffsetFrames	= scopeStimParams(nPulsesPerTrain:nPulsesPerTrain:end, 1);

minInterval = minel(stimOnsetFrames(2:end) - stimOffsetFrames(1:end - 1));

% Separate trials based on stimulated location
stimLocationIds = unique(nonzeros(scopeStimParams(:, 4)));

for i = stimLocationIds'
	dFTrial(i).dF = [];
end

for iTrial = 1:nTrials
	dFTrial(scopeStimParams(iTrial*nPulsesPerTrain, 4)).dF(size1(dFTrial(scopeStimParams(iTrial*nPulsesPerTrain, 4)).dF) + 1, :) = dF(stimOffsetFrames(iTrial) + 1:stimOffsetFrames(iTrial) + minInterval - 1);
end


% Calculate mean/variance of trial responses
for iStimLoc = stimLocationIds'
	dFMean(iStimLoc, 1:size2(dFTrial(iStimLoc).dF)) = mean(dFTrial(iStimLoc).dF, 1);
	dFVar(iStimLoc, 1:size2(dFTrial(iStimLoc).dF)) = var(dFTrial(iStimLoc).dF, 1);
end
