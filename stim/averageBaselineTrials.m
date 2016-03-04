%% averageBaselineTrials: Get the mean baseline activity of a trace
function [dFMeanBaseline, dFVarBaseline] = averageBaselineTrials(dF, stimLog, scopeStimParams)

% Get the minimal ITI among all trials, discard events occuring after minInterval for each trial.
nTrials			= size2(stimLog);
nPulsesPerTrain = size1(scopeStimParams)/nTrials;

stimProtocolStart 	= scopeStimParams(1, 1);
stimProtocolEnd 	= scopeStimParams(end, 1);

stimOnsetFrames 	= scopeStimParams(1:nPulsesPerTrain:end, 1);
stimOffsetFrames	= scopeStimParams(nPulsesPerTrain:nPulsesPerTrain:end, 1);

minInterval = minel(stimOnsetFrames(2:end) - stimOffsetFrames(1:end - 1));
framesInterval = minInterval - 1;

% Divide baseline trace into small trials with equal length = minInterval (stimulation ISI).
% First decide how many trials we can divide the baseline trace into.
nTrialsPreStim	= floor((stimProtocolStart - 1)/framesInterval);
nTrialsPostStim	= floor((size2(dF) - stimProtocolEnd - framesInterval)/framesInterval);

% First segment pre-stim trace into trials
for iTrial = 1:nTrialsPreStim
	dFTrialPreStim(iTrial, 1:framesInterval) = dF((1 + (iTrial - 1)*framesInterval):(iTrial*framesInterval));
end

for iTrial = 1:nTrialsPostStim
	dFTrialPostStim(iTrial, 1:framesInterval) = dF((stimProtocolEnd + 1 + iTrial*framesInterval):(stimProtocolEnd + (1 + iTrial)*framesInterval));
end

dFTrialBaseline = vertcat(dFTrialPreStim, dFTrialPostStim);

% Then calculate mean trial responses
dFMeanPreStim 	= mean(dFTrialPreStim, 1);
dFMeanPostStim 	= mean(dFTrialPostStim, 1);
dFMeanBaseline 	= mean(dFTrialBaseline, 1);

dFVarBaseline = var(dFTrialBaseline);
