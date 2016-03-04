%% getScopeStimparams
function scopeStimParams = getScopeStimParams(stimParams)

scopeStimParams = stimParams.scopeStimArtefact;
stimLog 		= stimParams.stimLog;
stimsPerTrain	= stimParams.stimsPerTrain;
mirrorPosList	= stimParams.mirrorPosList;

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
