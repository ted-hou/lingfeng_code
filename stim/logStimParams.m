sync_weak.time = sync_weak.time - sync_weak.time(1);
sync_weak.plot
sync_weak.parse(14, 2100);


stimParams_weak.mirrorPosList = {[1.8, 3.4], [-1.2, 0.4]};
stimParams_weak.scopeStimArtefact = sync_weak.p.ind.scopeStimArtefact;
stimParams_weak.stimLog = sync_weak.stimLog;
stimParams_weak.stimsPerTrain = 4;
stimParams_weak.scanFrameRate = HLF001.metaDataSI.SI4.scanFrameRate;



[n_strong, P_strong] = plotStimResults(dF_strong, 0.7,[], [], stimParams_strong, 'Strong', 0, 0.1);
[n_medium, P_medium] = plotStimResults(dF_medium, 0.7, [], [], stimParams_medium, 'Medium', 1, 0.1);
[n_weak, P_weak] = plotStimResults(dF_weak, 0.7, [], [], stimParams_weak, 'Weak', 2, 0.1);

[~, ~] = plotStimResults(dF_strong, 0.7, n_strong, P_strong, stimParams_strong, 'Strong', 0, 0.1);
[~, ~] = plotStimResults(dF_medium, 0.7, n_medium, P_medium, stimParams_medium, 'Medium', 1, 0.1);
[~, ~] = plotStimResults(dF_weak, 0.7, n_weak, P_weak, stimParams_weak, 'Weak', 2, 0.1);