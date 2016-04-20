function constrainedFoopsi = simpleDeconvolve(dF, stimParams)
dF = removeStimArtefact(dF, stimParams.scopeStimArtefact);
constrainedFoopsi = deconvolve(dF);

calcium = constrainedFoopsi.c;
spikes = constrainedFoopsi.deconvInSpikes;
