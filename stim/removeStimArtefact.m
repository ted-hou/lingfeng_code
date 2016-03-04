function dFNoArtefact = removeStimArtefact(dF, scopeStimArtefact)

dFNoArtefact = dF;

for iFrame = scopeStimArtefact'
    dFNoArtefact(:, iFrame) = mean(dF(:, [max(1, iFrame - 1), min(size(dF, 2), iFrame + 1)]), 2);
end
