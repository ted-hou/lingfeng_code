function [SVMModel, xTrain, xTest, yTrain, yTest, meanCorrect] = stimSvmTest(spikesReshaped, stimParams, trainingInterval);

nStimLocs 			= length(spikesReshaped);
nTrials				= size(spikesReshaped(1).stim, 3);
nTraces 			= size(spikesReshaped(1).stim, 1);
framesPerBin		= size(spikesReshaped(1).postStim, 2); % Ideally there should be one frame per bin.
binsPerTrial 		= size(spikesReshaped(1).postStim, 3);
scopeFramePeriod 	= 1/stimParams.scanFrameRate;

startTimeForTraining = trainingInterval(1);
timeStimOff 		 = size(spikesReshaped(1).stim, 2)*scopeFramePeriod;
endTimeForTraining 	 = timeStimOff + trainingInterval(2);

firstBinForTraining	= max(1, ceil(startTimeForTraining/(scopeFramePeriod*framesPerBin)));
lastBinForTraining  = floor(endTimeForTraining/(scopeFramePeriod*framesPerBin));
% e.g 30 trials as training set, 15 trials as test set
% Predictor params = activity of each cell
% Observation = population activity at a single frame
% Category = stim location
fractionAsTrainingSet 	= 1;
nTrialsAsTrainingSet 	= round(nTrials*fractionAsTrainingSet);
% nTrialsAsTestingSet 	= nTrials - nTrialsAsTrainingSet;
nTrialsAsTestingSet 	= nTrials;

randId					= randperm(nTrials);
iTrialsAsTrainingSet	= randId(1:nTrialsAsTrainingSet);
% iTrialsAsTestingSet		= randId(nTrialsAsTrainingSet + 1:end);
iTrialsAsTestingSet		= 1:nTrials;

xTrain = [];
xTest  = [];
yTrain = [];
yTest  = [];

for iStimLoc = 1:nStimLocs
	spikesPostStim 	= squeeze(mean(spikesReshaped(iStimLoc).postStim, 2));					% Convert to 3D array. If binned, get mean of bin.
	spikesStim 		= spikesReshaped(iStimLoc).stim;											% Include activity during stimulation. Only supports framesPerBin = 1;
	spikes 			= horzcat(spikesStim, spikesPostStim);

	spikesTrain 	= spikes(:, firstBinForTraining:lastBinForTraining, iTrialsAsTrainingSet);		% Randomly pull some trials for training; Only keep the first few seconds post-stim for training
	spikesTrain 	= reshape(spikesTrain, [nTraces, size(spikesTrain, 2)*size(spikesTrain, 3)]);	% Convert to 2D array. Discard trial ID.
	xTrain 			= vertcat(xTrain, spikesTrain');												% Predictor array
	yTrain 			= vertcat(yTrain, iStimLoc*ones(size(spikesTrain, 2), 1));						% Response array

	spikesTest 		= spikes(:, :, iTrialsAsTestingSet);																% Trials not used for training are used for testing
	spikesTest 		= reshape(spikesTest, [nTraces, size(spikesTest, 2)*size(spikesTest, 3)]);		% Convert to 2D array. Discard trial ID.
	xTest			= vertcat(xTest, spikesTest');													% Predictor array
	yTest 			= vertcat(yTest, iStimLoc*ones(size(spikesTest, 2), 1));						% Response array
end



SVMModel = fitcsvm(xTrain, yTrain, 'KernelFunction', 'rbf');
label = SVMModel.predict(xTrain); correct = label == yTrain; sum(correct)/numel(yTrain)
label = SVMModel.predict(xTest); correct = label == yTest; 
sum(correct)/numel(yTest)

% Do a prediction for each trial. Calculate prediction accuracy for each frame in trial.
correct 	= reshape(correct, [numel(correct)/(nTrialsAsTestingSet*nStimLocs), nTrialsAsTestingSet*nStimLocs]);
stdCorrect 	= squeeze(std(correct, 0, 2));
meanCorrect = squeeze(mean(correct, 2));
hCorrect	= ttest(correct', 0.5*ones(size(correct')), 'Alpha', 0.01);

t = (1:size(meanCorrect, 1))'*scopeFramePeriod;

figure('units', 'normalized', 'outerposition', [0, 0.7, 1, 0.3])
hold on;
xmarkers = t(find(hCorrect==1));
ymarkers = meanCorrect(find(hCorrect==1));
scatter(xmarkers, ymarkers, 15, 'k*');
shadedErrorBar(t, meanCorrect, stdCorrect, '-k', 1);
line([startTimeForTraining, startTimeForTraining], [-0.5 1.5], 'LineWidth', 2);
line([endTimeForTraining, endTimeForTraining], [-0.5 1.5], 'LineWidth', 2);
line([timeStimOff, timeStimOff], [-0.5, 1.5], 'LineWidth', 2, 'Color', 'k');
line([0, size(meanCorrect, 1)*scopeFramePeriod], [0.5, 0.5], 'LineWidth', 2);
hold off;
xlabel('Time after stim on (s)');
ylabel('Classification accuracy');

% Also plot two stimLocs separately just for fun
colorList = 'rgbc';
hold on;
correctReshaped = reshape(correct, [numel(correct)/(nTrialsAsTestingSet*nStimLocs), nTrialsAsTestingSet, nStimLocs]);
for iStimLoc = 1:nStimLocs
	plot(t, mean(correctReshaped(:, :, iStimLoc), 2), [colorList(iStimLoc), '--']);
end
hold off;

% Shuffle stimLoc labels
yTrainShuffled = yTrain(randperm(numel(yTrain)));

SVMModel = fitcsvm(xTrain, yTrainShuffled, 'KernelFunction', 'linear');
label = SVMModel.predict(xTrain); correct = label == yTrainShuffled; sum(correct)/numel(yTrainShuffled)
label = SVMModel.predict(xTest); correct = label == yTest; 
sum(correct)/numel(yTest)

% Do a prediction for each trial. Calculate prediction accuracy for each frame in trial.
correct 	= reshape(correct, [numel(correct)/(nTrialsAsTestingSet*nStimLocs), nTrialsAsTestingSet*nStimLocs]);
stdCorrect 	= squeeze(std(correct, 0, 2));
meanCorrect = squeeze(mean(correct, 2));
hCorrect	= ttest(correct', 0.5*ones(size(correct')), 'Alpha', 0.001);

t = (1:size(meanCorrect, 1))'*scopeFramePeriod;

figure('units', 'normalized', 'outerposition', [0, 0.4, 1, 0.3])
% plot((1:size(meanCorrect, 1))'*scopeFramePeriod, meanCorrect);
hold on
xmarkers = t(find(hCorrect==1));
ymarkers = meanCorrect(find(hCorrect==1));
scatter(xmarkers, ymarkers, 15, 'k*');
shadedErrorBar(t, meanCorrect, stdCorrect, '-k', 1);

line([startTimeForTraining, startTimeForTraining], [-0.5 1.5], 'LineWidth', 2);
line([endTimeForTraining, endTimeForTraining], [-0.5 1.5], 'LineWidth', 2);
line([timeStimOff, timeStimOff], [-0.5, 1.5], 'LineWidth', 2, 'Color', 'k');
line([0, size(meanCorrect, 1)*scopeFramePeriod], [0.5, 0.5], 'LineWidth', 2);
hold off
xlabel('Time after stim on (s)');
ylabel('Classification accuracy');

% Also plot two stimLocs separately just for fun
colorList = 'rgbc';
hold on;
correctReshaped = reshape(correct, [numel(correct)/(nTrialsAsTestingSet*nStimLocs), nTrialsAsTestingSet, nStimLocs]);
for iStimLoc = 1:nStimLocs
	plot(t, mean(correctReshaped(:, :, iStimLoc), 2), [colorList(iStimLoc), '--']);
end
hold off;
