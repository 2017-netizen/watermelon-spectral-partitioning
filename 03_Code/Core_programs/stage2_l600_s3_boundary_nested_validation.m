function stage2_l600_s3_boundary_nested_validation(subset)
% Repeated nested record-level CV for full-band versus S3 boundary variants.

if nargin < 1
    subset = 'combined';
end
subset = lower(char(subset));
validSubsets = {'combined', '2025', '2026'};
if ~ismember(subset, validSubsets)
    error('subset must be combined, 2025, or 2026.');
end

rng(20260712, 'twister');
set(0, 'DefaultFigureVisible', 'off');

scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(scriptDir));
programDir = fullfile(fileparts(scriptDir), 'Dependencies');
outDir = fullfile(packageRoot, '04_Model_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
inputFile = fullfile(packageRoot, '02_Data', 'L600_Data_Consolidated_EN_V02.xlsx');
addpath(programDir);

suffix = '';
if ~strcmp(subset, 'combined')
    suffix = ['_' subset];
end
logFile = fullfile(outDir, ['S3边界重复嵌套交叉验证' suffix '_运行日志.txt']);
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);
cleanupObj = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('Nested validation started: %s\n', datestr(now, 31));
fprintf('Data subset: %s\n', subset);
[wl, Y, Xraw, originalIndex] = load_l600_submission_ssc(inputFile, 'SSC_model_spectra');
if strcmp(subset, '2025')
    keep = originalIndex > 150;
    Y = Y(keep);
    Xraw = Xraw(keep, :);
    originalIndex = originalIndex(keep);
elseif strcmp(subset, '2026')
    keep = originalIndex <= 150;
    Y = Y(keep);
    Xraw = Xraw(keep, :);
    originalIndex = originalIndex(keep);
end

fullMask = 1:numel(wl);
rangeNames = {'Full_650_950', 'S3_Central_748_826', 'S3_Outer4_744_830', ...
    'S3_Inner4_752_822', 'S3_Outer10_738_836', 'S3_Inner10_758_816'};
rangeMasks = {fullMask, find(wl >= 748 & wl <= 826), ...
    find(wl >= 744 & wl <= 830), find(wl >= 752 & wl <= 822), ...
    find(wl >= 738 & wl <= 836), find(wl >= 758 & wl <= 816)};
rangeSpectra = cell(size(rangeMasks));
for rangeId = 1:numel(rangeMasks)
    % Preprocess each interval independently to preserve deployable range isolation.
    rangeSpectra{rangeId} = firstDerivativeLocal( ...
        snvLocal(Xraw(:, rangeMasks{rangeId})), 7);
end
featureNames = {'AllBands', 'CARS'};
pipelineNames = {'Full_AllBands', 'S3_Central_AllBands', 'S3_Outer4_AllBands', ...
    'S3_Inner4_AllBands', 'S3_Outer10_AllBands', 'S3_Inner10_AllBands', ...
    'Full_CARS', 'S3_Central_CARS', 'S3_Outer4_CARS', 'S3_Inner4_CARS', ...
    'S3_Outer10_CARS', 'S3_Inner10_CARS'};
numRanges = numel(rangeMasks);

cfg.outerRepeats = 5;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.carsRuns = 200;
cfg.svrC = [3 10 30 100];
cfg.svrEpsilon = [0.3 0.4 0.5];
cfg.svrScale = [3 5 7 10];
cfg.bootstrapReplicates = 2000;

n = numel(Y);
predCube = nan(n, numel(pipelineNames), cfg.outerRepeats);
foldRows = cell(0, 13);
metricRows = cell(0, 8);
pairRows = cell(0, 17);

for rep = 1:cfg.outerRepeats
    fprintf('\nOuter repeat %d/%d\n', rep, cfg.outerRepeats);
    outer = makeFoldsLocal(n, cfg.outerFolds, 1, 20260712 + 1000 * rep);
    bandCounts = nan(cfg.outerFolds, numel(pipelineNames));
    for foldId = 1:cfg.outerFolds
        tr = outer{foldId}.train;
        te = outer{foldId}.test;
        fprintf('  Fold %d/%d: train=%d test=%d\n', foldId, cfg.outerFolds, sum(tr), sum(te));
        for rangeId = 1:numRanges
            Xrange = rangeSpectra{rangeId};
            for featureId = 1:2
                if featureId == 1
                    selected = 1:size(Xrange, 2);
                    selectionCV = NaN;
                else
                    [selected, selectionCV] = selectCARSLocal( ...
                        Xrange(tr, :), Y(tr), cfg, 20260712 + 100000 * rep + 1000 * foldId);
                end
                pipelineId = pipelineIndex(rangeId, featureId, numRanges);
                [prediction, C, epsilon, scale, innerCV] = fitSVRLocal( ...
                    Xrange(tr, selected), Y(tr), Xrange(te, selected), cfg, ...
                    20260712 + 200000 * rep + 2000 * foldId + featureId);
                predCube(te, pipelineId, rep) = prediction;
                bandCounts(foldId, pipelineId) = numel(selected);
                foldRows(end+1, :) = {rep, foldId, pipelineNames{pipelineId}, ...
                    rangeNames{rangeId}, featureNames{featureId}, sum(tr), sum(te), ...
                    numel(selected), selectionCV, C, epsilon, scale, innerCV}; %#ok<AGROW>
            end
        end
    end

    for pipelineId = 1:numel(pipelineNames)
        prediction = predCube(:, pipelineId, rep);
        if any(~isfinite(prediction))
            error('Missing outer-fold predictions for %s, repeat %d.', pipelineNames{pipelineId}, rep);
        end
        met = predictionMetricsLocal(Y, prediction);
        metricRows(end+1, :) = {rep, pipelineNames{pipelineId}, ...
            mean(bandCounts(:, pipelineId)), std(bandCounts(:, pipelineId)), ...
            met.R2, met.RMSE, met.RPD, corr(Y, prediction)}; %#ok<AGROW>
    end

    for featureId = 1:2
        fullId = pipelineIndex(1, featureId, numRanges);
        fullPrediction = predCube(:, fullId, rep);
        fullMet = predictionMetricsLocal(Y, fullPrediction);
        for rangeId = 2:numRanges
            segmentId = pipelineIndex(rangeId, featureId, numRanges);
            segmentPrediction = predCube(:, segmentId, rep);
            segmentMet = predictionMetricsLocal(Y, segmentPrediction);
            [deltaR2, r2Low, r2High, deltaRMSE, rmseLow, rmseHigh] = ...
                bootstrapDeltaLocal(Y, segmentPrediction, fullPrediction, ...
                cfg.bootstrapReplicates, 20260712 + 3000 * rep + 100 * featureId + rangeId);
            pairRows(end+1, :) = {rep, rangeNames{rangeId}, featureNames{featureId}, ...
                fullMet.R2, segmentMet.R2, deltaR2, r2Low, r2High, ...
                fullMet.RMSE, segmentMet.RMSE, deltaRMSE, rmseLow, rmseHigh, ...
                deltaR2 > 0, deltaRMSE < 0, r2Low > 0, rmseHigh < 0}; %#ok<AGROW>
        end
    end
end

foldTable = cell2table(foldRows, 'VariableNames', {'Repeat', 'OuterFold', 'Pipeline', ...
    'RangeName', 'FeatureMethod', 'NTrain', 'NTest', 'BandCount', ...
    'SelectionCVRMSE', 'C', 'Epsilon', 'KernelScale', 'InnerCVRMSE'});
metricTable = cell2table(metricRows, 'VariableNames', {'Repeat', 'Pipeline', ...
    'MeanBandCount', 'SDBandCount', 'R2OOF', 'RMSEOOF', 'RPDOOF', 'Correlation'});
pairTable = cell2table(pairRows, 'VariableNames', {'Repeat', 'SegmentRange', 'FeatureMethod', ...
    'FullR2OOF', 'SegmentR2OOF', 'DeltaR2OOF', 'DeltaR2CI95Low', 'DeltaR2CI95High', ...
    'FullRMSEOOF', 'SegmentRMSEOOF', 'DeltaRMSEOOF', 'DeltaRMSECI95Low', ...
    'DeltaRMSECI95High', 'R2PointImproved', 'RMSEPointImproved', ...
    'R2CISupportsImprovement', 'RMSECISupportsImprovement'});

summaryRows = cell(0, 10);
for pipelineId = 1:numel(pipelineNames)
    q = metricTable(strcmp(metricTable.Pipeline, pipelineNames{pipelineId}), :);
    summaryRows(end+1, :) = {pipelineNames{pipelineId}, mean(q.MeanBandCount), ...
        std(q.MeanBandCount), mean(q.R2OOF), std(q.R2OOF), min(q.R2OOF), max(q.R2OOF), ...
        mean(q.RMSEOOF), std(q.RMSEOOF), mean(q.RPDOOF)}; %#ok<AGROW>
end
summaryTable = cell2table(summaryRows, 'VariableNames', {'Pipeline', 'MeanBandCount', ...
    'SDBandCount', 'MeanR2OOF', 'SDR2OOF', 'MinR2OOF', 'MaxR2OOF', ...
    'MeanRMSEOOF', 'SDRMSEOOF', 'MeanRPDOOF'});

predictionTable = table;
for rep = 1:cfg.outerRepeats
    predictionNames = [{'Repeat', 'OriginalSampleIndex', 'ActualSSC'}, ...
        strcat(pipelineNames, '_Prediction')];
    block = array2table([repmat(rep, n, 1), originalIndex, Y, predCube(:, :, rep)], ...
        'VariableNames', matlab.lang.makeValidName(predictionNames));
    predictionTable = [predictionTable; block]; %#ok<AGROW>
end

outXlsx = fullfile(outDir, ['S3边界重复嵌套交叉验证结果' suffix '.xlsx']);
if exist(outXlsx, 'file')
    delete(outXlsx);
end
writetable(summaryTable, outXlsx, 'Sheet', '汇总');
writetable(metricTable, outXlsx, 'Sheet', '重复验证指标');
writetable(pairTable, outXlsx, 'Sheet', '配对比较');
writetable(foldTable, outXlsx, 'Sheet', '外层折明细');
writetable(predictionTable, outXlsx, 'Sheet', '折外预测');

save(fullfile(outDir, ['S3边界重复嵌套交叉验证结果' suffix '.mat']), 'cfg', 'wl', 'Y', ...
    'originalIndex', 'predCube', 'summaryTable', 'metricTable', 'pairTable', 'foldTable', '-v7.3');
summaryFile = fullfile(outDir, ['S3边界重复嵌套交叉验证' suffix '_摘要.txt']);
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Repeated nested record-level cross-validation\n');
fprintf(fid, 'Data subset: %s\n', subset);
fprintf(fid, 'Outer design: %d repeats x %d folds; inner tuning: %d folds\n', ...
    cfg.outerRepeats, cfg.outerFolds, cfg.innerFolds);
fprintf(fid, 'All feature selection and SVR tuning occurred inside each outer training fold.\n\n');
for i = 1:height(summaryTable)
    fprintf(fid, '%s | bands %.2f +/- %.2f | R2 %.4f +/- %.4f | RMSE %.4f +/- %.4f | RPD %.4f\n', ...
        summaryTable.Pipeline{i}, summaryTable.MeanBandCount(i), summaryTable.SDBandCount(i), ...
        summaryTable.MeanR2OOF(i), summaryTable.SDR2OOF(i), ...
        summaryTable.MeanRMSEOOF(i), summaryTable.SDRMSEOOF(i), summaryTable.MeanRPDOOF(i));
end
fprintf(fid, '\nPoint-improvement counts across repeats:\n');
for rangeId = 2:numRanges
    for featureId = 1:2
        q = pairTable(strcmp(pairTable.SegmentRange, rangeNames{rangeId}) & ...
            strcmp(pairTable.FeatureMethod, featureNames{featureId}), :);
        fprintf(fid, '%s | %s: R2 %d/%d; RMSE %d/%d; both-CI-supported %d/%d\n', ...
            rangeNames{rangeId}, featureNames{featureId}, sum(q.R2PointImproved), height(q), ...
            sum(q.RMSEPointImproved), height(q), ...
            sum(q.R2CISupportsImprovement & q.RMSECISupportsImprovement), height(q));
    end
end
fclose(fid);

disp(summaryTable);
disp(pairTable);
fprintf('Nested validation outputs: %s\n', outDir);
fprintf('Nested validation completed: %s\n', datestr(now, 31));
end

function index = pipelineIndex(rangeId, featureId, numRanges)
index = (featureId - 1) * numRanges + rangeId;
end

function [selected, cvRmse] = selectCARSLocal(X, Y, cfg, seed)
rng(seed, 'twister');
maxLV = min([15, size(X, 2), size(X, 1) - 2]);
evalc('cars = carspls(X, Y, maxLV, cfg.innerFolds, ''center'', cfg.carsRuns, 1, 0, 1);');
selected = unique(cars.vsel(:)', 'stable');
if isempty(selected)
    error('CARS returned an empty variable subset.');
end
cvRmse = cars.RMSECV_min;
end

function [prediction, bestC, bestEpsilon, bestScale, bestRMSE] = fitSVRLocal(Xtr, Ytr, Xte, cfg, seed)
folds = makeFoldsLocal(numel(Ytr), cfg.innerFolds, 1, seed);
bestRMSE = inf;
bestC = NaN;
bestEpsilon = NaN;
bestScale = NaN;
for C = cfg.svrC
    for epsilon = cfg.svrEpsilon
        for scale = cfg.svrScale
            value = svrCVRMSELocal(Xtr, Ytr, C, epsilon, scale, folds);
            if value < bestRMSE
                bestRMSE = value;
                bestC = C;
                bestEpsilon = epsilon;
                bestScale = scale;
            end
        end
    end
end
model = fitrsvm(Xtr, Ytr, 'KernelFunction', 'gaussian', 'BoxConstraint', bestC, ...
    'Epsilon', bestEpsilon, 'KernelScale', bestScale, 'Standardize', true);
prediction = predict(model, Xte);
end

function rmse = svrCVRMSELocal(X, Y, C, epsilon, scale, folds)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    model = fitrsvm(X(fold.train, :), Y(fold.train), 'KernelFunction', 'gaussian', ...
        'BoxConstraint', C, 'Epsilon', epsilon, 'KernelScale', scale, 'Standardize', true);
    prediction = predict(model, X(fold.test, :));
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / count);
end

function folds = makeFoldsLocal(n, k, repeats, seed)
folds = cell(k * repeats, 1);
counter = 0;
for rep = 1:repeats
    rng(seed + rep - 1, 'twister');
    cv = cvpartition(n, 'KFold', k);
    for foldId = 1:k
        counter = counter + 1;
        folds{counter} = struct('train', training(cv, foldId), 'test', test(cv, foldId));
    end
end
end

function output = snvLocal(X)
rowMean = mean(X, 2);
rowStd = std(X, 0, 2);
rowStd(rowStd < eps) = 1;
output = (X - rowMean) ./ rowStd;
end

function output = firstDerivativeLocal(X, window)
[~, g] = sgolay(2, window);
half = floor(window / 2);
output = zeros(size(X));
for i = 1:size(X, 1)
    padded = [X(i, half:-1:1), X(i, :), X(i, end:-1:end-half+1)];
    derivative = conv(padded, g(:, 2), 'same');
    output(i, :) = derivative(half+1:end-half);
end
end

function met = predictionMetricsLocal(y, prediction)
met.R2 = 1 - sum((y - prediction) .^ 2) / sum((y - mean(y)) .^ 2);
met.RMSE = sqrt(mean((y - prediction) .^ 2));
met.RPD = std(y) / met.RMSE;
end

function [deltaR2, r2Low, r2High, deltaRMSE, rmseLow, rmseHigh] = bootstrapDeltaLocal(y, segmentPrediction, fullPrediction, repetitions, seed)
deltaR2 = r2ValueLocal(y, segmentPrediction) - r2ValueLocal(y, fullPrediction);
deltaRMSE = rmseValueLocal(y, segmentPrediction) - rmseValueLocal(y, fullPrediction);
rng(seed, 'twister');
n = numel(y);
r2Boot = nan(repetitions, 1);
rmseBoot = nan(repetitions, 1);
for i = 1:repetitions
    index = randi(n, n, 1);
    r2Boot(i) = r2ValueLocal(y(index), segmentPrediction(index)) - ...
        r2ValueLocal(y(index), fullPrediction(index));
    rmseBoot(i) = rmseValueLocal(y(index), segmentPrediction(index)) - ...
        rmseValueLocal(y(index), fullPrediction(index));
end
r2Bounds = prctile(r2Boot, [2.5 97.5]);
rmseBounds = prctile(rmseBoot, [2.5 97.5]);
r2Low = r2Bounds(1);
r2High = r2Bounds(2);
rmseLow = rmseBounds(1);
rmseHigh = rmseBounds(2);
end

function value = r2ValueLocal(y, prediction)
value = 1 - sum((y - prediction) .^ 2) / sum((y - mean(y)) .^ 2);
end

function value = rmseValueLocal(y, prediction)
value = sqrt(mean((y - prediction) .^ 2));
end
