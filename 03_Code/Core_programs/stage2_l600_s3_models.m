function stage2_l600_s3_models(mode)
% Fixed-split SSC comparisons; workflow selection uses the same record-level dataset.
% These estimates do not establish independent-fruit generalization.

if nargin < 1
    mode = 'all4test';
end
mode = validatestring(mode, {'all4test','fourrefine','s3boundary','s3test','smoke'});
isSmoke = strcmpi(mode, 'smoke');
isS3Test = strcmpi(mode, 's3test');
isAll4Test = strcmpi(mode, 'all4test');
isFourRefine = strcmpi(mode, 'fourrefine');
isS3Boundary = strcmpi(mode, 's3boundary');
isDetailed = isFourRefine || isS3Boundary;
rng(20260622, 'twister');
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
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

cfg = defaultConfig(isSmoke);
if isDetailed
    % Hyperparameters use training CV; the shared prediction set also informed workflow selection.
    cfg.bootstrapReplicates = 2000;
    cfg.carsRuns = 200;
    cfg.randomFrogIterations = 800;
    cfg.randomFrogChains = 5;
    cfg.cvRepeats = 2;
    cfg.finalCVRepeats = 10;
    cfg.svrC = [1 3 10 30 100];
    cfg.svrEpsilon = [0.10 0.20 0.30 0.40 0.50];
    cfg.svrScale = [1 2 3 5 7 10];
end
if isAll4Test
    resultStem = '四波段完整模型结果';
elseif isS3Test
    resultStem = 'S3初步模型结果';
elseif isFourRefine
    resultStem = 'S3_S4精细复核模型结果';
elseif isS3Boundary
    resultStem = 'S3边界敏感性模型结果';
elseif isSmoke
    resultStem = '冒烟测试模型结果';
else
    resultStem = '完整模型矩阵结果';
end
logFile = fullfile(outDir, [resultStem '_运行日志.txt']);
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);
cleanupObj = onCleanup(@() diary('off'));

fprintf('Stage 2 modelling started: %s\n', datestr(now, 31));
fprintf('Mode: %s\n', mode);
fprintf('Input: %s\n', inputFile);

[wl, Y, Xraw, ~] = load_l600_submission_ssc(inputFile, 'SSC_model_spectra');
if numel(wl) ~= size(Xraw, 2) || numel(Y) ~= size(Xraw, 1)
    error('Unexpected model data dimensions.');
end

tr = spxy(Xraw, Y, round(0.80 * size(Xraw, 1)));
te = setdiff(1:size(Xraw, 1), tr);
fprintf('Fixed SPXY split: training=%d, test=%d.\n', numel(tr), numel(te));

if isAll4Test
    rangeDefs = buildAll4Ranges(wl);
elseif isS3Test || isSmoke
    rangeDefs = buildS3TestRanges(wl);
elseif isFourRefine
    rangeDefs = buildFourRefineRanges(wl);
elseif isS3Boundary
    rangeDefs = buildS3BoundarySensitivityRanges(wl);
end

featureNames = {'AllBands', 'CARS', 'SPA', 'RandomFrog'};
modelNames = {'PLSR', 'Ridge', 'SVR', 'RandomForest'};
featureRunIndices = 1:numel(featureNames);
useSnvOnly = isDetailed;
if useSnvOnly
    modelNames = {'SVR'};
end
if isS3Boundary
    % All-band SVR tests the continuous interval itself; CARS is retained as
    % a secondary check of interval-by-feature-selection interaction.
    featureRunIndices = 1:2;
end

results = repmat(emptyResult(), 0, 1);
bandRows = cell(0, 8);
selectionRows = cell(0, 8);
runCount = 0;

for r = 1:height(rangeDefs)
    rangeMask = rangeDefs.StartIndex(r):rangeDefs.EndIndex(r);
    % Crop first so every segmented model can be reproduced from its own
    % wavelength interval without using out-of-range values in preprocessing.
    rangeRaw = Xraw(:, rangeMask);
    preps = preprocessData(rangeRaw, tr, te);
    if useSnvOnly
        preps = preps(2); % SNV_D1_7
    end
    fprintf('\nRange %d/%d: %s (%g-%g nm, %d variables)\n', ...
        r, height(rangeDefs), rangeDefs.RangeName{r}, rangeDefs.StartNm(r), ...
        rangeDefs.EndNm(r), numel(rangeMask));
    for p = 1:numel(preps)
        Xtr = preps(p).Xtr;
        Xte = preps(p).Xte;
        fprintf('  Preprocess: %s\n', preps(p).Name);
        selections = selectFeatureSets(Xtr, Y(tr), featureNames, cfg, ...
            20260622 + 1000 * r + 100 * p);

        for f = featureRunIndices
            sel = selections(f).Indices;
            if isempty(sel)
                sel = 1:size(Xtr, 2);
            end
            selectedWavelengths = wl(rangeMask(sel));
            for b = 1:numel(sel)
                bandRows(end+1, :) = {rangeDefs.RangeId(r), rangeDefs.RangeName{r}, ...
                    preps(p).Name, featureNames{f}, b, sel(b), ...
                    rangeMask(sel(b)), selectedWavelengths(b)}; %#ok<AGROW>
            end
            selectionRows(end+1, :) = {rangeDefs.RangeId(r), rangeDefs.RangeName{r}, ...
                preps(p).Name, featureNames{f}, numel(sel), selections(f).CVRMSE, ...
                selections(f).Detail, selections(f).Seed}; %#ok<AGROW>

            for m = 1:numel(modelNames)
                runCount = runCount + 1;
                fprintf('    %03d: %s + %s + %s\n', runCount, preps(p).Name, featureNames{f}, modelNames{m});
                try
                    [ptr, pte, param, cvRmse, numericParams] = fitModel( ...
                        modelNames{m}, Xtr(:, sel), Y(tr), Xte(:, sel), cfg, ...
                        20260622 + 10000 * r + 1000 * p + 100 * f + m);
                    met = metrics(Y(tr), ptr, Y(te), pte);
                    errorText = '';
                catch ME
                    ptr = nan(numel(tr), 1);
                    pte = nan(numel(te), 1);
                    met = struct('R2Train', NaN, 'R2Test', NaN, 'RMSETrain', NaN, ...
                        'RMSETest', NaN, 'RPDTrain', NaN, 'RPDTest', NaN, 'R2Gap', NaN);
                    param = '';
                    cvRmse = NaN;
                    numericParams = [NaN NaN NaN NaN];
                    errorText = ME.message;
                    warning('Model failed: %s', ME.message);
                end
                result = emptyResult();
                result.Key = sprintf('R%02d_P%02d_F%02d_M%02d', r, p, f, m);
                result.RangeId = rangeDefs.RangeId(r);
                result.RangeName = rangeDefs.RangeName{r};
                result.StartNm = rangeDefs.StartNm(r);
                result.EndNm = rangeDefs.EndNm(r);
                result.Preprocess = preps(p).Name;
                result.FeatureMethod = featureNames{f};
                result.Model = modelNames{m};
                result.BandCount = numel(sel);
                result.NTrain = numel(tr);
                result.NTest = numel(te);
                result.R2Train = met.R2Train;
                result.R2Test = met.R2Test;
                result.RMSETrain = met.RMSETrain;
                result.RMSETest = met.RMSETest;
                result.RPDTrain = met.RPDTrain;
                result.RPDTest = met.RPDTest;
                result.R2Gap = met.R2Gap;
                result.Stable = isfinite(met.R2Gap) && met.R2Gap <= 0.10;
                result.CVRMSE = cvRmse;
                result.Param1 = numericParams(1);
                result.Param2 = numericParams(2);
                result.Param3 = numericParams(3);
                result.Param4 = numericParams(4);
                result.Parameters = param;
                result.Error = errorText;
                result.TrainIndex = tr;
                result.TestIndex = te;
                result.ActualTrain = Y(tr);
                result.PredictedTrain = ptr;
                result.ActualTest = Y(te);
                result.PredictedTest = pte;
                result.SelectedLocalIndex = sel(:)';
                result.SelectedGlobalIndex = rangeMask(sel(:))';
                result.SelectedWavelengthNm = selectedWavelengths(:)';
                results(end+1) = result; %#ok<AGROW>
            end
        end
    end
end

metricsTable = resultsToTable(results);
metricsTable.Score = metricsTable.R2Test - 0.20 .* metricsTable.R2Gap - 0.01 .* metricsTable.RMSETest;
if isDetailed
    % Rank refined candidates by repeated training-set CV, never by test performance.
    metricsTable = sortrows(metricsTable, {'CVRMSE', 'BandCount'}, {'ascend', 'ascend'});
else
    metricsTable = sortrows(metricsTable, {'Stable', 'R2Test', 'RMSETest'}, {'descend', 'descend', 'ascend'});
end
predictionTable = predictionsToTable(results);
bandTable = cell2table(bandRows, 'VariableNames', {'RangeId', 'RangeName', 'Preprocess', ...
    'FeatureMethod', 'ImportanceRank', 'LocalVariableIndex', 'GlobalVariableIndex', 'WavelengthNm'});
selectionTable = cell2table(selectionRows, 'VariableNames', {'RangeId', 'RangeName', ...
    'Preprocess', 'FeatureMethod', 'BandCount', 'SelectionCVRMSE', 'SelectionDetail', 'Seed'});
pairTable = pairedComparisons(results, cfg.bootstrapReplicates);
splitTable = table((1:numel(Y))', ismember((1:numel(Y))', tr), ismember((1:numel(Y))', te), Y, ...
    'VariableNames', {'OriginalSampleIndex', 'InTrainingSet', 'InTestSet', 'ActualSSC'});
rangeTable = rangeDefs;

outXlsx = fullfile(outDir, [resultStem '.xlsx']);
if exist(outXlsx, 'file')
    delete(outXlsx);
end
writetable(rangeTable, outXlsx, 'Sheet', '光谱范围');
writetable(splitTable, outXlsx, 'Sheet', 'SPXY划分');
writetable(metricsTable, outXlsx, 'Sheet', '全部模型指标');
writetable(metricsTable(metricsTable.Stable, :), outXlsx, 'Sheet', '稳定模型');
writetable(selectionTable, outXlsx, 'Sheet', '变量筛选汇总');
writetable(bandTable, outXlsx, 'Sheet', '筛选波长');
writetable(pairTable, outXlsx, 'Sheet', '相对全波段配对比较');
writetable(predictionTable, outXlsx, 'Sheet', '预测明细');

summaryFile = fullfile(outDir, [resultStem '_摘要.txt']);
writeSummary(summaryFile, mode, cfg, rangeTable, metricsTable, pairTable, tr, te);
save(fullfile(outDir, [resultStem '.mat']), ...
    'cfg', 'rangeDefs', 'tr', 'te', 'results', 'metricsTable', 'selectionTable', ...
    'bandTable', 'pairTable', 'wl', 'Y', '-v7.3');
fprintf('\nCompleted %d model fits.\n', numel(results));
fprintf('Stable models: %d/%d.\n', sum(metricsTable.Stable), height(metricsTable));
if any(metricsTable.Stable)
    if isDetailed
        top = metricsTable(1, :);
        fprintf('Lowest training-CV model: %s | %s | %s | %s | CVRMSE=%.4f | Rp2=%.4f | RMSEP=%.4f\n', ...
            top.RangeName{1}, top.Preprocess{1}, top.FeatureMethod{1}, top.Model{1}, ...
            top.CVRMSE(1), top.R2Test(1), top.RMSETest(1));
    else
        top = metricsTable(find(metricsTable.Stable, 1, 'first'), :);
        fprintf('Top stable model: %s | %s | %s | %s | Rp2=%.4f | RMSEP=%.4f | Gap=%.4f\n', ...
            top.RangeName{1}, top.Preprocess{1}, top.FeatureMethod{1}, top.Model{1}, ...
            top.R2Test(1), top.RMSETest(1), top.R2Gap(1));
    end
end
fprintf('Outputs: %s\n', outDir);
fprintf('Stage 2 modelling completed: %s\n', datestr(now, 31));
end

function cfg = defaultConfig(isSmoke)
cfg.bootstrapReplicates = 1000;
cfg.carsRuns = 80;
cfg.randomFrogIterations = 300;
cfg.randomFrogChains = 2;
cfg.spaMaxVariables = 30;
cfg.cvFolds = 5;
cfg.cvRepeats = 1;
cfg.finalCVRepeats = 3;
cfg.ridgeLambdas = logspace(-5, 5, 25);
cfg.svrC = [1 10 100];
cfg.svrEpsilon = [0.1 0.3 0.5];
cfg.svrScale = [1 3 10];
cfg.rfLeaf = [1 5 10];
cfg.rfTreesTune = 100;
cfg.rfTreesFinal = 500;
cfg.rfTuneRepeats = 1;
if isSmoke
    cfg.bootstrapReplicates = 100;
    cfg.carsRuns = 8;
    cfg.randomFrogIterations = 15;
    cfg.randomFrogChains = 1;
    cfg.cvRepeats = 1;
    cfg.finalCVRepeats = 1;
    cfg.ridgeLambdas = logspace(-3, 3, 7);
    cfg.svrC = [1 10];
    cfg.svrEpsilon = [0.1 0.3];
    cfg.svrScale = [1 3];
    cfg.rfLeaf = [1 5];
    cfg.rfTreesTune = 30;
    cfg.rfTreesFinal = 80;
end
end

function rangeDefs = buildS3TestRanges(wl)
labels = {'Full_650_950'; 'S3_748_826'};
starts = [650; 748];
ends = [950; 826];
startIndex = zeros(numel(starts), 1);
endIndex = zeros(numel(starts), 1);
for i = 1:numel(starts)
    startIndex(i) = find(wl == starts(i), 1, 'first');
    endIndex(i) = find(wl == ends(i), 1, 'first');
    if isempty(startIndex(i)) || isempty(endIndex(i))
        error('S3 wavelength boundary is absent from the 2 nm grid.');
    end
end
rangeDefs = table((0:numel(starts)-1)', labels, startIndex, endIndex, starts, ends, ...
    endIndex - startIndex + 1, 'VariableNames', {'RangeId', 'RangeName', ...
    'StartIndex', 'EndIndex', 'StartNm', 'EndNm', 'VariableCount'});
end

function rangeDefs = buildAll4Ranges(wl)
labels = {'Full_650_950'; 'S1_650_688'; 'S2_690_746'; ...
    'S3_748_826'; 'S4_828_950'};
starts = [650; 650; 690; 748; 828];
ends = [950; 688; 746; 826; 950];
startIndex = zeros(numel(starts), 1);
endIndex = zeros(numel(starts), 1);
for i = 1:numel(starts)
    startIndex(i) = find(wl == starts(i), 1, 'first');
    endIndex(i) = find(wl == ends(i), 1, 'first');
    if isempty(startIndex(i)) || isempty(endIndex(i))
        error('Four-band wavelength boundary is absent from the 2 nm grid.');
    end
end
rangeDefs = table((0:numel(starts)-1)', labels, startIndex, endIndex, starts, ends, ...
    endIndex - startIndex + 1, 'VariableNames', {'RangeId', 'RangeName', ...
    'StartIndex', 'EndIndex', 'StartNm', 'EndNm', 'VariableCount'});
end

function rangeDefs = buildFourRefineRanges(wl)
labels = {'Full_650_950'; 'S3_748_826'; 'S4_828_950'};
starts = [650; 748; 828];
ends = [950; 826; 950];
rangeDefs = makeRangeTable(wl, labels, starts, ends, 'Refinement');
end

function rangeDefs = buildS3BoundarySensitivityRanges(wl)
% Test whether the S3 conclusion depends on a single pair of boundaries.
labels = {'Full_650_950'; 'S3_Central_748_826'; 'S3_Outer4_744_830'; ...
    'S3_Inner4_752_822'; 'S3_Outer10_738_836'; 'S3_Inner10_758_816'};
starts = [650; 748; 744; 752; 738; 758];
ends = [950; 826; 830; 822; 836; 816];
rangeDefs = makeRangeTable(wl, labels, starts, ends, 'S3 boundary-sensitivity');
end

function rangeDefs = makeRangeTable(wl, labels, starts, ends, context)
startIndex = zeros(numel(starts), 1);
endIndex = zeros(numel(starts), 1);
for i = 1:numel(starts)
    startIndex(i) = find(wl == starts(i), 1, 'first');
    endIndex(i) = find(wl == ends(i), 1, 'first');
    if isempty(startIndex(i)) || isempty(endIndex(i))
        error('%s wavelength boundary is absent from the 2 nm grid.', context);
    end
end
rangeDefs = table((0:numel(starts)-1)', labels, startIndex, endIndex, starts, ends, ...
    endIndex - startIndex + 1, 'VariableNames', {'RangeId', 'RangeName', ...
    'StartIndex', 'EndIndex', 'StartNm', 'EndNm', 'VariableCount'});
end

function preps = preprocessData(Xraw, tr, te)
preps = repmat(struct('Name', '', 'Xtr', [], 'Xte', []), 3, 1);
XtrRaw = Xraw(tr, :);
XteRaw = Xraw(te, :);
XtrAir = airPLS(XtrRaw, 1e8, 2, 0.1, 0.05, 20);
XteAir = airPLS(XteRaw, 1e8, 2, 0.1, 0.05, 20);

refAir = mean(XtrAir, 1);
preps(1).Name = 'airPLS_MSC_D1_7';
preps(1).Xtr = firstDerivative(mscApply(XtrAir, refAir), 7);
preps(1).Xte = firstDerivative(mscApply(XteAir, refAir), 7);

preps(2).Name = 'SNV_D1_7';
preps(2).Xtr = firstDerivative(snv(XtrRaw), 7);
preps(2).Xte = firstDerivative(snv(XteRaw), 7);

refRaw = mean(XtrRaw, 1);
preps(3).Name = 'MSC_SG15';
preps(3).Xtr = sgolayfilt(mscApply(XtrRaw, refRaw)', 2, 15)';
preps(3).Xte = sgolayfilt(mscApply(XteRaw, refRaw)', 2, 15)';
end

function selections = selectFeatureSets(X, Y, names, cfg, seed)
p = size(X, 2);
selections = repmat(struct('Indices', [], 'CVRMSE', NaN, 'Detail', '', 'Seed', seed), numel(names), 1);
selections(1).Indices = 1:p;
selections(1).CVRMSE = plsCVRMSE(X, Y, min([20 p size(X, 1) - 2]), makeFolds(numel(Y), cfg.cvFolds, 1, seed));
selections(1).Detail = 'All variables within the defined spectral range';

try
    rng(seed + 11, 'twister');
    maxLV = min([15, p, size(X, 1) - 2]);
    evalc('cars = carspls(X, Y, maxLV, cfg.cvFolds, ''center'', cfg.carsRuns, 1, 0, 1);');
    selections(2).Indices = unique(cars.vsel, 'stable');
    selections(2).CVRMSE = cars.RMSECV_min;
    selections(2).Detail = sprintf('CARS runs=%d; optIteration=%d; optLV=%d', cfg.carsRuns, cars.iterOPT, cars.optLV);
    selections(2).Seed = seed + 11;
catch ME
    selections(2).Indices = 1:p;
    selections(2).Detail = ['CARS fallback: ' ME.message];
end

try
    [spaIndices, spaRmse] = fastSPA(X, Y, cfg.spaMaxVariables, seed + 22);
    selections(3).Indices = spaIndices;
    selections(3).CVRMSE = spaRmse;
    selections(3).Detail = sprintf('SPA projection paths; selected=%d', numel(spaIndices));
    selections(3).Seed = seed + 22;
catch ME
    selections(3).Indices = 1:p;
    selections(3).Detail = ['SPA fallback: ' ME.message];
end

try
    [frogIndices, frogRmse, frogDetail] = aggregatedRandomFrog(X, Y, cfg, seed + 33);
    selections(4).Indices = frogIndices;
    selections(4).CVRMSE = frogRmse;
    selections(4).Detail = frogDetail;
    selections(4).Seed = seed + 33;
catch ME
    selections(4).Indices = 1:p;
    selections(4).Detail = ['RandomFrog fallback: ' ME.message];
end
for i = 1:numel(selections)
    selections(i).Indices = unique(selections(i).Indices(:)', 'stable');
    selections(i).Indices = selections(i).Indices(selections(i).Indices >= 1 & selections(i).Indices <= p);
    if isempty(selections(i).Indices)
        selections(i).Indices = 1:p;
    end
end
end

function [indices, bestRmse] = fastSPA(X, Y, maxVariables, seed)
n = size(X, 1);
p = size(X, 2);
mmax = min([maxVariables, p, n - 2]);
if mmax < 2
    indices = 1:p;
    bestRmse = NaN;
    return;
end
mu = mean(X, 1);
sd = std(X, 0, 1);
sd(sd < eps) = 1;
Z = (X - mu) ./ sd;
candidateSizes = unique([1:3:mmax, mmax]);
folds = makeFolds(n, 5, 1, seed);
bestRmse = inf;
bestPath = [];
for startVar = 1:p
    path = projections_qr(Z, startVar, mmax);
    path = path(path >= 1 & path <= p);
    if numel(path) < 2
        continue;
    end
    for m = candidateSizes(candidateSizes <= numel(path))
        selected = path(1:m);
        rmse = linearCVRMSE(X(:, selected), Y, folds);
        if rmse < bestRmse
            bestRmse = rmse;
            bestPath = selected;
        end
    end
end
if isempty(bestPath)
    error('SPA did not produce a valid variable subset.');
end
indices = unique(bestPath, 'stable');
end

function [indices, bestRmse, detail] = aggregatedRandomFrog(X, Y, cfg, seed)
p = size(X, 2);
maxLV = min([15, p, size(X, 1) - 2]);
probability = zeros(1, p);
for chain = 1:cfg.randomFrogChains
    rng(seed + chain, 'twister');
    evalc('frog = randomfrog_pls(X, Y, maxLV, ''center'', cfg.randomFrogIterations, min(5, p), ''regcoef'');');
    probability = probability + frog.probability;
end
probability = probability ./ cfg.randomFrogChains;
[~, rank] = sort(probability, 'descend');
candidateCounts = unique(min(p, [5 10 15 20 25 30 40]));
candidateCounts = candidateCounts(candidateCounts >= 2);
folds = makeFolds(numel(Y), cfg.cvFolds, 1, seed + 100);
rmseValues = inf(size(candidateCounts));
for i = 1:numel(candidateCounts)
    rmseValues(i) = plsCVRMSE(X(:, rank(1:candidateCounts(i))), Y, ...
        min([20 candidateCounts(i) size(X, 1) - 2]), folds);
end
[bestRmse, idx] = min(rmseValues);
indices = sort(rank(1:candidateCounts(idx)));
detail = sprintf('RandomFrog chains=%d; iterations=%d; selected=%d', ...
    cfg.randomFrogChains, cfg.randomFrogIterations, numel(indices));
end

function [ptr, pte, parameterText, finalCVRMSE, numericParams] = fitModel(name, Xtr, Ytr, Xte, cfg, seed)
switch name
    case 'PLSR'
        [ptr, pte, ncomp] = fitPLSR(Xtr, Ytr, Xte, cfg, seed);
        finalCVRMSE = plsCVRMSE(Xtr, Ytr, ncomp, makeFolds(numel(Ytr), cfg.cvFolds, cfg.finalCVRepeats, seed + 500));
        parameterText = sprintf('PLSR components=%d', ncomp);
        numericParams = [ncomp NaN NaN NaN];
    case 'Ridge'
        [ptr, pte, lambda] = fitRidge(Xtr, Ytr, Xte, cfg, seed);
        finalCVRMSE = ridgeCVRMSE(Xtr, Ytr, lambda, makeFolds(numel(Ytr), cfg.cvFolds, cfg.finalCVRepeats, seed + 500));
        parameterText = sprintf('Ridge lambda=%.8g', lambda);
        numericParams = [lambda NaN NaN NaN];
    case 'SVR'
        [ptr, pte, C, epsilon, scale] = fitSVR(Xtr, Ytr, Xte, cfg, seed);
        finalCVRMSE = svrCVRMSE(Xtr, Ytr, C, epsilon, scale, makeFolds(numel(Ytr), cfg.cvFolds, cfg.finalCVRepeats, seed + 500));
        parameterText = sprintf('RBF-SVR C=%g epsilon=%g kernelScale=%g', C, epsilon, scale);
        numericParams = [C epsilon scale NaN];
    case 'RandomForest'
        [ptr, pte, leaf, mtry] = fitRandomForest(Xtr, Ytr, Xte, cfg, seed);
        finalCVRMSE = rfCVRMSE(Xtr, Ytr, leaf, mtry, cfg.rfTreesTune, makeFolds(numel(Ytr), cfg.cvFolds, cfg.finalCVRepeats, seed + 500), seed + 700);
        parameterText = sprintf('RandomForest trees=%d leaf=%d mtry=%d', cfg.rfTreesFinal, leaf, mtry);
        numericParams = [cfg.rfTreesFinal leaf mtry NaN];
    otherwise
        error('Unknown model: %s', name);
end
end

function [ptr, pte, bestComp] = fitPLSR(Xtr, Ytr, Xte, cfg, seed)
maxComp = min([20, size(Xtr, 2), size(Xtr, 1) - 2]);
folds = makeFolds(numel(Ytr), cfg.cvFolds, cfg.cvRepeats, seed);
rmse = inf(maxComp, 1);
for c = 1:maxComp
    rmse(c) = plsCVRMSE(Xtr, Ytr, c, folds);
end
[~, bestComp] = min(rmse);
[~, ~, ~, ~, beta] = plsregress(Xtr, Ytr, bestComp);
ptr = [ones(size(Xtr, 1), 1), Xtr] * beta;
pte = [ones(size(Xte, 1), 1), Xte] * beta;
end

function rmse = plsCVRMSE(X, Y, ncomp, folds)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    nUse = min([ncomp, size(X(fold.train, :), 2), sum(fold.train) - 2]);
    if nUse < 1
        rmse = inf;
        return;
    end
    [~, ~, ~, ~, beta] = plsregress(X(fold.train, :), Y(fold.train), nUse);
    prediction = [ones(sum(fold.test), 1), X(fold.test, :)] * beta;
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / max(count, 1));
end

function [ptr, pte, bestLambda] = fitRidge(Xtr, Ytr, Xte, cfg, seed)
folds = makeFolds(numel(Ytr), cfg.cvFolds, cfg.cvRepeats, seed);
values = inf(numel(cfg.ridgeLambdas), 1);
for i = 1:numel(cfg.ridgeLambdas)
    values(i) = ridgeCVRMSE(Xtr, Ytr, cfg.ridgeLambdas(i), folds);
end
[~, idx] = min(values);
bestLambda = cfg.ridgeLambdas(idx);
[mu, sd, Ztr] = standardizeTraining(Xtr);
Zte = (Xte - mu) ./ sd;
beta = ridgeBeta(Ztr, Ytr, bestLambda);
ptr = [ones(size(Ztr, 1), 1), Ztr] * beta;
pte = [ones(size(Zte, 1), 1), Zte] * beta;
end

function rmse = ridgeCVRMSE(X, Y, lambda, folds)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    [mu, sd, Ztr] = standardizeTraining(X(fold.train, :));
    Zva = (X(fold.test, :) - mu) ./ sd;
    beta = ridgeBeta(Ztr, Y(fold.train), lambda);
    prediction = [ones(size(Zva, 1), 1), Zva] * beta;
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / max(count, 1));
end

function beta = ridgeBeta(X, Y, lambda)
A = [ones(size(X, 1), 1), X];
penalty = eye(size(A, 2));
penalty(1, 1) = 0;
beta = (A' * A + lambda * penalty) \ (A' * Y);
end

function [ptr, pte, bestC, bestEpsilon, bestScale] = fitSVR(Xtr, Ytr, Xte, cfg, seed)
folds = makeFolds(numel(Ytr), cfg.cvFolds, cfg.cvRepeats, seed);
best = inf;
bestC = NaN;
bestEpsilon = NaN;
bestScale = NaN;
for C = cfg.svrC
    for epsilon = cfg.svrEpsilon
        for scale = cfg.svrScale
            value = svrCVRMSE(Xtr, Ytr, C, epsilon, scale, folds);
            if value < best
                best = value;
                bestC = C;
                bestEpsilon = epsilon;
                bestScale = scale;
            end
        end
    end
end
model = fitrsvm(Xtr, Ytr, 'KernelFunction', 'gaussian', 'BoxConstraint', bestC, ...
    'Epsilon', bestEpsilon, 'KernelScale', bestScale, 'Standardize', true);
ptr = predict(model, Xtr);
pte = predict(model, Xte);
end

function rmse = svrCVRMSE(X, Y, C, epsilon, scale, folds)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    try
        model = fitrsvm(X(fold.train, :), Y(fold.train), 'KernelFunction', 'gaussian', ...
            'BoxConstraint', C, 'Epsilon', epsilon, 'KernelScale', scale, 'Standardize', true);
        prediction = predict(model, X(fold.test, :));
    catch
        rmse = inf;
        return;
    end
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / max(count, 1));
end

function [ptr, pte, bestLeaf, bestMtry] = fitRandomForest(Xtr, Ytr, Xte, cfg, seed)
p = size(Xtr, 2);
mtryCandidates = unique(max(1, min(p, [round(sqrt(p)), round(p / 3)])));
folds = makeFolds(numel(Ytr), cfg.cvFolds, cfg.rfTuneRepeats, seed);
best = inf;
bestLeaf = NaN;
bestMtry = NaN;
for leaf = cfg.rfLeaf
    for mtry = mtryCandidates
        value = rfCVRMSE(Xtr, Ytr, leaf, mtry, cfg.rfTreesTune, folds, seed + 100 * leaf + mtry);
        if value < best
            best = value;
            bestLeaf = leaf;
            bestMtry = mtry;
        end
    end
end
rng(seed + 9000, 'twister');
model = TreeBagger(cfg.rfTreesFinal, Xtr, Ytr, 'Method', 'regression', ...
    'MinLeafSize', bestLeaf, 'NumPredictorsToSample', bestMtry);
ptr = predict(model, Xtr);
pte = predict(model, Xte);
end

function rmse = rfCVRMSE(X, Y, leaf, mtry, trees, folds, seed)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    rng(seed + i, 'twister');
    model = TreeBagger(trees, X(fold.train, :), Y(fold.train), 'Method', 'regression', ...
        'MinLeafSize', leaf, 'NumPredictorsToSample', mtry);
    prediction = predict(model, X(fold.test, :));
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / max(count, 1));
end

function folds = makeFolds(n, k, repeats, seed)
folds = cell(k * repeats, 1);
counter = 0;
for r = 1:repeats
    rng(seed + r - 1, 'twister');
    cv = cvpartition(n, 'KFold', k);
    for f = 1:k
        counter = counter + 1;
        folds{counter} = struct('train', training(cv, f), 'test', test(cv, f));
    end
end
end

function rmse = linearCVRMSE(X, Y, folds)
sumSq = 0;
count = 0;
for i = 1:numel(folds)
    fold = folds{i};
    A = [ones(sum(fold.train), 1), X(fold.train, :)];
    beta = pinv(A) * Y(fold.train);
    prediction = [ones(sum(fold.test), 1), X(fold.test, :)] * beta;
    err = Y(fold.test) - prediction;
    sumSq = sumSq + sum(err .^ 2);
    count = count + numel(err);
end
rmse = sqrt(sumSq / max(count, 1));
end

function [mu, sd, Z] = standardizeTraining(X)
mu = mean(X, 1);
sd = std(X, 0, 1);
sd(sd < eps) = 1;
Z = (X - mu) ./ sd;
end

function output = snv(X)
output = (X - mean(X, 2)) ./ max(std(X, 0, 2), eps);
end

function output = mscApply(X, reference)
output = zeros(size(X));
for i = 1:size(X, 1)
    coefficients = polyfit(reference, X(i, :), 1);
    output(i, :) = (X(i, :) - coefficients(2)) ./ coefficients(1);
end
end

function output = firstDerivative(X, window)
[~, g] = sgolay(2, window);
half = floor(window / 2);
output = zeros(size(X));
for i = 1:size(X, 1)
    padded = [X(i, half:-1:1), X(i, :), X(i, end:-1:end-half+1)];
    derivative = conv(padded, g(:, 2), 'same');
    output(i, :) = derivative(half+1:end-half);
end
end

function out = metrics(Ytr, Ptr, Yte, Pte)
out.R2Train = 1 - sum((Ytr - Ptr) .^ 2) / sum((Ytr - mean(Ytr)) .^ 2);
out.R2Test = 1 - sum((Yte - Pte) .^ 2) / sum((Yte - mean(Yte)) .^ 2);
out.RMSETrain = sqrt(mean((Ytr - Ptr) .^ 2));
out.RMSETest = sqrt(mean((Yte - Pte) .^ 2));
out.RPDTrain = std(Ytr) / out.RMSETrain;
out.RPDTest = std(Yte) / out.RMSETest;
out.R2Gap = abs(out.R2Train - out.R2Test);
end

function result = emptyResult()
result = struct('Key', '', 'RangeId', NaN, 'RangeName', '', 'StartNm', NaN, 'EndNm', NaN, ...
    'Preprocess', '', 'FeatureMethod', '', 'Model', '', 'BandCount', NaN, ...
    'NTrain', NaN, 'NTest', NaN, 'R2Train', NaN, 'R2Test', NaN, ...
    'RMSETrain', NaN, 'RMSETest', NaN, 'RPDTrain', NaN, 'RPDTest', NaN, ...
    'R2Gap', NaN, 'Stable', false, 'CVRMSE', NaN, 'Param1', NaN, ...
    'Param2', NaN, 'Param3', NaN, 'Param4', NaN, 'Parameters', '', 'Error', '', ...
    'TrainIndex', [], 'TestIndex', [], 'ActualTrain', [], 'PredictedTrain', [], ...
    'ActualTest', [], 'PredictedTest', [], 'SelectedLocalIndex', [], ...
    'SelectedGlobalIndex', [], 'SelectedWavelengthNm', []);
end

function tableOut = resultsToTable(results)
n = numel(results);
key = string({results.Key})';
rangeId = [results.RangeId]';
rangeName = string({results.RangeName})';
startNm = [results.StartNm]';
endNm = [results.EndNm]';
preprocess = string({results.Preprocess})';
featureMethod = string({results.FeatureMethod})';
model = string({results.Model})';
bandCount = [results.BandCount]';
nTrain = [results.NTrain]';
nTest = [results.NTest]';
r2Train = [results.R2Train]';
r2Test = [results.R2Test]';
rmseTrain = [results.RMSETrain]';
rmseTest = [results.RMSETest]';
rpdTrain = [results.RPDTrain]';
rpdTest = [results.RPDTest]';
r2Gap = [results.R2Gap]';
stable = [results.Stable]';
cvRmse = [results.CVRMSE]';
param1 = [results.Param1]';
param2 = [results.Param2]';
param3 = [results.Param3]';
param4 = [results.Param4]';
parameters = string({results.Parameters})';
errors = string({results.Error})';
tableOut = table(key, rangeId, rangeName, startNm, endNm, preprocess, featureMethod, model, ...
    bandCount, nTrain, nTest, r2Train, r2Test, rmseTrain, rmseTest, rpdTrain, rpdTest, ...
    r2Gap, stable, cvRmse, param1, param2, param3, param4, parameters, errors, ...
    'VariableNames', {'ModelKey', 'RangeId', 'RangeName', 'StartNm', 'EndNm', 'Preprocess', ...
    'FeatureMethod', 'Model', 'BandCount', 'NTrain', 'NTest', 'R2Train', 'R2Test', ...
    'RMSETrain', 'RMSETest', 'RPDTrain', 'RPDTest', 'R2Gap', 'Stable', 'CVRMSE', ...
    'Param1', 'Param2', 'Param3', 'Param4', 'Parameters', 'Error'});
if n == 0
    error('No model results were produced.');
end
end

function tableOut = predictionsToTable(results)
rows = cell(0, 9);
for i = 1:numel(results)
    r = results(i);
    for j = 1:numel(r.TrainIndex)
        rows(end+1, :) = {r.Key, r.RangeName, r.Preprocess, r.FeatureMethod, r.Model, ...
            'Train', r.TrainIndex(j), r.ActualTrain(j), r.PredictedTrain(j)}; %#ok<AGROW>
    end
    for j = 1:numel(r.TestIndex)
        rows(end+1, :) = {r.Key, r.RangeName, r.Preprocess, r.FeatureMethod, r.Model, ...
            'Test', r.TestIndex(j), r.ActualTest(j), r.PredictedTest(j)}; %#ok<AGROW>
    end
end
tableOut = cell2table(rows, 'VariableNames', {'ModelKey', 'RangeName', 'Preprocess', ...
    'FeatureMethod', 'Model', 'Set', 'OriginalSampleIndex', 'Actual', 'Predicted'});
end

function paired = pairedComparisons(results, repetitions)
rows = cell(0, 17);
for i = 1:numel(results)
    candidate = results(i);
    if candidate.RangeId == 0 || ~isfinite(candidate.R2Test)
        continue;
    end
    baseIndex = find([results.RangeId] == 0 & strcmp({results.Preprocess}, candidate.Preprocess) & ...
        strcmp({results.FeatureMethod}, candidate.FeatureMethod) & strcmp({results.Model}, candidate.Model), 1, 'first');
    if isempty(baseIndex) || ~isfinite(results(baseIndex).R2Test)
        continue;
    end
    base = results(baseIndex);
    [deltaR2, r2Low, r2High, deltaRMSE, rmseLow, rmseHigh] = bootstrapDelta( ...
        candidate.ActualTest, candidate.PredictedTest, base.PredictedTest, repetitions, i);
    improved = candidate.R2Test > base.R2Test && candidate.RMSETest < base.RMSETest;
    rows(end+1, :) = {candidate.RangeName, candidate.Preprocess, candidate.FeatureMethod, ...
        candidate.Model, base.BandCount, candidate.BandCount, base.R2Test, candidate.R2Test, ...
        deltaR2, r2Low, r2High, base.RMSETest, candidate.RMSETest, deltaRMSE, rmseLow, rmseHigh, improved}; %#ok<AGROW>
end
paired = cell2table(rows, 'VariableNames', {'RangeName', 'Preprocess', 'FeatureMethod', ...
    'Model', 'FullBandCount', 'SegmentBandCount', 'FullR2Test', 'SegmentR2Test', ...
    'DeltaR2Test', 'DeltaR2CI95Low', 'DeltaR2CI95High', 'FullRMSETest', ...
    'SegmentRMSETest', 'DeltaRMSETest', 'DeltaRMSECI95Low', 'DeltaRMSECI95High', 'ImprovedBoth'});
end

function [deltaR2, r2Low, r2High, deltaRMSE, rmseLow, rmseHigh] = bootstrapDelta(y, segmentPrediction, fullPrediction, repetitions, seed)
deltaR2 = r2Value(y, segmentPrediction) - r2Value(y, fullPrediction);
deltaRMSE = rmseValue(y, segmentPrediction) - rmseValue(y, fullPrediction);
rng(20270000 + seed, 'twister');
n = numel(y);
r2Values = zeros(repetitions, 1);
rmseValues = zeros(repetitions, 1);
for b = 1:repetitions
    idx = randi(n, n, 1);
    r2Values(b) = r2Value(y(idx), segmentPrediction(idx)) - r2Value(y(idx), fullPrediction(idx));
    rmseValues(b) = rmseValue(y(idx), segmentPrediction(idx)) - rmseValue(y(idx), fullPrediction(idx));
end
r2CI = prctile(r2Values, [2.5 97.5]);
rmseCI = prctile(rmseValues, [2.5 97.5]);
r2Low = r2CI(1);
r2High = r2CI(2);
rmseLow = rmseCI(1);
rmseHigh = rmseCI(2);
end

function value = r2Value(y, prediction)
value = 1 - sum((y - prediction).^2) / sum((y - mean(y)).^2);
end

function value = rmseValue(y, prediction)
value = sqrt(mean((y - prediction).^2));
end

function writeSummary(filePath, mode, cfg, rangeDefs, metricsTable, pairTable, tr, te)
fid = fopen(filePath, 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'Stage 2 L600 segmented modelling summary\n\n');
fprintf(fid, 'Mode: %s\n', mode);
fprintf(fid, 'SPXY split: training=%d; test=%d\n', numel(tr), numel(te));
fprintf(fid, 'Selectors: CARS=%d iterations; RandomFrog=%d chains x %d iterations; SPA max=%d variables\n', ...
    cfg.carsRuns, cfg.randomFrogChains, cfg.randomFrogIterations, cfg.spaMaxVariables);
fprintf(fid, 'Model CV: %d-fold screening; final CV repeats=%d\n\n', cfg.cvFolds, cfg.finalCVRepeats);
fprintf(fid, 'Ranges:\n');
for i = 1:height(rangeDefs)
    fprintf(fid, '  %s: %g-%g nm (%d variables)\n', rangeDefs.RangeName{i}, ...
        rangeDefs.StartNm(i), rangeDefs.EndNm(i), rangeDefs.VariableCount(i));
end
fprintf(fid, '\nStable model count: %d/%d\n', sum(metricsTable.Stable), height(metricsTable));
top = metricsTable(metricsTable.Stable, :);
if ~isempty(top)
    top = top(1:min(20, height(top)), :);
    fprintf(fid, '\nTop stable models:\n');
    for i = 1:height(top)
        fprintf(fid, '%s | %s | %s | %s | Rp2=%.6f | RMSEP=%.6f | Gap=%.6f | bands=%d\n', ...
            top.RangeName{i}, top.Preprocess{i}, top.FeatureMethod{i}, top.Model{i}, ...
            top.R2Test(i), top.RMSETest(i), top.R2Gap(i), top.BandCount(i));
    end
end
fprintf(fid, '\nPaired segment versus full comparisons: %d\n', height(pairTable));
end
