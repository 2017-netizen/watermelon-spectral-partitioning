function summarize_four_band_results()
% Consolidate fixed-test, nested-validation, and thickness-effect evidence.

rng(20260817, 'twister');
scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(scriptDir));
outDir = fullfile(packageRoot, '04_Model_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
fixed = load(fullfile(outDir, '四波段完整模型结果.mat'), 'metricsTable');
nested = load(fullfile(outDir, '四波段重复嵌套交叉验证结果.mat'), ...
    'Y', 'originalIndex', 'predCube', 'summaryTable', 'metricTable', 'pairTable');

rangeNames = {'Full_650_950', 'S1_650_688', 'S2_690_746', 'S3_748_826', 'S4_828_950'};
bestRows = cell(0, 11);
cvRows = cell(0, 9);
matchedRows = cell(0, 9);
for r = 1:numel(rangeNames)
    stable = fixed.metricsTable(strcmp(fixed.metricsTable.RangeName, rangeNames{r}) & ...
        fixed.metricsTable.Stable, :);
    byTest = sortrows(stable, {'R2Test', 'RMSETest'}, {'descend', 'ascend'});
    best = byTest(1, :);
    bestRows(end + 1, :) = modelRow(best); %#ok<AGROW>

    byCV = sortrows(stable, {'CVRMSE', 'BandCount'}, {'ascend', 'ascend'});
    chosen = byCV(1, :);
    cvRows(end + 1, :) = {rangeNames{r}, chosen.Preprocess{1}, ...
        chosen.FeatureMethod{1}, chosen.Model{1}, chosen.BandCount, chosen.CVRMSE, ...
        chosen.R2Test, chosen.RMSETest, chosen.RPDTest}; %#ok<AGROW>

    matched = fixed.metricsTable(strcmp(fixed.metricsTable.RangeName, rangeNames{r}) & ...
        strcmp(fixed.metricsTable.Preprocess, 'SNV_D1_7') & ...
        strcmp(fixed.metricsTable.FeatureMethod, 'AllBands') & ...
        strcmp(fixed.metricsTable.Model, 'SVR'), :);
    matchedRows(end + 1, :) = {rangeNames{r}, matched.BandCount, matched.R2Train, ...
        matched.R2Test, matched.RMSETrain, matched.RMSETest, matched.RPDTest, ...
        matched.R2Gap, matched.CVRMSE}; %#ok<AGROW>
end
bestTable = cell2table(bestRows, 'VariableNames', {'光谱范围', '预处理', ...
    '变量策略', '回归器', '变量数', 'Rc2', 'Rp2', 'RMSEC', 'RMSEP', 'RPDp', 'R2差值'});
cvChoiceTable = cell2table(cvRows, 'VariableNames', {'光谱范围', '预处理', ...
    '变量策略', '回归器', '变量数', '校正集CV_RMSE', 'Rp2', 'RMSEP', 'RPDp'});
matchedFixedTable = cell2table(matchedRows, 'VariableNames', {'光谱范围', ...
    '变量数', 'Rc2', 'Rp2', 'RMSEC', 'RMSEP', 'RPDp', 'R2差值', '校正集CV_RMSE'});

Y = nested.Y;
meanPrediction = mean(nested.predCube, 3);
pipelineNames = nested.summaryTable.Pipeline;
aggregateRows = cell(0, 5);
for p = 1:numel(pipelineNames)
    met = metrics(Y, meanPrediction(:, p));
    aggregateRows(end + 1, :) = {pipelineNames{p}, met.R2, met.RMSE, met.RPD, ...
        corr(Y, meanPrediction(:, p))}; %#ok<AGROW>
end
aggregateTable = cell2table(aggregateRows, 'VariableNames', ...
    {'流程', '汇总折外R2', '汇总折外RMSE', '汇总折外RPD', '相关系数'});

thicknessScore = [NaN; 0.141869; 0.797076; 0.124221; 0.740192];
allBandR2 = nested.summaryTable.MeanR2OOF(1:5);
allBandRMSE = nested.summaryTable.MeanRMSEOOF(1:5);
carsR2 = [nested.summaryTable.MeanR2OOF(6); nested.summaryTable.MeanR2OOF(7:10)];
carsRMSE = [nested.summaryTable.MeanRMSEOOF(6); nested.summaryTable.MeanRMSEOOF(7:10)];
evidenceTable = table(rangeNames', thicknessScore, allBandR2, allBandRMSE, carsR2, carsRMSE, ...
    'VariableNames', {'光谱范围', '平均标准化厚度效应', '全变量嵌套R2', ...
    '全变量嵌套RMSE', 'CARS嵌套R2', 'CARS嵌套RMSE'});

comparisonPairs = [2 1; 3 1; 4 1; 5 1; 7 6; 8 6; 9 6; 10 6; 4 5; 4 6; 5 6];
comparisonNames = {'S1全变量减全波段全变量'; 'S2全变量减全波段全变量'; ...
    'S3全变量减全波段全变量'; 'S4全变量减全波段全变量'; ...
    'S1CARS减全波段CARS'; 'S2CARS减全波段CARS'; 'S3CARS减全波段CARS'; ...
    'S4CARS减全波段CARS'; 'S3全变量减S4全变量'; ...
    'S3全变量减全波段CARS'; 'S4全变量减全波段CARS'};
bootstrapReplicates = 5000;
bootstrapR2 = nan(bootstrapReplicates, size(comparisonPairs, 1));
bootstrapRMSE = nan(bootstrapReplicates, size(comparisonPairs, 1));
comparisonRows = cell(0, 11);
for c = 1:size(comparisonPairs, 1)
    candidateId = comparisonPairs(c, 1);
    referenceId = comparisonPairs(c, 2);
    candidateMet = metrics(Y, meanPrediction(:, candidateId));
    referenceMet = metrics(Y, meanPrediction(:, referenceId));
    for b = 1:bootstrapReplicates
        index = randi(numel(Y), numel(Y), 1);
        candidateBoot = metrics(Y(index), meanPrediction(index, candidateId));
        referenceBoot = metrics(Y(index), meanPrediction(index, referenceId));
        bootstrapR2(b, c) = candidateBoot.R2 - referenceBoot.R2;
        bootstrapRMSE(b, c) = candidateBoot.RMSE - referenceBoot.RMSE;
    end
    r2CI = quantile(bootstrapR2(:, c), [0.025 0.975]);
    rmseCI = quantile(bootstrapRMSE(:, c), [0.025 0.975]);
    comparisonRows(end + 1, :) = {comparisonNames{c}, pipelineNames{candidateId}, ...
        pipelineNames{referenceId}, candidateMet.R2 - referenceMet.R2, ...
        r2CI(1), r2CI(2), candidateMet.RMSE - referenceMet.RMSE, ...
        rmseCI(1), rmseCI(2), mean(bootstrapR2(:, c) > 0), ...
        mean(bootstrapRMSE(:, c) < 0)}; %#ok<AGROW>
end
comparisonTable = cell2table(comparisonRows, 'VariableNames', ...
    {'比较', '候选流程', '参照流程', 'DeltaR2', 'DeltaR2_CI95下限', ...
    'DeltaR2_CI95上限', 'DeltaRMSE', 'DeltaRMSE_CI95下限', ...
    'DeltaRMSE_CI95上限', 'R2改善频率', 'RMSE改善频率'});

predictionNames = [{'原始记录索引', 'SSC实测值'}, strcat(pipelineNames', '_平均折外预测')];
predictionTable = array2table([nested.originalIndex, Y, meanPrediction], ...
    'VariableNames', matlab.lang.makeValidName(predictionNames));

outXlsx = fullfile(outDir, '四波段关键比较结果.xlsx');
if exist(outXlsx, 'file')
    delete(outXlsx);
end
writetable(bestTable, outXlsx, 'Sheet', '各范围固定预测最佳');
writetable(cvChoiceTable, outXlsx, 'Sheet', '校正集CV候选');
writetable(matchedFixedTable, outXlsx, 'Sheet', '相同流程固定预测');
writetable(nested.summaryTable, outXlsx, 'Sheet', '嵌套验证汇总');
writetable(nested.metricTable, outXlsx, 'Sheet', '各次重复指标');
writetable(aggregateTable, outXlsx, 'Sheet', '平均折外预测指标');
writetable(evidenceTable, outXlsx, 'Sheet', '厚度效应与建模证据');
writetable(comparisonTable, outXlsx, 'Sheet', '关键配对比较');
writetable(predictionTable, outXlsx, 'Sheet', '平均折外预测明细');

save(fullfile(outDir, '四波段关键比较结果.mat'), 'bestTable', 'cvChoiceTable', ...
    'matchedFixedTable', 'aggregateTable', 'evidenceTable', 'comparisonTable', ...
    'predictionTable', 'bootstrapR2', 'bootstrapRMSE', 'bootstrapReplicates');

fid = fopen(fullfile(outDir, '四波段关键比较结果摘要.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Full band and S1-S4 comparison\n\n');
fprintf(fid, 'Best stable fixed-test models:\n');
for i = 1:height(bestTable)
    fprintf(fid, '%s | %s | %s | %s | bands=%d | Rp2=%.4f | RMSEP=%.4f | RPD=%.4f\n', ...
        bestTable{i, 1}{1}, bestTable{i, 2}{1}, bestTable{i, 3}{1}, ...
        bestTable{i, 4}{1}, bestTable{i, 5}, bestTable{i, 7}, bestTable{i, 9}, bestTable{i, 10});
end
fprintf(fid, '\nRepeated nested validation:\n');
for i = 1:height(nested.summaryTable)
    fprintf(fid, '%s | bands %.2f | R2 %.4f +/- %.4f | RMSE %.4f +/- %.4f | RPD %.4f\n', ...
        nested.summaryTable.Pipeline{i}, nested.summaryTable.MeanBandCount(i), ...
        nested.summaryTable.MeanR2OOF(i), nested.summaryTable.SDR2OOF(i), ...
        nested.summaryTable.MeanRMSEOOF(i), nested.summaryTable.SDRMSEOOF(i), ...
        nested.summaryTable.MeanRPDOOF(i));
end
fprintf(fid, '\nKey aggregate comparisons:\n');
for i = 1:height(comparisonTable)
    fprintf(fid, '%s | DeltaR2 %.4f [%.4f, %.4f] | DeltaRMSE %.4f [%.4f, %.4f]\n', ...
        comparisonTable{i, 1}{1}, comparisonTable{i, 4}, comparisonTable{i, 5}, ...
        comparisonTable{i, 6}, comparisonTable{i, 7}, comparisonTable{i, 8}, ...
        comparisonTable{i, 9});
end
end

function row = modelRow(model)
row = {model.RangeName{1}, model.Preprocess{1}, model.FeatureMethod{1}, ...
    model.Model{1}, model.BandCount, model.R2Train, model.R2Test, ...
    model.RMSETrain, model.RMSETest, model.RPDTest, model.R2Gap};
end

function result = metrics(y, prediction)
residual = y - prediction;
result.R2 = 1 - sum(residual .^ 2) / sum((y - mean(y)) .^ 2);
result.RMSE = sqrt(mean(residual .^ 2));
result.RPD = std(y) / result.RMSE;
end
