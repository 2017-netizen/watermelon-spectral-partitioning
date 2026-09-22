function summarize_s3_boundary_results()
% Consolidate S3 boundary validation, refined models, and wavelength overlap.

rng(20260817, 'twister');
scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(scriptDir));
outDir = fullfile(packageRoot, '04_Model_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

nested = load(fullfile(outDir, 'S3边界重复嵌套交叉验证结果.mat'), ...
    'Y', 'originalIndex', 'predCube', 'summaryTable', 'metricTable', 'pairTable');
refined = load(fullfile(outDir, 'S3_S4精细复核模型结果.mat'), ...
    'metricsTable', 'pairTable');
fixed = load(fullfile(outDir, '四波段完整模型结果.mat'), 'bandTable');

Y = nested.Y;
meanPrediction = mean(nested.predCube, 3);
pipelineNames = nested.summaryTable.Pipeline;
aggregateRows = cell(numel(pipelineNames), 5);
for p = 1:numel(pipelineNames)
    met = metricsLocal(Y, meanPrediction(:, p));
    aggregateRows(p, :) = {pipelineNames{p}, met.R2, met.RMSE, met.RPD, ...
        corr(Y, meanPrediction(:, p))};
end
aggregateTable = cell2table(aggregateRows, 'VariableNames', ...
    {'流程', '汇总折外R2', '汇总折外RMSE', '汇总折外RPD', '相关系数'});

fullPairs = [2 1; 3 1; 4 1; 5 1; 6 1; 8 7; 9 7; 10 7; 11 7; 12 7];
centralPairs = [3 2; 4 2; 5 2; 6 2; 9 8; 10 8; 11 8; 12 8];
comparisonPairs = [fullPairs; centralPairs];
comparisonType = [repmat({'相对全波段'}, size(fullPairs, 1), 1); ...
    repmat({'相对S3中心区间'}, size(centralPairs, 1), 1)];
comparisonTable = bootstrapComparisons(Y, meanPrediction, pipelineNames, ...
    comparisonPairs, comparisonType, 5000);

[overlapTable, highFrequencyTable] = featureOverlap(fixed.bandTable);

predictionNames = [{'原始记录索引', 'SSC实测值'}, strcat(pipelineNames', '_平均折外预测')];
predictionTable = array2table([nested.originalIndex, Y, meanPrediction], ...
    'VariableNames', matlab.lang.makeValidName(predictionNames));

outXlsx = fullfile(outDir, 'S3边界与特征重合综合结果.xlsx');
if exist(outXlsx, 'file')
    delete(outXlsx);
end
writetable(nested.summaryTable, outXlsx, 'Sheet', '边界嵌套验证汇总');
writetable(nested.metricTable, outXlsx, 'Sheet', '各次重复指标');
writetable(aggregateTable, outXlsx, 'Sheet', '平均折外预测指标');
writetable(comparisonTable, outXlsx, 'Sheet', '边界配对比较');
writetable(refined.metricsTable, outXlsx, 'Sheet', 'S3_S4精细复核');
writetable(refined.pairTable, outXlsx, 'Sheet', '固定预测配对比较');
writetable(overlapTable, outXlsx, 'Sheet', '特征波长区间分布');
writetable(highFrequencyTable, outXlsx, 'Sheet', '重复选择波长');
writetable(predictionTable, outXlsx, 'Sheet', '平均折外预测明细');

save(fullfile(outDir, 'S3边界与特征重合综合结果.mat'), 'aggregateTable', ...
    'comparisonTable', 'overlapTable', 'highFrequencyTable', 'predictionTable', '-v7.3');

fid = fopen(fullfile(outDir, 'S3边界与特征重合综合结果摘要.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'S3 boundary validation and full-band feature overlap\n\n');
fprintf(fid, 'Repeated nested validation: 5 repeats x 5 outer folds; 5-fold inner tuning.\n');
for i = 1:height(nested.summaryTable)
    fprintf(fid, '%s | bands %.2f | R2 %.4f +/- %.4f | RMSE %.4f +/- %.4f | RPD %.4f\n', ...
        nested.summaryTable.Pipeline{i}, nested.summaryTable.MeanBandCount(i), ...
        nested.summaryTable.MeanR2OOF(i), nested.summaryTable.SDR2OOF(i), ...
        nested.summaryTable.MeanRMSEOOF(i), nested.summaryTable.SDRMSEOOF(i), ...
        nested.summaryTable.MeanRPDOOF(i));
end
fprintf(fid, '\nAggregate paired comparisons based on mean out-of-fold predictions:\n');
for i = 1:height(comparisonTable)
    fprintf(fid, '%s | %s minus %s | DeltaR2 %.4f [%.4f, %.4f] | DeltaRMSE %.4f [%.4f, %.4f]\n', ...
        comparisonTable{i, 1}{1}, comparisonTable{i, 2}{1}, ...
        comparisonTable{i, 3}{1}, comparisonTable{i, 4}, ...
        comparisonTable{i, 5}, comparisonTable{i, 6}, ...
        comparisonTable{i, 7}, comparisonTable{i, 8}, comparisonTable{i, 9});
end
fprintf(fid, '\nFull-band selected wavelength distribution:\n');
for i = 1:height(overlapTable)
    fprintf(fid, '%s | %s | selected %d/%d (%.1f%%); interval share %.1f%%\n', ...
        overlapTable{i, 1}{1}, overlapTable{i, 2}{1}, ...
        overlapTable{i, 5}, overlapTable{i, 6}, ...
        100 * overlapTable{i, 7}, 100 * overlapTable{i, 4});
end
end

function comparisonTable = bootstrapComparisons(y, predictions, names, pairs, types, nBoot)
rows = cell(size(pairs, 1), 12);
for c = 1:size(pairs, 1)
    candidateId = pairs(c, 1);
    referenceId = pairs(c, 2);
    candidateMet = metricsLocal(y, predictions(:, candidateId));
    referenceMet = metricsLocal(y, predictions(:, referenceId));
    deltaR2 = nan(nBoot, 1);
    deltaRMSE = nan(nBoot, 1);
    for b = 1:nBoot
        index = randi(numel(y), numel(y), 1);
        a = metricsLocal(y(index), predictions(index, candidateId));
        z = metricsLocal(y(index), predictions(index, referenceId));
        deltaR2(b) = a.R2 - z.R2;
        deltaRMSE(b) = a.RMSE - z.RMSE;
    end
    r2CI = quantile(deltaR2, [0.025 0.975]);
    rmseCI = quantile(deltaRMSE, [0.025 0.975]);
    rows(c, :) = {types{c}, names{candidateId}, names{referenceId}, ...
        candidateMet.R2 - referenceMet.R2, r2CI(1), r2CI(2), ...
        candidateMet.RMSE - referenceMet.RMSE, rmseCI(1), rmseCI(2), ...
        mean(deltaR2 > 0), mean(deltaRMSE < 0), nBoot};
end
comparisonTable = cell2table(rows, 'VariableNames', {'比较类型', '候选流程', ...
    '参照流程', 'DeltaR2', 'DeltaR2_CI95下限', 'DeltaR2_CI95上限', ...
    'DeltaRMSE', 'DeltaRMSE_CI95下限', 'DeltaRMSE_CI95上限', ...
    'R2改善频率', 'RMSE改善频率', 'Bootstrap次数'});
end

function [overlapTable, highFrequencyTable] = featureOverlap(bandTable)
q = bandTable(strcmp(bandTable.RangeName, 'Full_650_950') & ...
    ismember(bandTable.FeatureMethod, {'CARS', 'SPA', 'RandomFrog'}), :);
methods = {'CARS', 'SPA', 'RandomFrog', 'CARS与RandomFrog合计'};
segments = {'S1'; 'S2'; 'S3'; 'S4'; 'S3与S4合计'};
variableCounts = [20; 29; 40; 62; 102];
rows = cell(0, 8);
for m = 1:numel(methods)
    if m < 4
        methodRows = q(strcmp(q.FeatureMethod, methods{m}), :);
    else
        methodRows = q(ismember(q.FeatureMethod, {'CARS', 'RandomFrog'}), :);
    end
    total = height(methodRows);
    for s = 1:numel(segments)
        mask = segmentMask(methodRows.WavelengthNm, segments{s});
        rows(end + 1, :) = {methods{m}, segments{s}, variableCounts(s), ...
            variableCounts(s) / 151, sum(mask), total, sum(mask) / total, ...
            numel(unique(methodRows.WavelengthNm(mask)))}; %#ok<AGROW>
    end
end
overlapTable = cell2table(rows, 'VariableNames', {'变量策略', '区间', ...
    '区间变量数', '全波段变量占比', '选择次数', '该策略总选择次数', ...
    '选择占比', '所选不重复波长数'});

q = q(ismember(q.FeatureMethod, {'CARS', 'RandomFrog'}), :);
wavelengths = unique(q.WavelengthNm);
rows = cell(0, 6);
for i = 1:numel(wavelengths)
    w = wavelengths(i);
    carsCount = sum(q.WavelengthNm == w & strcmp(q.FeatureMethod, 'CARS'));
    rfCount = sum(q.WavelengthNm == w & strcmp(q.FeatureMethod, 'RandomFrog'));
    totalCount = carsCount + rfCount;
    if totalCount >= 2
        rows(end + 1, :) = {w, segmentName(w), carsCount, rfCount, totalCount, ...
            totalCount / 6}; %#ok<AGROW>
    end
end
highFrequencyTable = cell2table(rows, 'VariableNames', {'波长_nm', '所属区间', ...
    'CARS选择次数', 'RandomFrog选择次数', '合计选择次数', '六条选择链频率'});
highFrequencyTable = sortrows(highFrequencyTable, {'合计选择次数', '波长_nm'}, ...
    {'descend', 'ascend'});
end

function mask = segmentMask(wavelengths, segment)
switch segment
    case 'S1'
        mask = wavelengths >= 650 & wavelengths <= 688;
    case 'S2'
        mask = wavelengths >= 690 & wavelengths <= 746;
    case 'S3'
        mask = wavelengths >= 748 & wavelengths <= 826;
    case 'S4'
        mask = wavelengths >= 828 & wavelengths <= 950;
    otherwise
        mask = wavelengths >= 748 & wavelengths <= 950;
end
end

function name = segmentName(wavelength)
if wavelength <= 688
    name = 'S1';
elseif wavelength <= 746
    name = 'S2';
elseif wavelength <= 826
    name = 'S3';
else
    name = 'S4';
end
end

function result = metricsLocal(y, prediction)
residual = y - prediction;
result.R2 = 1 - sum(residual .^ 2) / sum((y - mean(y)) .^ 2);
result.RMSE = sqrt(mean(residual .^ 2));
result.RPD = std(y) / result.RMSE;
end
