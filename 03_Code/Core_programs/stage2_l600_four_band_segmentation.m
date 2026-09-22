function stage2_l600_four_band_segmentation()
% Compare globally optimized partitions directly from the thickness input workbook.

rng(20260817, 'twister');
scriptDir = fileparts(mfilename('fullpath'));
packageRoot = fileparts(fileparts(scriptDir));
inputFile = fullfile(packageRoot, '02_Data', 'L600_Data_Consolidated_EN_V02.xlsx');
outDir = fullfile(packageRoot, '04_Model_results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sheetNames = {'Original peel slice', '1 mm peel', '5 mm peel'};
excludedIds = string.empty(1, 0);
[wl, groupSpectra, groupIds] = load_l600_submission_thickness(inputFile, 'Thickness_spectra', sheetNames);
sampleCounts = cellfun(@(x) size(x, 1), groupSpectra);
if ~isequal(sampleCounts, [99 150 142])
    warning('Expected final counts [99 150 142], received [%s].', num2str(sampleCounts));
end
wl = wl(:);
candidateK = (2:8)';
minPoints = 15;
smoothWindow = 11;
bootstrapReplicates = 500;

[loadingScore, effectScore] = componentScores(groupSpectra);
responses = {
    '融合响应', smoothNormalize(0.5 .* loadingScore + 0.5 .* effectScore, smoothWindow); ...
    '厚度效应', smoothNormalize(effectScore, smoothWindow); ...
    'PCA载荷', smoothNormalize(loadingScore, smoothWindow)
    };

candidateRows = cell(0, 9);
partitions = struct();
for r = 1:size(responses, 1)
    responseName = responses{r, 1};
    y = responses{r, 2};
    key = ['Response' num2str(r)];
    previousSSE = sum((y - mean(y)).^2);
    for i = 1:numel(candidateK)
        k = candidateK(i);
        p = bestPartitionDP(y, k, minPoints);
        p.deltaSSE = (previousSSE - p.sse) / previousSSE;
        previousSSE = p.sse;
        partitions.(key).(['K' num2str(k)]) = p;
        candidateRows(end + 1, :) = {responseName, k, p.sse, p.bic, p.aicc, ...
            p.r2, p.deltaSSE, formatBoundaries(wl, p.breaks), formatSegments(wl, p.breaks)}; %#ok<AGROW>
    end
end
candidateTable = cell2table(candidateRows, 'VariableNames', ...
    {'响应类型', '区间数', '残差平方和', 'BIC', 'AICc', '分段解释率', ...
    '相对前一方案SSE降幅', '边界', '区间'});

% Sensitivity to the arbitrary fusion weight and smoothing width.
weights = [0.25 0.50 0.75];
windows = [7 11 15 21];
sensitivityRows = cell(0, 8);
for w = weights
    for window = windows
        y = smoothNormalize((1 - w) .* loadingScore + w .* effectScore, window);
        allParts = cell(numel(candidateK), 1);
        allBic = zeros(numel(candidateK), 1);
        for i = 1:numel(candidateK)
            allParts{i} = bestPartitionDP(y, candidateK(i), minPoints);
            allBic(i) = allParts{i}.bic;
        end
        [~, bestIndex] = min(allBic);
        best = allParts{bestIndex};
        sensitivityRows(end + 1, :) = {w, window, best.k, best.bic, best.r2, ...
            formatBoundaries(wl, best.breaks), formatSegments(wl, best.breaks), ...
            allBic(candidateK == 4) - min(allBic)}; %#ok<AGROW>
    end
end
sensitivityTable = cell2table(sensitivityRows, 'VariableNames', ...
    {'厚度效应权重', '平滑窗口点数', 'BIC最优区间数', '最优BIC', ...
    '分段解释率', '边界', '区间', '四段BIC减最优BIC'});

% Stratified record bootstrap: recompute PCA, thickness effect and the fusion response.
maxBreaks = max(candidateK) - 1;
bootstrapBreaks = nan(bootstrapReplicates, numel(candidateK), maxBreaks);
bootstrapBic = nan(bootstrapReplicates, numel(candidateK));
bootstrapSelectedK = nan(bootstrapReplicates, 1);
oobRMSE = nan(bootstrapReplicates, numel(candidateK));
oobSelectedK = nan(bootstrapReplicates, 1);
bootstrapBreaksThickness = nan(bootstrapReplicates, numel(candidateK), maxBreaks);
bootstrapBicThickness = nan(bootstrapReplicates, numel(candidateK));
bootstrapSelectedKThickness = nan(bootstrapReplicates, 1);
oobRMSEThickness = nan(bootstrapReplicates, numel(candidateK));
oobSelectedKThickness = nan(bootstrapReplicates, 1);
for b = 1:bootstrapReplicates
    bootGroups = cell(size(groupSpectra));
    oobGroups = cell(size(groupSpectra));
    for g = 1:numel(groupSpectra)
        n = size(groupSpectra{g}, 1);
        sampledIndex = randi(n, n, 1);
        oobIndex = setdiff((1:n)', unique(sampledIndex));
        bootGroups{g} = groupSpectra{g}(sampledIndex, :);
        oobGroups{g} = groupSpectra{g}(oobIndex, :);
    end
    [bootLoading, bootEffect] = componentScores(bootGroups);
    bootResponse = smoothNormalize(0.5 .* bootLoading + 0.5 .* bootEffect, smoothWindow);
    bootThickness = smoothNormalize(bootEffect, smoothWindow);
    [oobLoading, oobEffect] = componentScores(oobGroups);
    oobResponse = smoothNormalize(0.5 .* oobLoading + 0.5 .* oobEffect, smoothWindow);
    oobThickness = smoothNormalize(oobEffect, smoothWindow);
    for i = 1:numel(candidateK)
        p = bestPartitionDP(bootResponse, candidateK(i), minPoints);
        bootstrapBic(b, i) = p.bic;
        bootstrapBreaks(b, i, 1:numel(p.breaks)) = p.breaks;
        oobRMSE(b, i) = partitionPredictionRMSE(bootResponse, oobResponse, p.breaks);
        pThickness = bestPartitionDP(bootThickness, candidateK(i), minPoints);
        bootstrapBicThickness(b, i) = pThickness.bic;
        bootstrapBreaksThickness(b, i, 1:numel(pThickness.breaks)) = pThickness.breaks;
        oobRMSEThickness(b, i) = partitionPredictionRMSE(bootThickness, oobThickness, pThickness.breaks);
    end
    [~, selectedIndex] = min(bootstrapBic(b, :));
    bootstrapSelectedK(b) = candidateK(selectedIndex);
    [~, oobSelectedIndex] = min(oobRMSE(b, :));
    oobSelectedK(b) = candidateK(oobSelectedIndex);
    [~, selectedThicknessIndex] = min(bootstrapBicThickness(b, :));
    bootstrapSelectedKThickness(b) = candidateK(selectedThicknessIndex);
    [~, oobSelectedThicknessIndex] = min(oobRMSEThickness(b, :));
    oobSelectedKThickness(b) = candidateK(oobSelectedThicknessIndex);
end

selectionFrequency = arrayfun(@(k) mean(bootstrapSelectedK == k), candidateK);
selectionFrequencyThickness = arrayfun(@(k) mean(bootstrapSelectedKThickness == k), candidateK);
selectionTable = table(candidateK, selectionFrequency, selectionFrequencyThickness, ...
    'VariableNames', {'区间数', '融合响应Bootstrap选择频率', '厚度效应Bootstrap选择频率'});

oobMean = mean(oobRMSE, 1, 'omitnan')';
oobSd = std(oobRMSE, 0, 1, 'omitnan')';
oobSe = oobSd ./ sqrt(sum(isfinite(oobRMSE), 1))';
oobSelectionFrequency = arrayfun(@(k) mean(oobSelectedK == k), candidateK);
[minimumMean, minimumIndex] = min(oobMean);
oneSeThreshold = minimumMean + oobSe(minimumIndex);
withinOneSe = oobMean <= oneSeThreshold;
oneSeK = candidateK(find(withinOneSe, 1, 'first'));
oobTable = table(candidateK, oobMean, oobSd, oobSe, ...
    oobMean - minimumMean, withinOneSe, oobSelectionFrequency, ...
    'VariableNames', {'区间数', '袋外RMSE均值', '袋外RMSE标准差', ...
    '袋外RMSE标准误', '相对最小RMSE差值', '处于一标准误差范围', '袋外最优频率'});

oobMeanThickness = mean(oobRMSEThickness, 1, 'omitnan')';
oobSdThickness = std(oobRMSEThickness, 0, 1, 'omitnan')';
oobSeThickness = oobSdThickness ./ sqrt(sum(isfinite(oobRMSEThickness), 1))';
oobSelectionFrequencyThickness = arrayfun(@(k) mean(oobSelectedKThickness == k), candidateK);
[minimumMeanThickness, minimumIndexThickness] = min(oobMeanThickness);
oneSeThresholdThickness = minimumMeanThickness + oobSeThickness(minimumIndexThickness);
withinOneSeThickness = oobMeanThickness <= oneSeThresholdThickness;
oneSeKThickness = candidateK(find(withinOneSeThickness, 1, 'first'));
oobThicknessTable = table(candidateK, oobMeanThickness, oobSdThickness, oobSeThickness, ...
    oobMeanThickness - minimumMeanThickness, withinOneSeThickness, ...
    oobSelectionFrequencyThickness, 'VariableNames', {'区间数', '袋外RMSE均值', ...
    '袋外RMSE标准差', '袋外RMSE标准误', '相对最小RMSE差值', ...
    '处于一标准误差范围', '袋外最优频率'});

comparisonRows = cell(0, 8);
comparisonRows(end + 1, :) = pairedComparison('融合响应', 3, 4, oobRMSE, candidateK); %#ok<AGROW>
comparisonRows(end + 1, :) = pairedComparison('融合响应', 4, 5, oobRMSE, candidateK); %#ok<AGROW>
comparisonRows(end + 1, :) = pairedComparison('厚度效应', 3, 4, oobRMSEThickness, candidateK); %#ok<AGROW>
comparisonRows(end + 1, :) = pairedComparison('厚度效应', 4, 5, oobRMSEThickness, candidateK); %#ok<AGROW>
comparisonTable = cell2table(comparisonRows, 'VariableNames', ...
    {'响应类型', '较少区间数', '较多区间数', 'RMSE平均降低值', ...
    'RMSE降低值中位数', '百分之2_5分位', '百分之97_5分位', '较多区间改善频率'});

stabilityRows = cell(0, 11);
stabilityNames = {'融合响应', '厚度效应'};
stabilityKeys = {'Response1', 'Response2'};
stabilityArrays = {bootstrapBreaks, bootstrapBreaksThickness};
for responseIndex = 1:2
    for i = 1:numel(candidateK)
        k = candidateK(i);
        originalBreaks = partitions.(stabilityKeys{responseIndex}).(['K' num2str(k)]).breaks;
        for j = 1:(k - 1)
            values = squeeze(stabilityArrays{responseIndex}(:, i, j));
            midpointNm = (wl(values) + wl(values + 1)) ./ 2;
            originalMidpoint = (wl(originalBreaks(j)) + wl(originalBreaks(j) + 1)) ./ 2;
            ci = quantile(midpointNm, [0.025 0.25 0.5 0.75 0.975]);
            withinSixNm = mean(abs(midpointNm - originalMidpoint) <= 6);
            stabilityRows(end + 1, :) = {stabilityNames{responseIndex}, k, j, ...
                originalMidpoint, ci(3), ci(1), ci(5), ci(2), ci(4), ...
                withinSixNm, std(midpointNm)}; %#ok<AGROW>
        end
    end
end
stabilityTable = cell2table(stabilityRows, 'VariableNames', ...
    {'响应类型', '区间数', '边界序号', '原始边界中点_nm', 'Bootstrap中位数_nm', ...
    '百分之2_5分位_nm', '百分之97_5分位_nm', '百分之25分位_nm', ...
    '百分之75分位_nm', '落在原边界正负6nm频率', '标准差_nm'});

outXlsx = fullfile(outDir, '分段数量重新分析.xlsx');
if exist(outXlsx, 'file')
    delete(outXlsx);
end
writetable(candidateTable, outXlsx, 'Sheet', '候选段数');
writetable(sensitivityTable, outXlsx, 'Sheet', '参数敏感性');
writetable(selectionTable, outXlsx, 'Sheet', 'Bootstrap段数频率');
writetable(oobTable, outXlsx, 'Sheet', '融合响应袋外验证');
writetable(oobThicknessTable, outXlsx, 'Sheet', '厚度效应袋外验证');
writetable(comparisonTable, outXlsx, 'Sheet', '相邻方案袋外比较');
writetable(stabilityTable, outXlsx, 'Sheet', 'Bootstrap边界稳定性');

scoreTable = table(wl, loadingScore, effectScore, responses{1, 2}, responses{2, 2}, ...
    'VariableNames', {'波长_nm', 'PCA载荷标准化得分', '厚度效应标准化得分', ...
    '平滑融合响应', '平滑厚度效应'});
writetable(scoreTable, outXlsx, 'Sheet', '波长响应');

save(fullfile(outDir, '分段数量重新分析.mat'), 'wl', 'candidateK', ...
    'groupSpectra', 'groupIds', 'sampleCounts', 'excludedIds', ...
    'loadingScore', 'effectScore', 'responses', 'candidateTable', 'sensitivityTable', ...
    'selectionTable', 'stabilityTable', 'bootstrapBic', 'bootstrapBreaks', ...
    'bootstrapSelectedK', 'bootstrapBicThickness', 'bootstrapBreaksThickness', ...
    'bootstrapSelectedKThickness', 'oobTable', 'oobRMSE', 'oobSelectedK', ...
    'oneSeK', 'oneSeThreshold', 'oobThicknessTable', 'oobRMSEThickness', ...
    'oobSelectedKThickness', 'oneSeKThickness', 'oneSeThresholdThickness', ...
    'comparisonTable', 'partitions', 'minPoints', 'smoothWindow', 'bootstrapReplicates');

fid = fopen(fullfile(outDir, '分段数量重新分析摘要.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '候选区间数：2-8\n');
fprintf(fid, '每段最少波长点数：%d（%g nm）\n', minPoints, 2 * (minPoints - 1));
fprintf(fid, 'Bootstrap次数：%d\n\n', bootstrapReplicates);
for r = 1:size(responses, 1)
    name = responses{r, 1};
    subset = candidateTable(strcmp(candidateTable{:, 1}, name), :);
    [~, ix] = min(subset{:, 4});
    fprintf(fid, '%s的BIC最优方案：%d段，边界%s\n', ...
        name, subset{ix, 2}, subset{ix, 8}{1});
end
fprintf(fid, '\n融合响应Bootstrap段数选择频率：\n');
for i = 1:height(selectionTable)
    fprintf(fid, '%d段：%.3f\n', selectionTable{i, 1}, selectionTable{i, 2});
end
fprintf(fid, '\n厚度效应Bootstrap段数选择频率：\n');
for i = 1:height(selectionTable)
    fprintf(fid, '%d段：%.3f\n', selectionTable{i, 1}, selectionTable{i, 3});
end
fprintf(fid, '\n袋外验证最小平均RMSE对应：%d段\n', candidateK(minimumIndex));
fprintf(fid, '一标准误差规则选择：%d段\n', oneSeK);
fprintf(fid, '厚度效应袋外验证最小平均RMSE对应：%d段\n', candidateK(minimumIndexThickness));
fprintf(fid, '厚度效应一标准误差规则选择：%d段\n', oneSeKThickness);
thickness34 = comparisonTable(strcmp(comparisonTable{:, 1}, '厚度效应') & comparisonTable{:, 2} == 3, :);
thickness45 = comparisonTable(strcmp(comparisonTable{:, 1}, '厚度效应') & comparisonTable{:, 2} == 4, :);
fprintf(fid, '\n厚度效应3段到4段RMSE平均降低：%.6f，95%%区间[%.6f, %.6f]\n', ...
    thickness34{1, 4}, thickness34{1, 6}, thickness34{1, 7});
fprintf(fid, '厚度效应4段到5段RMSE平均降低：%.6f，95%%区间[%.6f, %.6f]\n', ...
    thickness45{1, 4}, thickness45{1, 6}, thickness45{1, 7});
fprintf(fid, '宽波段层面的推荐区间数：4段\n');
end

function [loadingScore, effectScore] = componentScores(groups)
X = vertcat(groups{:});
Z = zscore(X, 0, 1);
[coeff, ~, latent] = pca(Z);
explained = 100 .* latent ./ sum(latent);
cumExplained = cumsum(explained);
nPC = find(cumExplained >= 85, 1, 'first');
if isempty(nPC)
    nPC = numel(explained);
end
weights = explained(1:nPC) ./ sum(explained(1:nPC));
weightedLoading = abs(coeff(:, 1:nPC)) * weights;

nGroups = numel(groups);
nWavelengths = size(X, 2);
means = zeros(nGroups, nWavelengths);
sds = zeros(nGroups, nWavelengths);
counts = zeros(nGroups, 1);
for g = 1:nGroups
    counts(g) = size(groups{g}, 1);
    means(g, :) = mean(groups{g}, 1);
    sds(g, :) = std(groups{g}, 0, 1);
end
pooledVariance = zeros(1, nWavelengths);
for g = 1:nGroups
    pooledVariance = pooledVariance + (counts(g) - 1) .* (sds(g, :) .^ 2);
end
pooledVariance = pooledVariance ./ max(sum(counts) - nGroups, 1);
pooledSd = sqrt(max(pooledVariance, eps));
thicknessEffect = (max(means, [], 1) - min(means, [], 1))' ./ pooledSd';
loadingScore = normalize01(weightedLoading);
effectScore = normalize01(thicknessEffect);
end

function y = smoothNormalize(x, window)
y = sgolayfilt(x(:), 2, window);
y = normalize01(y);
end

function y = normalize01(x)
x = x(:);
lo = min(x);
hi = max(x);
if hi <= lo
    y = zeros(size(x));
else
    y = (x - lo) ./ (hi - lo);
end
end

function result = bestPartitionDP(score, k, minPoints)
n = numel(score);
prefix = [0; cumsum(score(:))];
prefixSq = [0; cumsum(score(:) .^ 2)];
cost = inf(k, n);
previous = zeros(k, n);

for last = minPoints:n
    cost(1, last) = intervalSSE(prefix, prefixSq, 1, last);
end
for segment = 2:k
    minLast = segment * minPoints;
    for last = minLast:n
        firstPreviousEnd = (segment - 1) * minPoints;
        lastPreviousEnd = last - minPoints;
        for previousEnd = firstPreviousEnd:lastPreviousEnd
            value = cost(segment - 1, previousEnd) + ...
                intervalSSE(prefix, prefixSq, previousEnd + 1, last);
            if value < cost(segment, last)
                cost(segment, last) = value;
                previous(segment, last) = previousEnd;
            end
        end
    end
end

breaks = zeros(1, k - 1);
last = n;
for segment = k:-1:2
    breaks(segment - 1) = previous(segment, last);
    last = breaks(segment - 1);
end
sse = cost(k, n);
totalSSE = sum((score - mean(score)).^2);
parameterCount = 2 * k - 1;
bic = n * log(max(sse / n, eps)) + parameterCount * log(n);
aic = n * log(max(sse / n, eps)) + 2 * parameterCount;
aicc = aic + 2 * parameterCount * (parameterCount + 1) / ...
    max(n - parameterCount - 1, 1);

result = struct('k', k, 'breaks', breaks, 'sse', sse, 'bic', bic, ...
    'aicc', aicc, 'r2', 1 - sse / totalSSE, 'deltaSSE', NaN);
end

function rmse = partitionPredictionRMSE(trainResponse, testResponse, breaks)
starts = [1, breaks + 1];
ends = [breaks, numel(trainResponse)];
predicted = zeros(size(testResponse));
for i = 1:numel(starts)
    predicted(starts(i):ends(i)) = mean(trainResponse(starts(i):ends(i)));
end
rmse = sqrt(mean((testResponse - predicted) .^ 2));
end

function row = pairedComparison(responseName, lowerK, higherK, rmseMatrix, candidateK)
lowerIndex = find(candidateK == lowerK, 1);
higherIndex = find(candidateK == higherK, 1);
difference = rmseMatrix(:, lowerIndex) - rmseMatrix(:, higherIndex);
ci = quantile(difference, [0.025 0.5 0.975]);
row = {responseName, lowerK, higherK, mean(difference), ci(2), ci(1), ci(3), ...
    mean(difference > 0)};
end

function value = intervalSSE(prefix, prefixSq, firstIndex, lastIndex)
n = lastIndex - firstIndex + 1;
sumValue = prefix(lastIndex + 1) - prefix(firstIndex);
sumSquares = prefixSq(lastIndex + 1) - prefixSq(firstIndex);
value = max(sumSquares - (sumValue ^ 2) / n, 0);
end

function textValue = formatBoundaries(wl, breaks)
parts = strings(1, numel(breaks));
for i = 1:numel(breaks)
    parts(i) = sprintf('%g/%g nm', wl(breaks(i)), wl(breaks(i) + 1));
end
textValue = char(strjoin(parts, '; '));
end

function textValue = formatSegments(wl, breaks)
starts = [1, breaks + 1];
ends = [breaks, numel(wl)];
parts = strings(1, numel(starts));
for i = 1:numel(starts)
    parts(i) = sprintf('%g-%g nm', wl(starts(i)), wl(ends(i)));
end
textValue = char(strjoin(parts, '; '));
end
