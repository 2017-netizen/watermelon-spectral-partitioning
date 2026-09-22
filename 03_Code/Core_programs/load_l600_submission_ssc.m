function [wl, Y, Xraw, originalIndex] = load_l600_submission_ssc(inputFile, sheetName)
% Load the manuscript SSC matrix from the single consolidated workbook.

T = readtable(inputFile, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
originalIndex = T{:, 1};
Y = T{:, 2};
Xraw = T{:, 3:end};
labels = string(T.Properties.VariableNames(3:end));
wl = str2double(erase(labels, ' nm'))';

if any(isnan(wl)) || size(Xraw, 2) ~= numel(wl)
    error('The SSC wavelength columns in the consolidated workbook are invalid.');
end
end
