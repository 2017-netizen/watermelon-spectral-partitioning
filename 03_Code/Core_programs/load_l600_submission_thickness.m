function [wl, groupSpectra, groupIds] = load_l600_submission_thickness(inputFile, sheetName, treatmentNames)
% Load the three peel-thickness groups from the single consolidated workbook.

T = readtable(inputFile, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
labels = string(T.Properties.VariableNames(3:end));
wl = str2double(erase(labels, ' nm'))';
treatment = string(T{:, 1});
recordId = string(T{:, 2});
groupSpectra = cell(1, numel(treatmentNames));
groupIds = cell(1, numel(treatmentNames));

for g = 1:numel(treatmentNames)
    keep = treatment == string(treatmentNames{g});
    groupSpectra{g} = T{keep, 3:end};
    groupIds{g} = recordId(keep);
end

if any(isnan(wl)) || any(cellfun(@isempty, groupSpectra))
    error('The peel-thickness data in the consolidated workbook are invalid.');
end
end
