function run_reproduction(stage)
% Reproduce the reported analysis. Start with run_reproduction('check').
if nargin < 1, stage = 'check'; end
stage = validatestring(stage, {'check','segmentation','fixed','refine','nested','boundary','summaries','all'});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'03_Code','Core_programs'));
addpath(fullfile(root,'03_Code','Dependencies'));
addpath(fullfile(root,'tests'));
if strcmp(stage,'check'), check_release(); return; end
if any(strcmp(stage,{'segmentation','all'})), stage2_l600_four_band_segmentation(); end
if any(strcmp(stage,{'fixed','all'})), stage2_l600_s3_models('all4test'); end
if any(strcmp(stage,{'refine','all'})), stage2_l600_s3_models('fourrefine'); end
if any(strcmp(stage,{'nested','all'})), stage2_l600_four_band_nested_validation('combined'); end
if any(strcmp(stage,{'boundary','all'})), stage2_l600_s3_boundary_nested_validation('combined'); end
if any(strcmp(stage,{'summaries','all'}))
    summarize_four_band_results();
    summarize_s3_boundary_results();
end
end
