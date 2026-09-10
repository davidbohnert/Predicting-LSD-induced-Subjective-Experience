function [inputs, code] = cpm_analysis_dependencies(cfg, pipelineDirs, extraInputs)
%CPM_ANALYSIS_DEPENDENCIES Content dependencies for reproducible checkpoints.
% Include behavior fallbacks so adding a preferred file invalidates a run.
if nargin < 3, extraInputs = {}; end
inputs = {};
folders = {fullfile(cfg.input_root,'behav')};
for i = 1:numel(folders)
    entries = dir(fullfile(folders{i}, '*.mat'));
    for j = 1:numel(entries)
        inputs{end+1} = fullfile(entries(j).folder,entries(j).name); %#ok<AGROW>
    end
end
inputs = [inputs, {cfg.model.covariate_file, cfg.model.grouping_file}, ...
    reshape(cellstr(string(extraInputs)),1,[])];
inputs = unique(inputs, 'stable');
root = cfg.repository_root;
folders = [{fullfile(root,'scripts','common')}, ...
    reshape(cellstr(string(pipelineDirs)),1,[])];
code = {};
for i = 1:numel(folders)
    entries = dir(fullfile(folders{i}, '*.m'));
    for j = 1:numel(entries)
        if startsWith(entries(j).name, {'render_','verify_'})
            continue
        end
        code{end+1} = fullfile(entries(j).folder,entries(j).name); %#ok<AGROW>
    end
end
code = [code, {fullfile(root,'scripts','utilities','participant_bootstrap_ci.m'), ...
    fullfile(root,'scripts','utilities','benjamini_hochberg.m')}];
end
