function build_network_interaction_workbook(outputPath)
%BUILD_NETWORK_INTERACTION_WORKBOOK Create compact paper-facing XLSX tables.

if nargin < 1 || isempty(outputPath)
    scriptDir = fileparts([mfilename('fullpath') '.m']);
    repositoryRoot = fileparts(fileparts(scriptDir));
    outputPath = fullfile(repositoryRoot, 'outputs', 'network_interactions');
end
outputPath = char(outputPath);
pairPath = fullfile(outputPath, 'network_pair_results.csv');
summaryPath = fullfile(outputPath, 'outcome_summary.csv');
resultPath = fullfile(outputPath, 'network_interaction_results.mat');
workbookPath = fullfile(outputPath, 'network_interaction_results.xlsx');

assert(isfile(pairPath) && isfile(summaryPath) && isfile(resultPath), ...
    'CPM:NetworkInteraction:WorkbookInputsMissing', ...
    'CSV/MAT inputs must exist before workbook creation.');
pairResults = readtable(pairPath, 'TextType', 'string');
outcomeSummary = readtable(summaryPath, 'TextType', 'string');
loaded = load(resultPath, 'analysisSettings');
settings = loaded.analysisSettings;

pairColumns = {'Outcome', 'PairIndex', 'Network1', 'Network2', ...
    'NetworkPair', 'ObservedConsensusEdges', 'PossibleEdges', ...
    'NetworkPairDensity', 'ObservedExpectedRatio', 'CountPermutationP', ...
    'CountFDRQ', 'CountFDRSignificant', 'SparseObservedCount'};
assert(all(ismember(pairColumns, pairResults.Properties.VariableNames)), ...
    'CPM:NetworkInteraction:WorkbookPairColumns', ...
    'Network-pair CSV does not contain the required count-inference fields.');
pairResults = pairResults(:, pairColumns);

parameterRows = {
    'Analysis', settings.analysis;
    'Model', 'Covariate-adjusted LSD-PCB difference CPM, positive network';
    'Outcomes', strjoin(settings.outcomes, ', ');
    'Permutations per outcome', num2str(settings.no_iterations);
    'Workers requested', num2str(settings.parallel_workers_requested);
    'Workers active', num2str(settings.parallel_workers_active);
    'Cross-validation folds', num2str(settings.no_folds);
    'Feature-selection threshold', sprintf('p < %.3g', settings.edge_p_threshold);
    'Consensus criterion', sprintf('at least %d/%d folds', ...
        settings.consensus_folds_required, settings.no_folds);
    'Covariates', settings.covariates;
    'Subject subsets', 'GDE/VRS rows 1:67; MEQ30 rows 20:67';
    'Pair inferential statistic', 'Raw consensus-edge count (NCA style)';
    'Pair empirical p-value', settings.empirical_test;
    'Multiple comparisons', settings.multiple_comparison;
    'Displayed effect size', settings.display_effect_size;
    'Total-mask test', settings.total_mask_test;
    'Zero-edge null masks', settings.zero_edge_rule;
    'Fold RNG', settings.fold_rng;
    'Permutation RNG', settings.permutation_rng;
    'Permutation method', settings.permutation_method;
    };

if isfile(workbookPath), delete(workbookPath); end
writetable(pairResults, workbookPath, 'Sheet', 'Network Pairs');
writetable(outcomeSummary, workbookPath, 'Sheet', 'Outcome Summary', ...
    'WriteMode', 'overwritesheet');
writecell([{'Parameter', 'Value'}; parameterRows], workbookPath, ...
    'Sheet', 'Parameters', 'WriteMode', 'overwritesheet');

fprintf('Workbook data written: %s\n', workbookPath);
end
