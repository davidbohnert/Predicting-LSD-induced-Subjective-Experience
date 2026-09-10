function test_table_numbering
%TEST_TABLE_NUMBERING Verify renumbered parcellation exports and legacy alias.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'pipelines', 'parcellation_sensitivity'));
runRoot = tempname;
mkdir(runRoot);
cleanup = onCleanup(@() rmdir(runRoot, 's'));
outcomes = ["GDE";"MEQ30";"VRS";"OBN";"BDE";"AED"];
for atlas = ["original", "no_cerebellum", "with_cerebellum"]
    folder = fullfile(runRoot, atlas);
    mkdir(fullfile(folder, 'internal_validation'));
    mkdir(fullfile(folder, 'cross_drug_loso'));
    internal = table("LSD_" + outcomes, (1:6)'/10, (1:6)'/100, ...
        (11:16)', (21:26)', 'VariableNames', ...
        {'Analysis','r_positive','p_positive','mse_positive','positive_consensus_edges'});
    writetable(internal, fullfile(folder, 'internal_validation', ...
        'internal_validation_results.csv'));
    cross = table(repelem(outcomes(1:3),2), ...
        repmat(["other_psych";"all_amphs"],3,1), repmat(.01,6,1), ...
        (1:6)'/10, (1:6)'/100, repmat([50;46],3,1), ...
        'VariableNames', {'Behavior Scale','Testing Set','P Threshold', ...
        'r_pos','p_pos','N Test'});
    writetable(cross, fullfile(folder, 'cross_drug_loso', ...
        'cross_drug_loso_results.csv'));
    if atlas == "with_cerebellum"
        cerebellar = table(outcomes,(1:6)',(1:6)'/100,zeros(6,1), ...
            'VariableNames', {'Outcome','CerebellarEdges', ...
            'CerebellarEdgeFraction','WithinCerebellumEdges'});
        writetable(cerebellar, fullfile(folder, ...
            'cerebellar_consensus_edge_summary.csv'));
    end
end
legacyPath = fullfile(runRoot,'table_S12_internal.csv');
fid = fopen(legacyPath,'w'); fprintf(fid,'preserve existing result'); fclose(fid);
[a,b] = build_table_s13(runRoot);
assert(height(a)==18 && height(b)==18);
assert(isequal(a.R,repmat((1:6)'/10,3,1)));
assert(isequal(a.N,repmat([67;48;67;67;67;67],3,1)));
assert(isequal(b.NTestSessions,repmat([50;46],9,1)));
assert(isfile(fullfile(runRoot,'table_S13_internal.csv')));
assert(isfile(fullfile(runRoot,'table_S13_cross_drug.csv')));
[oldA,oldB] = build_table_s12(runRoot);
assert(isequaln(a,oldA) && isequaln(b,oldB));
assert(strcmp(fileread(legacyPath),'preserve existing result'));
fprintf('Table S13 exports and legacy builder compatibility passed.\n');
end
