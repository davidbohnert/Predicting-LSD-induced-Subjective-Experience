function run_all_tests
%RUN_ALL_TESTS Run deterministic tests that do not launch paper analyses.
testDir = fileparts(mfilename('fullpath'));
addpath(testDir);
tests = {@test_publication_structure, @test_covariate_selection, ...
    @test_reporting, @test_table_numbering, @test_optimized_equivalence, @test_loso_equivalence, ...
    @test_cv_wrappers, @test_checkpoint_workflow};
for index = 1:numel(tests)
    fprintf('\n[%d/%d] %s\n', index, numel(tests), func2str(tests{index}));
    tests{index}();
end
fprintf('\nAll publication tests passed.\n');
end
