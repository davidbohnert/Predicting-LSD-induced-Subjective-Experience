function test_reporting
%TEST_REPORTING Verify BH correction and fixed-prediction BCa output.
testDir = fileparts(mfilename('fullpath'));
root = fileparts(testDir);
addpath(fullfile(root,'scripts','utilities'));
p = [.01 .04 .03 .20];
q = benjamini_hochberg(p);
assert(max(abs(q-[.04 .053333333333333 .053333333333333 .20])) < 1e-12);
rng(41,'twister');
observed = randn(40,1);
predicted = observed + randn(40,1);
interval = participant_bootstrap_ci(observed,predicted,200,1300040);
assert(all(isfinite(interval.bca)) && interval.bca(1)<interval.bca(2));
rng(1300040,'twister');
expectedDraws = randi(numel(observed), 200, numel(observed));
assert(isequal(interval.draws, expectedDraws), ...
    'Bootstrap draws do not reproduce the authoritative matrix ordering.');
fprintf('Publication reporting utilities passed.\n');
end
