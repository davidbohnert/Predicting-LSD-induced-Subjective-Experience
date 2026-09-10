function test_covariate_selection
%TEST_COVARIATE_SELECTION Verify study-reference coding in common samples.
testDir = fileparts(mfilename('fullpath'));
root = fileparts(testDir);
addpath(fullfile(root, 'scripts', 'common'));

studies = [ones(3,1); 2*ones(4,1); 3*ones(5,1)];
rng(77,'twister');
indicator = [studies==2, studies==3, mod((1:12)',2), ...
    randn(12,1), randn(12,1)];
[full, names, reference] = cpm_select_covariates(indicator, studies);
assert(reference==1 && isequal(names(:), ...
    {'study_2';'study_3';'sex';'age';'mean_FD'}));
assert(rank([ones(12,1),full])==6);

keep = studies~=1;
[twoStudy, names, reference] = ...
    cpm_select_covariates(indicator(keep,:), studies(keep));
assert(reference==2 && isequal(names(:), ...
    {'study_3';'sex';'age';'mean_FD'}));
assert(rank([ones(sum(keep),1),twoStudy])==5);

keep = studies==3;
[oneStudy, names, reference] = ...
    cpm_select_covariates(indicator(keep,:), studies(keep));
assert(reference==3 && ~any(startsWith(names,'study_')));
assert(rank([ones(sum(keep),1),oneStudy])==size(oneStudy,2)+1);
fprintf('Full-rank covariate selection passed.\n');
end
