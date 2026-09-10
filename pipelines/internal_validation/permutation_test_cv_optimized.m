function [results, details] = permutation_test_cv_optimized(varargin)
%PERMUTATION_TEST_CV_OPTIMIZED Internal-validation interface to shared CPM.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'scripts', 'common'));
[results, details] = cpm_permutation_cv(varargin{:});
end
