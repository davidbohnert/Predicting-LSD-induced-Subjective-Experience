function test_cv_wrappers
%TEST_CV_WRAPPERS Preserve output/seed conventions across the shared engine.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
rng(18,'twister');
mats = randn(6,6,24);
for i = 1:24
    a = mats(:,:,i); a = (a+a')/2; a(1:7:end) = 0; mats(:,:,i) = a;
end
behavior = randn(24,1);
settings = paper_analysis_config;
settings = settings.rng;
settings.cv_seed = 321;
[a,da] = permutation_test_cv_optimized('synthetic','LSD_VRS',mats, ...
    behavior,3,[.1 .2],[],'',{},2,'Pearson',0,settings);
[b,db] = permutation_test_cv_partition('synthetic','LSD_VRS',mats, ...
    behavior,3,[.1 .2],[],'',{},2,321,false,'Pearson',0,settings);
[c,dc] = permutation_test_cv_sample_size('synthetic','LSD_VRS',mats, ...
    behavior,3,[.1 .2],[],'',{},2,321,false,0,settings);
assert(isequaln(a,b) && isequaln(b,c));
assert(isequaln(da.perm_r,db.perm_r) && isequaln(db.perm_r,dc.perm_r));
assert(db.cv_seed==321 && dc.cv_seed==321);
fprintf('Shared CV wrapper outputs and seed propagation passed.\n');
end
