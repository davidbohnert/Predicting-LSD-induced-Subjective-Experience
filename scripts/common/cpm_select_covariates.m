function [covars, columnNames, referenceStudy] = ...
        cpm_select_covariates(covarsIndicator, studyLabels)
%CPM_SELECT_COVARIATES Build a full-rank adjustment matrix for a sample.
%
% [COVARS, NAMES, REFERENCE] = CPM_SELECT_COVARIATES(INDICATORS, STUDIES)
% accepts rows from covars_indicator.mat and the corresponding study labels
% from column 1 of covars.mat.  The regression matrix always uses indicator
% columns, never the numeric study labels.  The lowest represented study is
% the local reference category.  Constant or linearly redundant nuisance
% columns are omitted deterministically so that [1 COVARS] is full rank.

validateattributes(covarsIndicator, {'numeric'}, ...
    {'2d','real','finite','ncols',5});
studyLabels = studyLabels(:);
assert(size(covarsIndicator, 1) == numel(studyLabels), ...
    'CPM:Covariates:Alignment', ...
    'Covariate rows and study labels must align.');
assert(all(ismember(studyLabels, [1 2 3])), ...
    'CPM:Covariates:StudyLabels', ...
    'Study labels must contain only 1, 2, or 3.');
assert(all(ismember(covarsIndicator(:, 1:3), [0 1]), 'all'), ...
    'CPM:Covariates:IndicatorCoding', ...
    'study_2, study_3, and sex must be binary.');
assert(isequal(covarsIndicator(:, 1), double(studyLabels == 2)) && ...
       isequal(covarsIndicator(:, 2), double(studyLabels == 3)), ...
    'CPM:Covariates:StudyMismatch', ...
    'Study indicators do not agree with the grouping labels.');

representedStudies = unique(studyLabels, 'sorted');
referenceStudy = representedStudies(1);
candidateColumns = [];
candidateNames = strings(0, 1);
for study = representedStudies(:)'
    if study == referenceStudy
        continue
    end
    assert(study > 1, 'CPM:Covariates:ReferenceCoding', ...
        'A non-reference study lacks an indicator column.');
    candidateColumns(end + 1) = study - 1; %#ok<AGROW>
    candidateNames(end + 1, 1) = "study_" + study; %#ok<AGROW>
end
candidateColumns = [candidateColumns 3 4 5];
candidateNames = [candidateNames; "sex"; "age"; "mean_FD"];

design = ones(numel(studyLabels), 1);
kept = false(1, numel(candidateColumns));
for index = 1:numel(candidateColumns)
    proposed = [design covarsIndicator(:, candidateColumns(index))]; %#ok<AGROW>
    if rank(proposed) > rank(design)
        design = proposed;
        kept(index) = true;
    end
end
covars = covarsIndicator(:, candidateColumns(kept));
columnNames = cellstr(candidateNames(kept));
assert(rank([ones(size(covars, 1), 1) covars]) == size(covars, 2) + 1, ...
    'CPM:Covariates:RankDeficient', ...
    'The selected nuisance-variable design is not full rank.');
end
