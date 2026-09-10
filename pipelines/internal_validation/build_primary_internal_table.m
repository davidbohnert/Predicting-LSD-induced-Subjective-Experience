function paperTable = build_primary_internal_table( ...
        rawResults, details, noBootstraps)
%BUILD_PRIMARY_INTERNAL_TABLE Table 1 statistics from fixed OOF predictions.
if nargin < 3 || isempty(noBootstraps), noBootstraps = 10000; end
assert(height(rawResults) == 6 && numel(details) == 6, ...
    'CPM:Reporting:PrimaryRowCount', ...
    'The primary table requires exactly one result for each of six outcomes.');

preferredOrder = {'GDE','MEQ30','VRS','OBN','BDE','AED'};
outcome = strings(6,1);
n = zeros(6,1); r = zeros(6,1); p = zeros(6,1);
ciLower = zeros(6,1); ciUpper = zeros(6,1);
mse = zeros(6,1); mseP = zeros(6,1); edges = zeros(6,1);
covariateColumns = strings(6,1); studyReference = zeros(6,1);
for row = 1:6
    index = find(strcmp({details.outcome}, preferredOrder{row}));
    assert(isscalar(index), 'CPM:Reporting:OutcomeMatch', ...
        'Could not uniquely match %s.', preferredOrder{row});
    engine = details(index).engine_details;
    observed = engine.observed_behavior;
    predicted = engine.observed_predictions.positive;
    assert(isvector(predicted), 'CPM:Reporting:PredictionShape', ...
        'Primary predictions must contain one threshold.');
    interval = participant_bootstrap_ci(observed, predicted, ...
        noBootstraps, 1300000 + numel(observed));
    outcome(row) = preferredOrder{row};
    n(row) = numel(observed);
    r(row) = rawResults.r_positive(index);
    p(row) = rawResults.p_positive(index);
    ciLower(row) = interval.bca(1);
    ciUpper(row) = interval.bca(2);
    mse(row) = rawResults.mse_positive(index);
    mseP(row) = rawResults.mse_p_positive(index);
    edges(row) = rawResults.positive_consensus_edges(index);
    covariateColumns(row) = strjoin(details(index).covariate_columns, ', ');
    studyReference(row) = details(index).study_reference;
end
q = benjamini_hochberg(p);
paperTable = table(outcome,n,r,p,q,ciLower,ciUpper,mse,mseP,edges, ...
    covariateColumns,studyReference, 'VariableNames', ...
    {'Outcome','N','R','P','Q','BCa95Lower','BCa95Upper','MSE','MSEP', ...
     'PositiveConsensusEdges','CovariateColumns','StudyReference'});
end
