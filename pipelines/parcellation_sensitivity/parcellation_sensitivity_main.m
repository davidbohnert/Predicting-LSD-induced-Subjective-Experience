function outputs = parcellation_sensitivity_main(cfg, atlasName, options)
%PARCELLATION_SENSITIVITY_MAIN One atlas condition for Table S13.
%
% Run this entry point separately for "original", "no_cerebellum", and
% "with_cerebellum".  Each invocation performs internal validation and
% subject-wise cross-drug LOSO sequentially and writes to its own directory.
arguments
    cfg (1,1) struct
    atlasName (1,1) string
    options.no_iterations (1,1) double = NaN
    options.workers (1,1) double = NaN
    options.output_path (1,1) string = ""
    options.connectome_directory (1,1) string = ""
end
atlasName = validatestring(atlasName, ...
    {'original','no_cerebellum','with_cerebellum'});
if isnan(options.no_iterations), iterations = cfg.runtime.main_iterations;
else, iterations = options.no_iterations; end
if isnan(options.workers), workers = cfg.runtime.workers;
else, workers = options.workers; end
if options.output_path == ""
    outputPath = fullfile(cfg.output_root, 'parcellation_sensitivity', atlasName);
else
    outputPath = char(options.output_path);
end

switch atlasName
    case 'original'
        preset = "primary";
        suffix = "";
        defaultDirectory = fullfile(cfg.input_root, 'connectomes');
    case 'no_cerebellum'
        preset = "altparcel_no_cerebellum";
        suffix = "_altparcel_no_cerebellum";
        defaultDirectory = cfg.altparcel.no_cerebellum_directory;
    case 'with_cerebellum'
        preset = "altparcel_with_cerebellum";
        suffix = "_altparcel_with_cerebellum";
        defaultDirectory = cfg.altparcel.with_cerebellum_directory;
end
if options.connectome_directory == ""
    connectomeDirectory = defaultDirectory;
else
    connectomeDirectory = char(options.connectome_directory);
end
assert(isfolder(connectomeDirectory), 'CPM:Parcellation:MissingDirectory', ...
    'Connectome directory does not exist: %s', connectomeDirectory);

internalPath = fullfile(outputPath, 'internal_validation');
crossDrugPath = fullfile(outputPath, 'cross_drug_loso');
internal = internal_validation_main(cfg, preset, ...
    no_iterations=iterations, workers=workers, output_path=internalPath, ...
    connectome_directory=connectomeDirectory);
crossDrug = cross_drug_loso_main(cfg, no_iterations=iterations, ...
    workers=workers, output_path=crossDrugPath, ...
    scales={'GDE','MEQ30','VRS'}, thresholds=0.01, ...
    connectome_directory=connectomeDirectory, connectome_suffix=suffix);

if strcmp(atlasName, 'with_cerebellum')
    edgeSummary = summarize_cerebellar_edges( ...
        internalPath, cfg.altparcel.with_cerebellum_metadata);
    writetable(edgeSummary, fullfile(outputPath, ...
        'cerebellar_consensus_edge_summary.csv'));
else
    edgeSummary = table;
end
outputs = struct('internal',internal, 'cross_drug',{crossDrug}, ...
    'cerebellar_edges',edgeSummary, 'atlas',atlasName, ...
    'output_path',outputPath);
end

function summary = summarize_cerebellar_edges(internalPath, metadataPath)
assert(isfile(metadataPath), 'CPM:Parcellation:MissingAtlasMetadata', ...
    'The 439-node structure metadata file is missing: %s', metadataPath);
parcelMetadata = readtable(metadataPath, 'TextType','string');
assert(isequal(parcelMetadata.Properties.VariableNames, ...
    {'ParcelIndex','Structure'}) && height(parcelMetadata) == 439 && ...
    isequal(sort(parcelMetadata.ParcelIndex), (1:439)'), ...
    'CPM:Parcellation:AtlasMetadataSchema', ...
    'Atlas metadata must describe each parcel index 1:439 exactly once.');
cerebellarNodes = parcelMetadata.ParcelIndex( ...
    strcmpi(parcelMetadata.Structure, 'Cerebellum'));
assert(numel(cerebellarNodes) == 7, ...
    'CPM:Parcellation:CerebellarParcelCount', ...
    'Expected exactly seven cerebellar parcels in the 439-node atlas.');
loaded = load(fullfile(internalPath, 'internal_validation_results.mat'), ...
    'details');
details = loaded.details;
rows = cell(numel(details), 5);
for index = 1:numel(details)
    consensus = details(index).consensus{1}.observed_pos_mask;
    assert(size(consensus,1) == height(parcelMetadata), ...
        'CPM:Parcellation:NodeCount', ...
        'The consensus mask and atlas metadata have different node counts.');
    [nodeA,nodeB] = find(triu(consensus,1));
    nodeAIsCerebellar = ismember(nodeA, cerebellarNodes);
    nodeBIsCerebellar = ismember(nodeB, cerebellarNodes);
    involvesCerebellum = nodeAIsCerebellar | nodeBIsCerebellar;
    total = numel(nodeA);
    cerebellar = nnz(involvesCerebellum);
    rows(index,:) = {details(index).outcome, total, cerebellar, ...
        cerebellar / max(total,1), ...
        nnz(nodeAIsCerebellar & nodeBIsCerebellar)};
end
summary = cell2table(rows, 'VariableNames', ...
    {'Outcome','PositiveConsensusEdges','CerebellarEdges', ...
     'CerebellarEdgeFraction','WithinCerebellumEdges'});
end
