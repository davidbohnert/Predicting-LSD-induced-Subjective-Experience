function signature = cpm_checkpoint_guard(outputPath, settings, inputFiles, codeFiles)
%CPM_CHECKPOINT_GUARD Reject incompatible resumes before writing run outputs.
% Workers and presentation options are excluded by callers.
signature.version = 1;
signature.settings = settings;
signature.inputs = hashes(inputFiles, false);
signature.code = hashes(codeFiles, true);
path = fullfile(outputPath, 'checkpoint_signature.mat');
if isfile(path)
    saved = load(path, 'signature');
    assert(isfield(saved,'signature') && isequaln(saved.signature, signature), ...
        'CPM:Checkpoint:Incompatible', ...
        'Inputs, analysis settings, or code changed. Use a fresh output directory.');
else
    existing = dir(outputPath);
    names = string({existing.name});
    existing = existing(~ismember(names, [".","..",".DS_Store"]));
    assert(isempty(existing), 'CPM:Checkpoint:Unverified', ...
        'Existing outputs have no compatibility signature. Use a fresh output directory.');
    cpm_atomic_save(path, struct('signature',signature));
end
end

function result = hashes(files, includeName)
files = cellstr(string(files));
result = cell(numel(files),2);
for i = 1:numel(files)
    assert(isfile(files{i}), 'CPM:Checkpoint:MissingFile', ...
        'Missing analysis dependency: %s', files{i});
    if includeName
        [~,name,ext] = fileparts(files{i});
        result{i,1} = [name ext];
    else
        result{i,1} = i;
    end
    result{i,2} = file_sha256(files{i});
end
end
