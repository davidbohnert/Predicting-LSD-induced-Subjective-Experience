function recordPath = save_analysis_run_record(outputPath, recordPath, record)
%SAVE_ANALYSIS_RUN_RECORD Atomically save a compact analysis-run manifest.
%
%   recordPath = save_analysis_run_record(outputPath, '', record) creates a
%   timestamped run_record_*.txt file. Passing that path again updates the
%   same file, for example when a run changes from RUNNING to COMPLETE.
%
%   The canonical numerical results remain in their normal TXT/XLSX files;
%   this manifest records their paths and the parameters needed to identify
%   and reconstruct the run without duplicating the result tables.
%   Run records were added after submission and are operational metadata;
%   this utility is not part of CPM estimation or permutation inference.

    arguments
        outputPath (1, :) char
        recordPath (1, :) char
        record (1, 1) struct
    end

    if ~isfolder(outputPath)
        mkdir(outputPath);
    end

    if isempty(recordPath)
        runId = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
        recordPath = fullfile(outputPath, ['run_record_' runId '.txt']);
    end

    record.record_format = 'CPM_ANALYSIS_RUN_V1';
    record.record_updated_at = timestamp_now();
    if ~isfield(record, 'matlab_release')
        record.matlab_release = version('-release');
    end
    if ~isfield(record, 'matlab_version')
        record.matlab_version = version;
    end
    if isfield(record, 'entry_script') && ...
            ~isfield(record, 'entry_script_sha256') && ...
            isfile(record.entry_script)
        existingHash = read_existing_hash(recordPath);
        if isempty(existingHash)
            existingHash = sha256_file(record.entry_script);
        end
        record.entry_script_sha256 = existingHash;
    end
    record = attach_current_file_hash(record, 'results_text');
    record = attach_current_file_hash(record, 'results_workbook');

    temporaryPath = [tempname(outputPath) '.txt'];
    fileId = fopen(temporaryPath, 'wt');
    assert(fileId ~= -1, 'CPM:RunRecord:CannotOpen', ...
        'Could not create run record: %s', temporaryPath);

    try
        names = fieldnames(record);
        for fieldIdx = 1:numel(names)
            name = names{fieldIdx};
            fprintf(fileId, '%s=%s\n', name, encode_value(record.(name)));
        end
        fclose(fileId);
        fileId = -1;

        [moved, message] = movefile(temporaryPath, recordPath, 'f');
        assert(moved, 'CPM:RunRecord:CannotReplace', ...
            'Could not save run record %s: %s', recordPath, message);
    catch exception
        if fileId ~= -1
            fclose(fileId);
        end
        if isfile(temporaryPath)
            delete(temporaryPath);
        end
        rethrow(exception);
    end
end

function text = encode_value(value)
    if ischar(value)
        text = value;
    elseif isstring(value)
        text = strjoin(cellstr(value(:)), ',');
    elseif isnumeric(value) || islogical(value)
        text = mat2str(value);
    elseif iscell(value)
        parts = cellfun(@encode_value, value, 'UniformOutput', false);
        text = strjoin(parts, ',');
    elseif isdatetime(value)
        text = char(value);
    else
        error('CPM:RunRecord:UnsupportedValue', ...
            'Run-record values must be text, numeric, logical, cell, or datetime.');
    end

    text = strrep(text, char(13), ' ');
    text = strrep(text, newline, ' ');
end

function record = attach_current_file_hash(record, pathField)
    if ~isfield(record, pathField)
        return;
    end
    pathName = record.(pathField);
    if ~(ischar(pathName) && isfile(pathName))
        return;
    end
    hashField = [pathField '_sha256'];
    record.(hashField) = sha256_file(pathName);
end

function digestText = read_existing_hash(recordPath)
    digestText = '';
    if ~isfile(recordPath)
        return;
    end
    existingText = fileread(recordPath);
    match = regexp(existingText, ...
        '(?m)^entry_script_sha256=([0-9a-fA-F]{64})$', ...
        'tokens', 'once');
    if ~isempty(match)
        digestText = lower(match{1});
    end
end

function digestText = sha256_file(pathName)
    if ispc
        quoted = ['"' strrep(pathName, '"', '""') '"'];
        commands = {['certutil -hashfile ' quoted ' SHA256']};
    else
        singleQuote = char(39);
        escaped = [singleQuote '"' singleQuote '"' singleQuote];
        quoted = [singleQuote strrep(pathName, singleQuote, escaped) singleQuote];
        commands = {['sha256sum ' quoted], ['shasum -a 256 ' quoted], ...
            ['openssl dgst -sha256 ' quoted]};
    end
    digestText = '';
    for index = 1:numel(commands)
        [status, output] = system(commands{index});
        if status ~= 0, continue; end
        match = regexp(output, '(?i)([0-9a-f]{64})', 'tokens', 'once');
        if ~isempty(match), digestText = lower(match{1}); break; end
    end
    assert(~isempty(digestText), 'CPM:RunRecord:CannotHash', ...
        'No supported SHA-256 command was available for %s.', pathName);
end

function value = timestamp_now()
    value = char(datetime('now', ...
        'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
end
