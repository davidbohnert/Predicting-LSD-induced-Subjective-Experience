function digest = file_sha256(pathName)
%FILE_SHA256 Portable, streaming SHA-256 digest for a local file.
assert(isfile(pathName), 'CPM:Hash:MissingFile', ...
    'Cannot hash missing file: %s', pathName);
quoted = quote_path(pathName);
if ispc
    commands = {['certutil -hashfile ' quoted ' SHA256']};
else
    commands = {['sha256sum ' quoted], ['shasum -a 256 ' quoted], ...
        ['openssl dgst -sha256 ' quoted]};
end
digest = '';
for index = 1:numel(commands)
    [status, output] = system(commands{index});
    if status ~= 0, continue; end
    match = regexp(output, '(?i)([0-9a-f]{64})', 'tokens', 'once');
    if ~isempty(match), digest = lower(match{1}); break; end
end
assert(~isempty(digest), 'CPM:Hash:Unavailable', ...
    'No supported SHA-256 command was available.');
end

function quoted = quote_path(pathName)
if ispc
    quoted = ['"' strrep(pathName, '"', '""') '"'];
else
    singleQuote = char(39);
    escaped = [singleQuote '"' singleQuote '"' singleQuote];
    quoted = [singleQuote strrep(pathName, singleQuote, escaped) singleQuote];
end
end
