function cpm_atomic_save(path, payload)
%CPM_ATOMIC_SAVE Replace a MAT checkpoint only after its complete write.
folder = fileparts(path);
if ~isfolder(folder), mkdir(folder); end
temporary = [tempname(folder) '.mat'];
cleanup = onCleanup(@() remove_temporary(temporary));
save(temporary, '-struct', 'payload', '-v7.3');
[ok,message] = movefile(temporary, path, 'f');
assert(ok, 'CPM:Checkpoint:Write', '%s', message);
clear cleanup
end

function remove_temporary(path)
if isfile(path), delete(path); end
end
