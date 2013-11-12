function dataStructDLMWrite(S, fname)
%writes a nirs/physiology/31P standard data structure to disk

fid = fopen(fname, 'wt');


headers = {'time', 'elapsed (s)'};
headers = [headers, S.headers];
nHeaders = length(headers);
for n = 1:nHeaders
    fprintf(fid, '%s\t', headers{n});
end
fprintf(fid, '\n');
fclose(fid);

out = S.t;
out = [out, S.elapsed];
out = [out, S.data];

dlmwrite(fname, out, '-append', 'delimiter', '\t', 'precision', '%.6f');
