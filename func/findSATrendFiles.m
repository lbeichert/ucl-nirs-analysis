function S = findSATrendFiles(pigletDir)

systDir = [pigletDir, filesep, 'systemic'];
systemicfolders = getFolders(systDir);
trendFiles = [];
%look for trend output files
for n = 1:length(systemicfolders)
    subFiles = getFiles(systemicfolders{n});
    [trendFiles{n}, ~, index{n}] = stringToken(subFiles, '(T|t)rend(.+).txt');
end

if isempty(trendFiles)
    S.present = false;
else
    S.present = true;
    S.files = trendFiles;
end
    