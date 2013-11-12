function [systFiles] = findNISystemicFiles(pigletDir)
systDir = [pigletDir, filesep, 'systemic'];
systemicFolderFiles = getFiles(systDir);


systFiles.mat = stringToken(systemicFolderFiles, 'systemic_(.+).mat');


if isempty(systFiles.mat)
    systFiles.present=false;
else
    systFiles.present=true;
end