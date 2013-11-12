function [poxFiles] = findPulseOxFiles(pigletDir)
systDir = [pigletDir, filesep, 'systemic'];
systemicFolderFiles = getFiles(systDir);

poxFiles.dat = stringToken(systemicFolderFiles, 'pulseox(.+).dat');
poxFiles.mat = stringToken(systemicFolderFiles, 'pulseox(.+).mat');
if all([isempty(poxFiles.dat), isempty(poxFiles.mat)])
    poxFiles.mat = stringToken(systemicFolderFiles, 'po_(.+).mat');
end

if all([isempty(poxFiles.dat), isempty(poxFiles.mat)])
    poxFiles.present=false;
else
    poxFiles.present=true;
end