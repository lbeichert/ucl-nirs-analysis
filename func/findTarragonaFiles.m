function [tdf, fileinfo] = findTarragonaFiles(pigletDir)
	%find tarragona data files inside the parent piglet directory

	nirsDir = [pigletDir, filesep, 'nirs'];
	nirsFolderFiles = getFiles(nirsDir);
	nirsFilesInfo = dir(nirsDir);

	%find tarragona data files in the nirs directory
	[tdf, ~, index] = stringToken(nirsFolderFiles, '(.+).tdf');

	fileinfo = nirsFilesInfo(index);