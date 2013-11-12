clear

datadir = getappdata(0, 'pigletdatadir');
studyDir = [datadir, filesep, 'post_conditioning'];
outputDir = [studyDir, filesep, 'output'];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

csch = uclColourScheme();
scrsz = get(0,'ScreenSize');

usedSystHeaders =  {'BP3 Mean';
                    'BP3 Rate';
                    'Temperature';
                    'BP3 Systolic';
                    'IPC Pressure';
                    'Laser Doppler 1'};
                
fileSuffix = 'eeg';

excelExpLogFile = [studyDir, filesep, 'eeg_log.xlsx'];
%% analyse
[pigData, pigDataNTB] = nirsAnalyseAndSync(studyDir,excelExpLogFile, usedSystHeaders, fileSuffix);

