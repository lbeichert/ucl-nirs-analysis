clear
addpath([filesep,'func']);
% 
datadir = getappdata(0, 'pigletdatadir');
studyDir = [datadir, filesep, 'argon-dex'];
outputDir = [studyDir, filesep, 'output',filesep,'31p'];

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%csch = uclColourScheme();
%scrsz = get(0,'ScreenSize');

usedSystHeaders =  {'BP3 Mean';
                    'BP3 Rate';
                    'Temperature';
                    'BP3 Systolic';
                    'IPC Pressure';
                    'Laser Doppler 1'};
                
fileSuffix = '31p';

excelExpLogFile = [studyDir, filesep, '31p_log.xlsx'];
%% analyse
[pigData, pigDataNTB] = nirsAnalyseAndSync(studyDir,excelExpLogFile, usedSystHeaders, fileSuffix);
[phosData, pigDataPTB] = phosAnalyseAndSync(studyDir, excelExpLogFile, pigDataNTB);

save([outputDir,filesep,'pigData'], 'pigData', 'pigDataNTB');
save([outputDir,filesep,'pigDataPTB'], 'pigDataNTB');