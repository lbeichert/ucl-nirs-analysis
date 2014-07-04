function [phosData, pigDataPTB] = phosAnalyseAndSync(studyDir, excelExpLogFile, pigDataNTB)
% PHOSANALYSEANDSYNC loads 31P-data from the insult and synchronises it
% with NIRS and systemic data

[status, sheets, format] = xlsfinfo(excelExpLogFile);

if isempty(status)
    error('not a valid excel spreadsheet');
    return
end


NSheets = length(sheets);

for cs = 1:NSheets
    
    [~,txtStudyLog,rawStudyLog] = xlsread(excelExpLogFile,sheets{cs});
    pigNum = [rawStudyLog{3:end,1}];
    expDate = {rawStudyLog{3:end,2}};
    exposure = [rawStudyLog{3:end,3}]; %exposure time
    nirsStart = [rawStudyLog{3:end,4}];
    refFiles = {rawStudyLog{3:end,5}};
    nirsFiles = {rawStudyLog{3:end,6}};
    systFiles = {rawStudyLog{3:end,7}};
    poFiles = {rawStudyLog{3:end,8}};
    group = {rawStudyLog{3:end,9}};
    excludeAfter = [rawStudyLog{3:end,10}];
    notes = {rawStudyLog{3:end,11}};
    phosFiles = {rawStudyLog{3:end,12}};
    
    NPigs = length(pigNum);
    pigs = 1:length(pigNum);
    
    for cp = pigs
        subj{cp} = ['LWP', num2str(pigNum(cp))];
        pigletDir{cp} = [studyDir, filesep, subj{cp}];
        phosFiles{cp} = [pigletDir{cp}, filesep, 'phosphorus', filesep, phosFiles{cp}];
    end
    
    % load phosphorus data
    for cp = pigs
        phosData(cp) = load31PData(phosFiles{cp});
    end
    
    % resample with pigDataNTB
    for cp = pigs
        pigDataPTB(cp) = resamplePigDataPTB(pigDataNTB(cp), phosData(cp));
    end
    
end
end