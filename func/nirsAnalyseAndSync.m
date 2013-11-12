function [pigData, pigDataNTB] = nirsAnalyseAndSync(studyDir, excelExpLogFile, usedSystHeaders, fileSuffix)

% processes NIRS and physiological data from the insult, synchronises the
% data and produces graphs for initial viewing/analysis
%%

[status, sheets, format] = xlsfinfo(excelExpLogFile);

if isempty(status)
    error('not a valid excel spreadsheet');
    return
end



NSheets = length(sheets);
cs=1;
cp=1;

%%
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
    
    NPigs = length(pigNum);
    pigs = 1:length(pigNum);
    
    for cp = pigs
        subj{cp} = ['LWP', num2str(pigNum(cp))];
        pigletDir{cp} = [studyDir, filesep, subj{cp}];
        
        refFiles{cp} = [pigletDir{cp}, filesep, 'nirs', filesep, refFiles{cp}];
        if ~isnan(nirsFiles{cp})
            nirsFiles{cp} = [pigletDir{cp}, filesep, 'nirs',filesep, nirsFiles{cp}];
        end
        if ~isnan(systFiles{cp})
            systFiles{cp} = [pigletDir{cp}, filesep, 'systemic',filesep, systFiles{cp}];
            
        end
        if ~isnan(poFiles{cp})
            poFiles{cp} = [pigletDir{cp}, filesep, 'systemic',filesep, poFiles{cp}];
        end
    end
    
    for cp = pigs
        disp([sheets{cs}, subj{cp}]);
        
        
        pigData(cp) = loadPigData(refFiles{cp}, nirsFiles{cp}, systFiles{cp}, poFiles{cp}, nirsStart(cp), exposure(cp), expDate{cp}, usedSystHeaders); %load in all corresponding data
        close all;
        if ~isempty( pigData(cp).pulseox)
            pigData(cp) = despikePulseOx(pigData(cp));
        end
        
        if ~isempty(pigData(cp).nirs)
            pigDataNTB(cp) = resamplePigData(pigData(cp), pigletDir{cp}, subj{cp}, fileSuffix, true); %resample pig data to the piglet
        end
    end
    
    
    
end