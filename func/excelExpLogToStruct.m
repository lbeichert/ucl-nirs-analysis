function expLog = excelExpLogToStruct(excelExpLogFile)
%converts an excel study log file into a structure with the following
%fields:
%pigNum - LWP number
%expDate - experiment start date (date of day 1 of the experiment)
%exposure - NIRS exposure duration in seconds
%nirsStart - nirs start time of the 
%refFiles - nirs reference file name
%nirsFiles - nirs data file name
%systFiles - systemic data file name
%poFiles - pulse oximetry data file name
%group
%excludeAfter - seconds to exclude after

studyDir = fileparts(excelExpLogFile);
[status, sheets, ~] = xlsfinfo(excelExpLogFile);

if isempty(status)
    error('not a valid excel spreadsheet');
    return
end



NSheets = length(sheets);


for cs = 1:NSheets
    
    [~,txtStudyLog,rawStudyLog] = xlsread(excelExpLogFile,sheets{cs});
    
    
    expLog(cs).pigNum = [rawStudyLog{3:end,1}];
    expLog(cs).expDate = {rawStudyLog{3:end,2}};
    expLog(cs).exposure = [rawStudyLog{3:end,3}]; %exposure time
    expLog(cs).nirsStart = [rawStudyLog{3:end,4}];
    expLog(cs).refFiles = {rawStudyLog{3:end,5}};
    expLog(cs).nirsFiles = {rawStudyLog{3:end,6}};
    expLog(cs).systFiles = {rawStudyLog{3:end,7}};
    expLog(cs).poFiles = {rawStudyLog{3:end,8}};
    expLog(cs).group = {rawStudyLog{3:end,9}};
    expLog(cs).excludeAfter = [rawStudyLog{3:end,10}];
    expLog(cs).notes = {rawStudyLog{3:end,11}};
    
    if size(txtStudyLog,2) > 11
        for n = 1:size(txtStudyLog,2) - 11
            expLog(cs).(txtStudyLog{2,11+n}) = [rawStudyLog{3:end,11+n}];
        end
    end
    
    NPigs = length(expLog(cs).pigNum);
    pigs = 1:length(expLog(cs).pigNum);
    pigNum = expLog(cs).pigNum;
    
    for cp = pigs
        expLog(cs).subj{cp} = ['LWP', num2str(pigNum(cp))];
        expLog(cs).pigletDir{cp} = [studyDir, filesep, expLog(cs).subj{cp}];
        
        nirsDir = [expLog(cs).pigletDir{cp}, filesep, 'nirs'];
        systDir = [expLog(cs).pigletDir{cp}, filesep, 'systemic'];
        
        if ~isnan(expLog(cs).refFiles{cp})
            expLog(cs).refFiles{cp} = [nirsDir, filesep, expLog(cs).refFiles{cp}];
        end
        if ~isnan(expLog(cs).nirsFiles{cp})
            expLog(cs).nirsFiles{cp} = [nirsDir, filesep, expLog(cs).nirsFiles{cp}];
        end
        if ~isnan(expLog(cs).systFiles{cp})
            expLog(cs).systFiles{cp} = [systDir, filesep, expLog(cs).systFiles{cp}];
            
        end
        if ~isnan(expLog(cs).poFiles{cp})
            expLog(cs).poFiles{cp} = [systDir, filesep, expLog(cs).poFiles{cp}];
        end
    end
end