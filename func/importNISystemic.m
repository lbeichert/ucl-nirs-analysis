function S = importNISystemic(systFile)
%loads in systemic data files recorded using the NI system
%data is stored as a mat file

temp = load(systFile);
startTimeStr = datestr([temp.recInfo.expDate, temp.recInfo.startTime]);
tempdata = [temp.elapsed', temp.data(:,2:end)];
tempdata = interpolateToTimebase(tempdata, 1); %interpolate to 1 sample per second
S.date = startTimeStr(1:11);
S.start = startTimeStr(13:end);
S.elapsed = tempdata(:,1);
S.t = createDateNumFromElapsed(S.date, S.start, S.elapsed)';

calDate = '27/05/2013'; %date that calibration values were added
oldAnCal = [1,1,18,1,1,1,1,1];

%check code revision
if isfield(temp.recInfo, 'coderev')
    if temp.recInfo.coderev >= 0.1 %saved data is up-to-date and has conversion factors etc.
        
        
        S.headers = {temp.recInfo.anChanName{:}, temp.recInfo.digChanName{:}};
        S.chUnits = {temp.recInfo.anChanUnits{:}, repmat({'digital'}, 1, length(temp.recInfo.digChanName))};
        nAn = length(temp.recInfo.anChanName);
        
    end
else %old data, needs to be converted to physical units
    [~, anChanName, anChanCal, anChanOff, anChanUnits, ~, digChanName] = getNIChanSetup; %get acquisition parameters
    nAn = length(anChanName);
    S.headers = {anChanName{:}, digChanName{:}};
    S.chUnits = [anChanUnits{:}, repmat({'digital'},1,length(digChanName))];
    
    if datenum(S.date) >= datenum(calDate, 'dd/mm/yyyy')
        %data is in correct units
        
    else
        %data needs to be converted
        
        %convert everything back to voltages
        
        anData = tempdata(:,2:nAn+1);
        nt = size(anData,1); %number of timepoints
        anData = anData./repmat(oldAnCal,nt,1);
        
        %perform conversion
        anData = anData.*repmat(anChanCal,nt,1) + repmat(anChanOff,nt,1);
        
        tempdata(:,2:nAn+1) = anData; %put back into tempdata
    end
end

%round all digital channels following interpolation
tempdata(:,2+nAn:end) = round(tempdata(:,2+nAn:end));



S.data = tempdata(:,2:end);