function [out, S] = importPulseOx(poxFiles)
% loads in pulse oximetry data files
% pulse oximetry files are in one of three formats:
% 1.     .dat, tabular text data, where waveform data is also sampled (only
%        a few LWP's).
% 2.     .mat old style, from the standalone pulse ox logging script.
% 3.     .mat new style, from the combined physiological monitoring logging
%        script
% This script is hopefully smart enought to determine what format the data
% is in...
% Usage
% poxData = importPulseOx(poxFiles)
% [poxData S] = importPulse(poxFiles)
% arguments:
%   poxFiles:   cell array of pulse oximetry files
% returns
%   poxData:    An m-by-3 array of all the pulse oximetry files
%               synchronised to the same start time.  Column 1 is the
%               serial data number of the sample, column 2 is the heartrate
%               and column 3 is the oxygen saturation.
%   S:          a structure with the raw data for each file, containing the
%               following fields:
%       data:       heartrate and saturation array of samples for that pulse
%                   oximetry file
%       date:       the date the pulse oximetry data begins on
%       start:      the start time of the pulse oximetry data
%       elapsed:    vector of elapsed time since the start
%       t:          vector of serial date numbers (absolute times) for each sample

%load in pulse oximetry data

if ~isstruct(poxFiles)
    fname = poxFiles;
    poxFiles = [];
    [~, ~, ext] = fileparts(fname);
    if strcmp(ext, '.mat')
        poxFiles.mat{1} = fname;
    end
    if strcmp(ext, '.dat')
        poxFiles.dat{1} = fname;
    end
end
    
if isempty(poxFiles.mat) && length(poxFiles.mat) ~= length(poxFiles.dat) %load in the .dat text files
    DELIMITER = '\t';
    HEADERLINES = 2;
    NFiles = length(poxFiles.dat);
    for n = 1:NFiles  
        %Import the file
        temp = importdata(poxFiles.dat{n}, DELIMITER, HEADERLINES);
        %interpolate data to every second
        tEnd = floor(temp.data(end,1));
        tElapsed = 0:tEnd;
        
        %create correct serial date numbers for all samples
        
        dateStr = temp.textdata{1,1}(16:26);
        S(n).date = dateStr;
        S(n).start = temp.textdata{3,1};
        S(n).elapsed = tElapsed;
        S(n).t = createDateNumFromElapsed(dateStr,S(n).start,tElapsed);
        S(n).headers = {'SpO2', 'PO Heartrate'};
        S(n).chUnits = {'% Saturation', 'BPM'};
        
        %determine unique values (some times are duplicates, this isn't a
        %problem because this is at a much faster rate than the 1s sampling
        %of the heartrates and spO2, as we will interpolate them anyway.
        [unTime, ia] = unique(temp.data(:,1));
        unHr = temp.data(ia,2);
        unSat = temp.data(ia,3);
        S(n).data(:,1) = interp1(unTime, unSat, tElapsed);  %interpolate to 1 sample per second
        S(n).data(:,2) = interp1(unTime, unHr, tElapsed);
        
       
        
    end
else %load in any mat files
    if strfind(poxFiles.mat{1}, 'pulseox') %if a dat file also exists then it will be the old script
        NFiles = length(poxFiles.mat);
        for n = 1:NFiles
            temp = load(poxFiles.mat{n});
            tempdata = [temp.elapsed2', temp.sat', temp.heartRate'];
            tempdata = interpolateToTimebase(tempdata, 1);  %interpolate to once per second
            startTimeVec = temp.timestamp2(1,:);
            startTimeStr = datestr(startTimeVec);
            
            S(n).date = startTimeStr(1:11);
            S(n).start = startTimeStr(13:end);
            S(n).elapsed = tempdata(:,1);                                  
            S(n).t = createDateNumFromElapsed(S(n).date, S(n).start, S(n).elapsed);
            S(n).data = tempdata(:,2:3);
            S(n).headers = {'SpO2', 'PO Heartrate'};
            S(n).chUnits = {'% Saturation', 'BPM'};
            

        end
    elseif strfind(poxFiles.mat{1}, 'po')
        NFiles = length(poxFiles.mat);
        for n = 1:NFiles
            temp = load(poxFiles.mat{n});
            tempdata = temp.poData(:,2:end);
            % turn -1's into NaNs
            noData = tempdata==-1;
            tempdata(noData) = NaN;
            
            tempdata = interpolateToTimebase(tempdata,1);
            startTimeStr = datestr(temp.poData(1,1));
            
            S(n).date = startTimeStr(1:11);
            S(n).start = startTimeStr(13:end);
            S(n).elapsed = tempdata(:,1);
            S(n).t = createDateNumFromElapsed(S(n).date, S(n).start, S(n).elapsed);
            S(n).data = tempdata(:,2:3);
            S(n).headers = {'SpO2', 'PO Heartrate'};
            S(n).chUnits = {'% Saturation', 'BPM'};
        end
    end
    
end

%Place all data into the same vector
poxData = arrangeTimeBlocks(S);
out.data = poxData(:,2:end);
out.t = poxData(:,1);
out.elapsed = (out.t - out.t(1))*60*60*24;
out.headers = S(1).headers;
out.chUnits = S(1).chUnits;




    
