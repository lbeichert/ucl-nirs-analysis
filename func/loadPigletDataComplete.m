function pigData = loadPigletDataComplete(pigletDir, nirsStartTime, nirsExposure, expDate, expEnd)
% Loads in piglet data (nirs and systemic)

%look for tarragona nirs files
pigData.files.nirs = findTarragonaFiles(pigletDir);

%look for pulseox files
pigData.files.pulseox = findPulseOxFiles(pigletDir);

%Look for SA Instrument Files
pigData.files.trend = findSATrendFiles(pigletDir);

% Look for NI systemic data
pigData.files.syst = findNISystemicFiles(pigletDir);

nirsStartTimeFull = datenum(expDate, 'dd/mm/yyyy') + nirsStartTime;
expEnd = datenum(expDate, 'dd/mm/yyyy') + expEnd;
%% Load In Data
% 1. Systemic

if pigData.files.trend.present
    a=0;
    trendCols = [4,5,3];
    for n = 1:length(pigData.files.trend.files)
        for m = 1:length(pigData.files.trend.files{n})
            a=a+1;
            trendFile = pigData.files.trend.files{n}{m};
            trends(a) = loadSATrendFile(trendFile);
        end
    end
    trendData = arrangeTimeBlocks(trends); %synchronises all of the trend data and places in one large array
    
    %only use the useful trend data
    usedData = trendData(:,[1,trendCols]); %first 7 columns
    
    pigData.syst = trends;
    pigData.usedsyst.data = usedData(:,2:end);
    pigData.usedsyst.t = usedData(:,1);
    pigData.usedsyst.elapsed = (pigData.usedsyst.t - pigData.usedsyst.t(1))*60*60*24;
    pigData.usedsyst.headers = trends(1).headers(trendCols);
end

if pigData.files.syst.present
    
    for n = 1:length(pigData.files.syst.mat)
        syst(n) = importNISystemic(pigData.files.syst.mat{n});
        %determine NIRS start times from the NIRS gate channel
        nirsGateCh = (strcmp(syst(n).headers,'NIRS Gate'));
        nirsGate{n} = round(syst(n).data(:,nirsGateCh));
        nirsGateTimes{n} = findGateTimes(nirsGate{n}, syst(n).t); %determine the start times of the NIRS acquisition
        
    end
    
    systData = arrangeTimeBlocks(syst); %synchronise all physiological data if multiple files exist
    
    pigData.syst = syst;
    
    %only use useful trend data (BP Mean, BP Rate, Temperature, BP
    %Systolic - first four channels)
    systCols = [2,3,4,5];
    usedData = systData(:,1:5);
    pigData.usedsyst.data = usedData(:,2:end);
    pigData.usedsyst.t = usedData(:,1);
    pigData.usedsyst.elapsed = (pigData.usedsyst.t - pigData.usedsyst.t(1))*60*60*24;
    pigData.usedsyst.headers = syst(1).headers(systCols-1);
    
    %consolidate nirs start times
    nirsGateTimes = [nirsGateTimes{:}];
    %find the closest gate time that corresponds to nirsStarTime (if there
    %are multiple gate times)
    nirsStartTimeFull = nirsGateTimes(dsearchn(nirsGateTimes, nirsStartTimeFull));
    
end






%% 2. Pulse Oximetry
pigData.pulseox = [];
if pigData.files.pulseox.present
    pigData.pulseox = importPulseOx(pigData.files.pulseox);
    
end

%% 3. NIRS
%determine ref files
refFiles = stringToken(pigData.files.nirs, 'ref(.+)');
if length(refFiles) > 1
    optionsText = [pigletDir,'\nPlease select the reference spectrum and press ENTER\n'];
    
    for i = 1:length(refFiles)
        optionsText = [optionsText, '(',num2str(i), ')', refFiles{i}, '\n'];
    end
    
    bRefSelect = false;
    while(~bRefSelect)
        refIndex = str2num(input(optionsText, 's'));
        if ~isempty(refIndex)
            bRefSelect = true;
        end
    end
else
    refIndex = 1;
end
%determine nirs data Files
nirsFiles = stringToken(pigData.files.nirs, '^((?!ref).)*$');
if length(nirsFiles) > 1
    optionsText = [pigletDir,'\nPlease select the data files (numerically, separated by commas if more than 1) and press ENTER\n'];
    for i = 1:length(nirsFiles)
        optionsText = [optionsText, '(',num2str(i), ')', nirsFiles{i}, '\n'];
    end
    bNIRSSelect = false;
    while(~bNIRSSelect)
        nirsIndex = str2num(input(optionsText, 's'));
        if ~isempty(nirsIndex)
            bNIRSSelect = true;
        end
    end
else
    nirsIndex = 1;
end

for n = 1:length(nirsIndex)
    bStartTime = true;
    if nirsStartTime == 0
        bStartTime = false;
        
    end
    
    nirs(n) = importTarragona(refFiles{refIndex}, nirsFiles{nirsIndex(n)}, nirsStartTimeFull, nirsExposure, expEnd, bStartTime);
end

% %read in tarragona data
% lambdaMin = 780;
% lambdaMax = 900;
% bAutoPath = true;
% nNirsRef = 5; %number of reference scans prior to acquiring NIRS spectra
%
%
% for n = 1:length(nirsIndex)
%     %read intensities
%     refInts{n} = tarragona_intensityread(refFiles{refIndex}, 0, 1, 2);
%     nirsInts{n} = tarragona_intensityread(nirsFiles{nirsIndex(n)}, 0, 1, 2);
%
%     %process data
%     [nirs(n).raw, nirs(n).results, nirs(n).resultsFile] = Broadband_Tarragona(refFiles{refIndex},1, refInts{n}.num_spectra, nirsFiles{nirsIndex(n)},1, nirsInts{n}.num_spectra, [lambdaMin, lambdaMax], bAutoPath, 0);
%
%     if nirsStartTime == 0
%         %work back from the experiment end time
%         for i = length(nirs(n).results.Time):-1:1
%             nirs(n).results.Time(i) = addtodate(usedData(end,1), -nirs.results.Time(i)*nirsExposure, 'second');
%         end
%         nirs(n).results.Time = fliplr(nirs(n).results.Time);
%     else
%         %nirs start time is given, so use that
%         for i = 1:length(nirs(n).results.Time)
%             nirs(n).results.Time(i) = addtodate(datenum(trends(1).date) + nirsStartTime, (nirs.results.Time(i)+nNirsRef)*nirsExposure, 'second');
%         end
%     end
% end

pigData.nirs = nirs;

