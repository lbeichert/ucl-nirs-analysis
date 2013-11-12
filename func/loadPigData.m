function pigData = loadPigData(refFile, NIRSFile, SystFile, PoxFile, nirsStartTime, nirsExposure, expDate, systHeaders)
%loads in systemic and NIRS piglet data
%arguments
% refFile       filename of the NIRS reference file
% NIRSFile      filename of the NIRS data file
% SystFile      filename of the systemic data file (NI data acquisition
% only)
% PoxFile       filename of the pulse oximetry file
% nirsStartTime start time of the NIRS acquisition as a decimal fraction of
% 24 hours
% nirsExposure  NIRS exposure time in seconds
% expDate       experiment date in dd/mm/yyyy format
% systHeaders   cell array of systemic channel headers to copy to the usedData structure



nirsStartTimeFull = datenum(expDate, 'dd/mm/yyyy') + nirsStartTime;


%% Systemic
if ~isnan(SystFile)
    syst = importNISystemic(SystFile);
    nirsGateTime = getNirsGate(syst, nirsStartTimeFull);
    syst = smoothData(syst, 10, 'sgolay'); %smooth data with a 10s window

    %determine used column headers

    for n = 1:length(systHeaders)
        systCols(n) = find(strcmp(syst.headers, systHeaders{n}));
    end

    usedsyst = syst;
    usedsyst.data = syst.data(:,systCols);
    usedsyst.headers = systHeaders;
    usedsyst.chUnits = syst.chUnits(systCols);

    pigData.syst = syst;
    pigData.usedsyst = usedsyst;
else
    pigData.syst = [];
    pigData.usedsyst = [];
    nirsGateTime = nirsStartTimeFull;
end
%% Pulse oximetry
if ~isnan(PoxFile)
    pigData.pulseox = importPulseOx(PoxFile);
else
    pigData.pulseox = [];
end
    


%% NIRS
if ~isnan(NIRSFile)
    pigData.nirs = importTarragona(refFile, NIRSFile, nirsGateTime, nirsExposure, [], true);
else
    pigData.nirs=[];
end
