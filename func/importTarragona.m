function S = importTarragona(refFile, nirsFile, startTime, exposure, expEnd, bStartTime)
%loads in tarragona NIRS files and performs basic analysis 
%Usage:
%   S = importTarragona(refFile, nirsFile, startTime, exposure)
%arguments
%   refFile:    reference spectra filename
%   nirsFile:   nirs acquisition filename
%   startTime:  the time the nirs acquisition is initiated as a serial
%               date number
%   exposure:   the exposure length in seconds

%read in tarragona data
lambdaMin = 780;
lambdaMax = 900;
bAutoPath = true; 
nNirsRef = 5; %number of reference scans prior to acquiring NIRS spectra

% construct filename to save data to
[pth, name, ~] = fileparts(nirsFile); %determine filename
savefile = [pth, filesep, name, '.mat'];
%S = [];

if exist(savefile, 'file') %if the processed file exists load it in
    temp = load(savefile);
    S = temp.S;
else %otherwise call broadband tarragona
    
    refInts = tarragona_intensityread(refFile, 0,1,2);
    nirsInts = tarragona_intensityread(nirsFile, 0,1,2);
    [raw results resultsFile] = Broadband_Tarragona(refFile,1,refInts.num_spectra,nirsFile,1,nirsInts.num_spectra, [lambdaMin, lambdaMax], bAutoPath, 0,false,false);
                               
    %create times
    elapsed = (1:nirsInts.num_spectra)*exposure;
    
    if ~bStartTime
        firstAcqTime = expEnd - elapsed(end)/(60*60*24);
    else        
        firstAcqTime = addtodate(startTime, nNirsRef*exposure, 'second'); %time for the first acquisition (start time + reference time, usually 5 exposures)
    end
    firstAcqTimeStr = datestr(firstAcqTime);
    
    S.date = firstAcqTimeStr(1:11);
    S.start = firstAcqTimeStr(13:end);
    S.elapsed = elapsed;
    S.t = createDateNumFromElapsed(S.date, S.start, S.elapsed);
    S.raw = raw;
    S.results = results;
    S.exposure = exposure;
    S.headers = {'Hb','HbO2','CtOx','HbT','HbDiff'};
    S.chUnits = {'umol/L', 'umol/L', 'umol/L', 'umol/L', 'umol/L'};
    save(savefile, 'S');
    
end
    
    

