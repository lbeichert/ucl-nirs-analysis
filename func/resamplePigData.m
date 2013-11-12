function S = resamplePigData(pigData, pigletDir, subj, fileString, bNIRSTimesOnly)
% interpolate systemic and pulseox data to the same time increment as the nirs data (eg. every 60s)

tStart = [];
tEnd = [];
if nargin < 4
    fileString = [];
end

if ~isempty(pigData.usedsyst)
    syst = pigData.usedsyst;
    tStart = [tStart, syst.t(1)];
    tEnd = [tEnd, syst.t(end)];
    
systCols = length(syst.headers);
end

if ~isempty(pigData.pulseox)
    pox = pigData.pulseox;
    tStart = [tStart, pox.t(1)];
    tEnd = [tEnd, pox.t(end)];
    poxCols = length(pox.headers);
end

if ~isempty(pigData.nirs)
    nirs = pigData.nirs;
    tStart = [tStart, nirs.t(1)];
    tEnd = [tEnd, nirs.t(end)];
    nirsStart = nirs.t(1);
    nirsCols = length(nirs.headers);
end

nirsExposure = nirs.exposure;
dExp = nirsExposure/(60*60*24);

if bNIRSTimesOnly
    tStart = nirs.t(1);
    tEnd = nirs.t(end);
else
    tStart = min(tStart); %determine the start time of all the data
    tEnd = max(tEnd); %determine the end time of all the data
end


NExpBeforeNirsStart = (nirsStart - tStart)/dExp;
tStart = nirsStart - floor(NExpBeforeNirsStart)*dExp;
NExpAfterNirsStart = (tEnd - nirsStart)/dExp;
tEnd = nirsStart + floor(NExpAfterNirsStart)*dExp;

totalSeconds = (tEnd - tStart)*60*60*24;
starTimetStr = datestr(tStart);

elapsedNirs = round(linspace(0, totalSeconds, length(nirs.t)));
%elapsedNirs = 0:nirsExposure:totalSeconds;
expTimeNirs = createDateNumFromElapsed(starTimetStr(1:11), starTimetStr(13:end), elapsedNirs);



 out = [];
headers = {'time', 'elapsed (s)'};
units = {'s', 's'};
out = zeros(length(expTimeNirs), 2);
out(:,1) = expTimeNirs;
out(:,2) = elapsedNirs;
if ~isempty(pigData.usedsyst)
    % 1. interpolate systemic
    for n = 1:systCols
        out(:,2+n) = interp1(syst.t, syst.data(:,n), expTimeNirs);
        headers = [headers, syst.headers(n)];
        units = [units, syst.chUnits(n)];
    end
end


% 2. interpolate pulse ox data
if ~isempty(pigData.pulseox)
    poxInterp = ones(length(expTimeNirs), 2)*(-1);
    if isfield(pigData, 'pulseox')
        poxInterp(:,1) = interp1(pox.t, pox.data(:,1), expTimeNirs);
        poxInterp(:,2) = interp1(pox.t, pox.data(:,2), expTimeNirs);
        headers = [headers, pox.headers];
        units = [units, pox.chUnits];
    end
    out = [out, poxInterp];
end
% 3. pad out the nirs data

nirsOut = ones(length(expTimeNirs), nirsCols)*(-1);
nirsStartIndex = dsearchn(expTimeNirs', nirsStart);
nirsAcq = length(nirs.t);
nirsOut(nirsStartIndex:nirsStartIndex + nirsAcq-1, :) = nirs.results.Concentrations';
headers = [headers, nirs.headers];
units = [units, nirs.chUnits];


out = [out, nirsOut]; %concatenate

S.data = out(:,3:end);
S.t = out(:,1);
S.elapsed = out(:,2);
S.headers = headers(3:end);
S.units = units(3:end);
S.subj = subj;

save([pigletDir, filesep, 'NTB_', subj, fileString, '.mat'], 'S');

fname = [pigletDir, filesep, 'NTB_', subj, fileString, '.txt'];
fid = fopen(fname, 'wt');

for n = 1:length(headers)
    fprintf(fid, '%s\t', headers{n});
end
fprintf(fid, '\n');
fclose(fid);

dlmwrite(fname, out, '-append', 'delimiter', '\t', 'precision', '%.6f');
