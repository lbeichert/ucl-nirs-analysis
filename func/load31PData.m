function S = load31PData(insultAnalysisExcelFile)

[status, sheets, ~] = xlsfinfo(insultAnalysisExcelFile);

if isempty(status)
    error('not a valid excel spreadsheet');
    return
end

NSheets = length(sheets);

[~,txt,raw] = xlsread(insultAnalysisExcelFile,1);

hdr = {txt{6,:}};

%determine columns for the following

dateCol = find(strcmp(hdr,'date')); %date
scanCol = find(strcmp(hdr, 'scan')); %scan
dayCol = find(strcmp(hdr,'day')); %day
timeCol = find(strcmp(hdr,'time')); %time
PCreppCol = find(strcmp(hdr, 'PCr/epp')); %PCr/epp
NTPeppCol = find(strcmp(hdr, 'NTP/epp')); %NTP/epp
PieppCol = find(strcmp(hdr, 'Pi/epp')); %Pi/epp
deoccCol = find(strcmp(hdr, 'StartDEOCCL')); %de-occlusion time

scans ={raw{:,scanCol}};
range = 7:find(strcmp(scans, 'day2'))-1;
N = length(range);

%sort out times
dates = num2str([raw{range, dateCol}]');
dateformat = 'yyyymmdd';
times = num2str([raw{range, timeCol}]');
timeformat = 'hhmmss';

dnum = datenum(dates, dateformat);
dstr = datestr(dnum);
h=1:2;
m=3:4;
s=5:6;
cl = repmat(':', N,1);
times2 = [times(:,h), cl, times(:,m), cl, times(:,s)];
fullDateTime = [dstr, repmat(' ',N,1), times2];
% ft = fractionalTime([dstr, repmat(' ',N,1), times2]);
% calculate fractional time by hand
for j=1:size(times, 1)
    ft(j) = (str2double(times(j,h))*3600 + str2double(times(j,m))*60 + str2double(times(j,s)))/(24*60*60);
end

S.t = dnum + ft'; %serial date number
S.elapsed = zeros(N,1);
for n = 1:N
    S.elapsed(n) = etime(datevec(S.t(n)), datevec(S.t(1)));
end

deocc = [raw{range, deoccCol}]';
PCr = [raw{range, PCreppCol}]';
NTP = [raw{range, NTPeppCol(1)}]';
Pi = [raw{range, PieppCol}]';

%insult end time
insultEndEl = dsearchn(deocc, 0);
insultEnd = S.t(insultEndEl);

S.headers = {'NTP/epp', 'PCr/epp', 'Pi/epp'};
S.data = [NTP, PCr, Pi];
S.date = fullDateTime(1,1:11);
S.start = fullDateTime(1,13:end);
S.units = {'', '', ''};
S.insultEnd = insultEnd;
