function S = loadSATrendFile(fileName)
%loads in a SA Instruments trend file
%Usage
% S = loadSATrendFile(fileName)
%arguments:
%   filename:   filename of the trend file (csv text)
%returns a structure S with fields:
%   data:     An m-by-n array of trend data with m samples and n trends
%   date:       the date the trend data begins on
%   start:      the start time of the trend data
%   elapsed:    a m element vector of seconds elapsed for each sample since
%               the start of the trend acquisition
%   t:          a m element vector of serial date numbers for each sample
%   headers:    n element cell array of column headers for the trends.

%load in the trend data
trendData = importTrendOutput(fileName);
trendInfo = trendData.textdata{1,1};

%% parse the header data to determine the start and stop times

trendInfo = regexprep(trendInfo, '"', ''); %remove " from the string if any exist

% Start Time
tkn = regexp(trendInfo, '(.+)Date, (.+)Output, (.+), Start, (.+), Stop, (.+)', 'tokens');
start = regexp(trendInfo, 'Start, (\d{2}:)?(\d{2}:)?(\d{2})', 'tokens');
start = [start{:}{:}];

% Stop Time
stop = regexp(trendInfo, 'Stop, (\d{2}:)?(\d{2}:)?(\d{2})', 'tokens');
stop = [stop{:}{:}];

% Start date
expdate = regexp(trendInfo, 'Date,\s(\d{2}[\/])?(\d{2}[\/])?(\d{2})', 'tokens');
expdate = [expdate{:}{:}];

%% Timing 
%create a vector a datenums based on the start time,
%and the elapsed seconds (first column of the trend data)
trendData.timeNum = createDateNumFromElapsed(expdate, start, trendData.data(:,1));


S.t = trendData.timeNum; %time of each sample as a serial date number
S.date = expdate;
S.start = start;
S.elapsed = trendData.data(:,1); %elapsed time in seconds from the start of the trend file
S.data = trendData.data(:,2:end); %trend data

%determine the column headers
colheaders = trendData.textdata{3};
colheaders = regexprep(colheaders, '"', ''); %strip "
colheaders = (textscan(colheaders, '%s', 'Delimiter', ','));
if all(size(colheaders) == [1,1])
    colheaders = colheaders{:};
end

S.headers = colheaders;



