function S = phosDataFillGap(pData)
%fills in the gap between the baseline data and insult date with regular
%timepoints 

t31P = pData.t;
%fill in the gaps
t31PInsultStart = t31P(11);
dt31P = t31P(2) - t31P(1); %time between 31P acquisitions
t31PNew = t31P(1:10);
t31PGap = (t31P(10)+dt31P:dt31P:t31P(11))';
t31PNew = [t31PNew; t31PGap;t31P(11:end)];
elInsultStart = dsearchn(t31PNew,t31PInsultStart);
N31P = length(t31PNew);
for n=1:N31P
    elapsed31P(n) = etime(datevec(t31PNew(n)), datevec(t31PNew(1)));
end

data = nan(N31P,3);
data(1:10,:) = pData.data(1:10,:);
data(elInsultStart:end,:) = pData.data(11:end,:);
S = pData;
S.headers = pData.headers;
S.data = data;
S.date = pData.date;
S.start = pData.start;
S.elapsed = elapsed31P';
S.t = t31PNew;

