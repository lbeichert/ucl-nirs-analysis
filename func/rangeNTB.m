function S = rangeNTB(NTBData, tStart, tEnd, methodStr)


elapsed = NTBData.elapsed;
t = NTBData.t;

if strcmp(methodStr, 'elapsed')
     startEl = dsearchn(elapsed, tStart)
    endEl = dsearchn(elapsed, tEnd)
end
if strcmp(methodStr, 'time')
    startEl = dsearchn(t, tStart);
    endEl = dsearchn(t, tEnd);
end


range = startEl:endEl;
S = NTBData;
S.t = NTBData.t(range);
S.data = NTBData.data(range,:);
S.elapsed = elapsed(range) - elapsed(startEl);

%% Re-normalise the NIRS data

nirsHeaders = {'Hb','HbO2','CtOx','HbT','HbDiff'};
nirsChan = false(size(S.headers));
%determine nirs channel numbers
for n = 1:length(nirsHeaders)
    nirsChan = nirsChan | strcmp(S.headers, nirsHeaders{n}); 
end
normFactor =  nirsChan.*S.data(1,:); %create normalisation factors for the nirs channels
S.data = S.data - repmat(normFactor,length(S.elapsed),1); %normalise
