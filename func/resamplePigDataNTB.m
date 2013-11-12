function S = resamplePigDataNTB(NTBData,pData)
%resamples the piglet NIRS date to the phosphorus data times


t = NTBData.t;
elapsed = NTBData.elapsed;

nCol = length(pData.headers);

N = length(t);

data = zeros(N,nCol);
for cc = 1:nCol
    data(:,cc) = interp1(pData.t, pData.data(:,cc), t);
end


%combine all data
S = NTBData;
S.data = [NTBData.data,data];
S.headers = [NTBData.headers, pData.headers];
S.units = [NTBData.units, pData.units];

