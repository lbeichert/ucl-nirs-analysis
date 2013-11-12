function S = resamplePigDataPTB(NTBData,pData)
%resamples the piglet NIRS date to the phosphorus data times


t = pData.t;
elapsed = pData.elapsed;

nCol = length(NTBData.headers);

N = length(t);

data = zeros(N,nCol);
for cc = 1:nCol
    data(:,cc) = interp1(NTBData.t, NTBData.data(:,cc), t);
end


%combine all data
S = pData;
S.data = [data, pData.data];
S.headers = [NTBData.headers, pData.headers];
S.units = [NTBData.units, pData.units];


    