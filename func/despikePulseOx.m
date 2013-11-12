function S = despikePulseOx(pigData)
%interpolates between points where pulse oximetry data is not present, used
%to despike where there are '-1' values.  Only valid for short periods
%where there is no information.

pox = pigData.pulseox;

nPoxData = size(pox.data,2);
stdTuning = 0;

for n=1:nPoxData
    data = pox.data(:,n);
    D = data > 10;
    data = interp1(pox.t(D), data(D), pox.t); %interpolate at all times, using only the good data as input
    
%     D = abs(del2(data)) >= rms(del2(data)) + stdTuning*std(del2(data));
%     D = ~D;
%     data = interp1(pox.t(D), data(D), pox.t);
    pox.data(:,n) = data;
end

pigData.pulseox = pox;
S = pigData;