function S = denoiseData(dataStruct, channel, noiseThresh, smoothWindow)

if isempty(channel)
    channel = dataStruct.headers;
end

if ~iscell(channel)
    tmp = channel;
    clear channel;
    channel{1} = tmp;
end

span = smoothWindow / diff(dataStruct.elapsed(1:2));
for n = 1:length(channel)
    
    [temp, dc] = getDataTrace(dataStruct, channel{n});
    tempsm = smooth(temp, span);
    
    %calculate del2
    noiseEst = del2(tempsm, dataStruct.elapsed);
   
    
    stdnoise = nanstd(noiseEst);
    
    D = abs(noiseEst) < noiseThresh*stdnoise;
    

    
    temp2 = interp1(dataStruct.elapsed(D), temp(D), dataStruct.elapsed);
    
    dataStruct.data(:,dc) = temp2;
    
end

S = dataStruct;