function out = interpolateToTimebase(data,timebase)
%data is an array, first column is times, second onwards is data
%timebase and first column must be in the same units

nCols = size(data,2) - 1;
tStart = floor(data(1,1));
tEnd = floor(data(end,1));

tElapsed = tStart:timebase:tEnd;

[unTime, ia] = unique(data(:,1)); %determine all the unique times (otherwise the interpolation won't work)

for n = 1:nCols
    unData(:,n) = data(ia,n+1);
end

for n = 1:nCols
  intData(:,n) = interp1(unTime, unData(:,n), tElapsed);
end

out= [tElapsed', intData];





