function [tstart, tfinish] = findGateTimes(gate,t)
%determines the positive and negative transition points of digital gating data, and returns
%the times that these occur (start, finish times)

gatediff = round(diff(gate)); %calculate the differential, and round so there is only one value per transition

pos = gatediff > 0;
neg = gatediff < 0;

tstart = t(pos);
tfinish = t(neg);