function out = arrangeTimeBlocks(S)
%creates a single dataset from the discrete data sets and places them all
%on the same timebase
%sample interval must be equal for all blocks

NBlocks = length(S); %number of data blocks
ts=[];

for n=1:NBlocks
    if isrow(S(n).t)
        ts = [ts,S(n).t]; %concatenate times
    else
        ts = [ts, S(n).t'];
    end
end
startTime = min(ts);
endTime = max(ts);
interval = ts(2) - ts(1);

totalTime = startTime:interval:endTime;
length(totalTime);

nData = size(S(1).data,2); %number of data columns
out = zeros(length(totalTime), nData+1);
out(:,1) = totalTime;

for n = 1:NBlocks
    %determine start time of the block
    blockStart = S(n).t(1);
    blockLength = length(S(n).t);
    
    blockStartIndex = find(totalTime>blockStart - interval/2 & totalTime < blockStart + interval/2, 1,'first');
    blockIndex = blockStartIndex:blockStartIndex+blockLength-1;
    %length(blockIndex);
    %length(S(n).data);
    out(blockIndex,2:end) = S(n).data;
end

fixindex = find(out(:,1)==0);
if ~isempty(fixindex)
    for n = 1:fixindex
        out(fixindex,1) = out(fixindex-1,1) + interval;
    end
end
    


