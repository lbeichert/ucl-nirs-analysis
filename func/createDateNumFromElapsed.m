function t = createDateNumFromElapsed(startDate,startTime,tElapsed);
tStrStart = [startDate, ' ', startTime];
tStart = datenum(tStrStart);
for nt = 1:length(tElapsed)
    t(nt) = addtodate(tStart, tElapsed(nt), 'second');
end

