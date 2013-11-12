function S = smoothData(dataStruct, smoothWindow, method)

if nargin < 3
    method = 'moving';
end
data = dataStruct.data;


span = smoothWindow / diff(dataStruct.elapsed(1:2));

for n = 1:length(dataStruct.headers)
    data(:,n) = smooth(data(:,n),span, method);
end

S = dataStruct;
S.data = data;
