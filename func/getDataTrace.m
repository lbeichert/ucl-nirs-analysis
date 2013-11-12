function [out, dc] = getDataTrace(dataStruct, headerString)

dc =find(strcmp(dataStruct.headers, headerString));

out = dataStruct.data(:,dc);