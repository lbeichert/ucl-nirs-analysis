function gateTime = getNirsGate(syst, nirsStartTimeFull)

nirsGateCh = (strcmp(syst.headers,'NIRS Gate'));
nirsGate = round(syst.data(:,nirsGateCh));
nirsGateTimes = findGateTimes(nirsGate, syst.t);


gateTime = nirsGateTimes(dsearchn(nirsGateTimes, nirsStartTimeFull));