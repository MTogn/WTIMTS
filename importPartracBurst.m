%The role of this function is to import the along-beam burst velocity
%records and the times (in datenums) corresponding to each timestep. Any
%filtering for bad data, spikes and excess tilt should be done within the
%function call, so that the output is a best approximation of the actual
%along-beam velocity time series.
function [burstEnsembleNos,burstDatenums,burstBeamVelocities,burstDepth] = importPartracBurst(dataLocation,burstCtr)

%Simply loading all the data from the saved Partrac bursts should create
%the output variables.
load([dataLocation 'burstData' int2str(burstCtr)]);

%We filter the raw data for NaNs (or values of -32678, equal to -8000 in
%hex, which is the 'bad data' value in the RDI Workhorse firmware), and we
%normalise by 1000 since the raw data is in mms-1 rather than ms-1.
burstBeamVelocities.beam1 = NaNFilterv2(burstBeamVelocities.beam1);
burstBeamVelocities.beam2 = NaNFilterv2(burstBeamVelocities.beam2);
burstBeamVelocities.beam3 = NaNFilterv2(burstBeamVelocities.beam3);
burstBeamVelocities.beam4 = NaNFilterv2(burstBeamVelocities.beam4);

end