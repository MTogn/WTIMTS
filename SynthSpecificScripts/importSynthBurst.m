%The role of this function is to import the along-beam burst velocity
%records and the times (in datenums) corresponding to each timestep.
function burstBeamVelocities = importSynthBurst(dataLocation,burstCtr)

%Simply loading all the data from the saved bursts should create the output
%variables.
load([dataLocation 'burst' int2str(burstCtr)]);

end