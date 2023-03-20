%Convert virtADCP data into the format needed for ADCP processing
load('C:\Users\michael\Documents\WTIMTS\virtADCP\output08Mar23.mat')
paramStruc.beamAngle = deg2rad(ADCPData.beamAngle);
paramStruc.anisoParam = 0.1684;
paramStruc.blankDist = ADCPData.blankDistance;
paramStruc.binVertSize = ADCPData.dz;
paramStruc.sampFreq = Wave.sampFreq;

burstBeamVelocities.beam1 = flowData.beam1.beamVel';
burstBeamVelocities.beam2 = flowData.beam2.beamVel';
burstBeamVelocities.beam3 = flowData.beam3.beamVel';
burstBeamVelocities.beam4 = flowData.beam4.beamVel';

% addpath('C:\Users\michael\Documents\WTIMTS\analysisCode');
testSynthADCPTKE = calcBurst4BeamTKE(burstBeamVelocities,paramStruc);
testSynthSpecFilterStoppedTKE = calcBurst4BeamTKE(specFilterStoppedVelocities,paramStruc);
testSynthSpecFilterPassedTKE = calcBurst4BeamTKE(specFilterPassedVelocities,paramStruc);