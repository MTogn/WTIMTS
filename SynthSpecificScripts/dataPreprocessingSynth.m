%Initialise some key parameters
paramStruc.beamAngle = 25*pi/180; paramStruc.anisoParam = 0.1684;
paramStruc.blankDist = 1.6; paramStruc.binVertSize = 1;
paramStruc.sampFreq = 2;
paramStruc.dataLocation = 'C:\Users\michael\Documents\WTIMTS\virtADCP\synthADCPBursts22May23_irreg\';

%This loads the burstMeanDepths, burstMaxBins and wholeRecordDatenums
%into the workspace from the file in which they're stored
load([paramStruc.dataLocation 'synthRecordData.mat']);