%For the original WADZ data set, preprocessing tasks are gathered into the
%function WADZPreprocessing, which includes:
% - Loading and formatting tilt and depth data
[tiltDataLong,demozoneRawDepth] = WADZPreprocessing();

%Initialise some key parameters
paramStruc.beamAngle = 25*pi/180; paramStruc.anisoParam = 0.1684;
paramStruc.blankDist = 1.89; paramStruc.binVertSize = 0.6;
paramStruc.sampFreq = 2;
paramStruc.dataLocation = 'C:\Users\michael\Documents\ADCP\DEMOZONE\';