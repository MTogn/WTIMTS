%For the WADZ north data set, preprocessing tasks are gathered into the
%function WADZPreprocessing_north, which includes:
% - Loading and formatting tilt and depth data
[tiltDataLong,demozoneRawDepth] = WADZPreprocessing_north();

%Initialising some key parameters
paramStruc.beamAngle = 25*pi/180; paramStruc.anisoParam = 0.1684;
%WADZ north data is collected with a Sentinel V50 500kHz ADCP, for which
%the inclined beam angle to the vertical is 25 degrees.
paramStruc.blankDist = 2.07; paramStruc.binVertSize = 0.6;
%Blanking distance is taken from IMARDIS record: https://portal.imardis.org/metadata/imardis:metadata:c0a6eadd-9c89-474e-8f7d-4052dc73bb29
%In this case "blanking distance" is the height above seabed + distance to first bin
paramStruc.sampFreq = 2;
paramStruc.dataLocation = 'C:\Users\michael\Documents\ADCP\NWDZ_north\burstData\';
%If mean flow properties must be calculated, a file containing heading data corresponding to each ensemble in each burst should also be
%present. If you already have burst-mean data for each burst, this can be 
paramStruc.rawHeadingDataFile = 'completeNWDZHeadingData_North.mat';