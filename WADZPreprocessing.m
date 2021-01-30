%Burst velocities at times of excess tilt should be disregarded. For the
%WADZ data, tilt data is stored in a file in the immediate parent directory
%of the directory containing the burst data. This is too large a file to
%load and unload for each burst, so we load and preprocess the tilt data
%here to strip it down to as small a size as possible. The data replicates
%the ensemble numbers and time/date data included in the burst files
%containing the velocity records; only one of these is needed to find the
%portion of the whole-deployment tilt data corresponding to each shorter-
%duration burst, so we discard the time/date.
function tiltDataLong = WADZPreprocessing();

    if ~exist('tiltDataLong')
        tiltDataLong = load('C:\Users\michael\Documents\ADCP\DEMOZONE\completeNWDZTiltData');
        tiltDataLong = tiltDataLong(:,[1 9:10]);
    end

end