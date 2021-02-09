%The role of this function is to import the along-beam burst velocity
%records and the times (in datenums) corresponding to each timestep. Any
%filtering for bad data, spikes and excess tilt should be done within the
%function call, so that the output is a best approximation of the actual
%along-beam velocity time series.
function [burstEnsembleNos,burstDatenums,burstVelocityData] = importWADZBurst(burstCtr,tiltDataLong)

WADZDataAbsoluteLocation = 'C:\Users\michael\Documents\ADCP\DEMOZONE\burstData\';
burstData = load([WADZDataAbsoluteLocation 'burstData' int2str(burstCtr)]);
burstData = burstData.burstData;

% - For the WADZ data, the first column contains the ensemble number and
%the next 7 contain the time data.
% - Rows occur in groups of five; each group of five rows corresponds to a
%single record/time.
% - Within a group of five rows, row 1 contains beam 1 data, row 2 contains
%beam 2 data etc.

dataRecordLength = size(burstData,1);
burstEnsembleNos = burstData(1:5:dataRecordLength,1);
burstYear = burstData(1:5:dataRecordLength,2);
burstMonth = burstData(1:5:dataRecordLength,3);
burstDay = burstData(1:5:dataRecordLength,4);
burstHour = burstData(1:5:dataRecordLength,5);
burstMinute = burstData(1:5:dataRecordLength,6);
burstSecond = burstData(1:5:dataRecordLength,7) + burstData(1:5:dataRecordLength,8)/100;
burstDatenums = datenum(burstYear,burstMonth,burstDay,burstHour,burstMinute,burstSecond);

burstVelocityData.beam1 = burstData(1:5:dataRecordLength,9:end);
burstVelocityData.beam2 = burstData(2:5:dataRecordLength,9:end);
burstVelocityData.beam3 = burstData(3:5:dataRecordLength,9:end);
burstVelocityData.beam4 = burstData(4:5:dataRecordLength,9:end);
burstVelocityData.beam5 = burstData(5:5:dataRecordLength,9:end);

%We filter the raw data for NaNs (or values of -32678, equal to -8000 in
%hex, which is the 'bad data' value in the RDI Workhorse firmware), and we
%normalise by 1000 since the raw data is in mms-1 rather than ms-1.
burstVelocityData.beam1 = NaNFilterv2(burstVelocityData.beam1); burstVelocityData.beam1 = burstVelocityData.beam1/1000;
burstVelocityData.beam2 = NaNFilterv2(burstVelocityData.beam2); burstVelocityData.beam2 = burstVelocityData.beam2/1000;
burstVelocityData.beam3 = NaNFilterv2(burstVelocityData.beam3); burstVelocityData.beam3 = burstVelocityData.beam3/1000;
burstVelocityData.beam4 = NaNFilterv2(burstVelocityData.beam4); burstVelocityData.beam4 = burstVelocityData.beam4/1000;
burstVelocityData.beam5 = NaNFilterv2(burstVelocityData.beam5); burstVelocityData.beam5 = burstVelocityData.beam5/1000;

%We must also filter for times of excess tilt. For the WADZ data, tilt data
%is stored in a file in the immediate parent directory of the directory
%containing the burst data. I think it is likely better to keep the
%whole-record tilt data in the workspace rather than load the file
%containing the data with each call to this function, but this may need to
%be revised.
%The tilt data replicates the ensemble numbers and time/date
%data included in the burst files containing the velocity records; only one
%of these is needed to find the portion of the whole-deployment tilt data
%corresponding to each shorter-duration burst, so we discard the time/date
if ~exist('tiltDataLong')
    tiltDataLong = load('C:\Users\michael\Documents\ADCP\DEMOZONE\completeNWDZTiltData');
    tiltDataLong = tiltDataLong(:,[1 9:10]);
end

burstVelocityData.beam1 = excessTiltFilter(burstVelocityData.beam1,...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),2),...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),3));
burstVelocityData.beam2 = excessTiltFilter(burstVelocityData.beam2,...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),2),...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),3));
burstVelocityData.beam3 = excessTiltFilter(burstVelocityData.beam3,...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),2),...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),3));
burstVelocityData.beam4 = excessTiltFilter(burstVelocityData.beam4,...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),2),...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),3));
burstVelocityData.beam5 = excessTiltFilter(burstVelocityData.beam5,...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),2),...
                                            tiltDataLong(burstEnsembleNos(1):burstEnsembleNos(end),3));

%Function end
end