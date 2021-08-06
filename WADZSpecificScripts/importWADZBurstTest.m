%The role of this function is to import the along-beam burst velocity
%records and the times (in datenums) corresponding to each timestep. Any
%filtering for bad data, spikes and excess tilt should be done within the
%function call, so that the output is a best approximation of the actual
%along-beam velocity time series.
function [burstDatenums,burstVelocityData] = importWADZBurstTest(burstCtr,tiltDataLong)

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
burstEnsembleNo = burstData(1:5:dataRecordLength,1);
burstYear = burstData(1:5:dataRecordLength,2);
burstMonth = burstData(1:5:dataRecordLength,3);
burstDay = burstData(1:5:dataRecordLength,4);
burstHour = burstData(1:5:dataRecordLength,5);
burstMinute = burstData(1:5:dataRecordLength,6);
burstSecond = burstData(1:5:dataRecordLength,7) + burstData(1:5:dataRecordLength,8)/100;
burstDatenums = datenum(burstYear,burstMonth,burstDay,burstHour,burstMinute,burstSecond);

burstVelocityData.beam1Raw = burstData(1:5:dataRecordLength,9:end);
burstVelocityData.beam2Raw = burstData(2:5:dataRecordLength,9:end);
burstVelocityData.beam3Raw = burstData(3:5:dataRecordLength,9:end);
burstVelocityData.beam4Raw= burstData(4:5:dataRecordLength,9:end);
burstVelocityData.beam5Raw = burstData(5:5:dataRecordLength,9:end);

%We filter the raw data for NaNs (or values of -32678, equal to -8000 in
%hex, which is the 'bad data' value in the RDI Workhorse firmware), and we
%normalise by 1000 since the raw data is in mms-1 rather than ms-1.
burstVelocityData.beam1Denanned = NaNFilterv2(burstVelocityData.beam1Raw); burstVelocityData.beam1Denanned = burstVelocityData.beam1Denanned/1000;
burstVelocityData.beam2Denanned = NaNFilterv2(burstVelocityData.beam2Raw); burstVelocityData.beam2Denanned = burstVelocityData.beam2Denanned/1000;
burstVelocityData.beam3Denanned = NaNFilterv2(burstVelocityData.beam3Raw); burstVelocityData.beam3Denanned = burstVelocityData.beam3Denanned/1000;
burstVelocityData.beam4Denanned = NaNFilterv2(burstVelocityData.beam4Raw); burstVelocityData.beam4Denanned = burstVelocityData.beam4Denanned/1000;
burstVelocityData.beam5Denanned = NaNFilterv2(burstVelocityData.beam5Raw); burstVelocityData.beam5Denanned = burstVelocityData.beam5Denanned/1000;

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

burstVelocityData.beam1 = excessTiltFilter(burstVelocityData.beam1Denanned,...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),2),...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),3));
burstVelocityData.beam2 = excessTiltFilter(burstVelocityData.beam2Denanned,...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),2),...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),3));
burstVelocityData.beam3 = excessTiltFilter(burstVelocityData.beam3Denanned,...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),2),...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),3));
burstVelocityData.beam4 = excessTiltFilter(burstVelocityData.beam4Denanned,...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),2),...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),3));
burstVelocityData.beam5 = excessTiltFilter(burstVelocityData.beam5Denanned,...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),2),...
                                            tiltDataLong(burstEnsembleNo(1):burstEnsembleNo(end),3));

%Function end
end