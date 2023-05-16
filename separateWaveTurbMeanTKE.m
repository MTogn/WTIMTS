%This function finds the mean of a whole-record TKE array, split into the
%part of the mean that is attributable to turbulence only and the part that
%is attributable to wave action. The criterion for distinguishing between
%the two is when the EC of the first EOF is <0.
function [meanTKEWave,meanTKETurb] = separateWaveTurbMeanTKE(TKE,EC_1stEOF);

meanTKETotal = mean(TKE,1,'omitnan');
meanTKETurb = mean(TKE(EC_1stEOF < 0,:),1,'omitnan');
meanTKEWave = meanTKETotal - meanTKETurb;

end