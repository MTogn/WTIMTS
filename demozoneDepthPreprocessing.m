function [burstDepths,burstMaxBins] = demozoneDepthPreprocessing(depthParameters,burstEnsNoLims,burstStartIndex,burstEndIndex);

demozoneRawDepth = load('C:\Users\michael\Documents\ADCP\DEMOZONE\demozoneRawDepth');
demozoneRawDepth = demozoneRawDepth(:,[1 9]);

burstDepths = nan(burstEndIndex,1);
burstMaxBins = burstDepths;

for burstCtr = burstStartIndex:burstEndIndex
   burstStartEns = burstEnsNoLims(burstCtr,1); burstEndEns = burstEnsNoLims(burstCtr,2);
   burstDepths(burstCtr) = mean(demozoneRawDepth(burstStartEns:burstEndEns,2));
   burstMaxBins(burstCtr) = 1 + floor((burstDepths(burstCtr)*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);
end