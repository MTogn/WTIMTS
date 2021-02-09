function [burstDepths,burstMaxBins] = demozoneDepthPreprocessing(depthParameters,burstEnsNoLims,burstStartIndex,burstEndIndex);

demozoneRawDepth = load('C:\Users\michael\Documents\ADCP\DEMOZONE\demozoneRawDepth')

burstDepths = nan(burstEndIndex - burstStartIndex + 1,:);
burstMaxBins = burstDepths;

for burstCtr = burstStartIndex:burstEndIndex
    

end