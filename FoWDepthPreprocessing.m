function [burstMeanDepth,burstMaxBin] = FoWDepthPreprocessing(depthParameters,burstRawDepth);

burstMeanDepth = mean(burstRawDepth,2);
burstMaxBin = 78;%1 + floor((burstMeanDepth*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);

end