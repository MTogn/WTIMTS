function [burstMeanDepth,burstMaxBin] = partracDepthPreprocessing(depthParameters,burstRawDepth);

burstMeanDepth = mean(burstRawDepth,2);
burstMaxBin = 1 + floor((burstMeanDepth*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);

end