function [burstMeanDepth,burstSurfBin] = FoWDepthPreprocessing(depthParameters,burstRawDepth);

burstMeanDepth = mean(burstRawDepth,2);
burstSurfBin = 1 + floor((burstMeanDepth*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);%78;

end