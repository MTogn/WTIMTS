function [burstDepth,burstMaxBin] = demozoneDepthPreprocessing_north(depthParameters,demozoneRawDepth,burstEnsNoLims,burstCtr);

burstStartEns = burstEnsNoLims(1); burstEndEns = burstEnsNoLims(2);
burstDepth = mean(demozoneRawDepth(burstStartEns:burstEndEns,2));
burstMaxBin = 1 + floor((burstDepth*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);

end