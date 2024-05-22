function [burstDepth,burstMaxBin] = demozoneDepthPreprocessing_north(depthParameters,demozoneRawDepth,burstEnsNoLims,burstCtr);

burstStartEns = burstEnsNoLims(1); burstEndEns = burstEnsNoLims(2);
if size(demozoneRawDepth,1) ~= 2 & size(demozoneRawDepth,2) == 2,
    demozoneRawDepth = demozoneRawDepth';
end
burstDepth = mean(demozoneRawDepth(2,burstStartEns:burstEndEns));
burstMaxBin = 1 + floor((burstDepth*cos(depthParameters.beamAngle) - depthParameters.blankDist)/depthParameters.binVertSize);

end