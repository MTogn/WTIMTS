function [waveBeamVelocity,nonwaveBeamVelocity] = spectralFilter(beamVelocity,sampFreq,filterFunction)

waveBeamVelocity = nan(size(beamVelocity)); nonwaveBeamVelocity = waveBeamVelocity;
for depthCtr = 1:size(beamVelocity,2)
    [tempWSST,tempFreqVec] = wsst(beamVelocity(:,depthCtr),sampFreq);
    %Call filter function here
    [waveTransformComponent,nonwaveTransformComponent] = filterFunction(tempWSST,tempFreqVec);
    waveBeamVelocity(:,depthCtr) = iwsst(waveTransformComponent);
    nonwaveBeamVelocity(:,depthCtr) = iwsst(nonwaveTransformComponent);
end

end