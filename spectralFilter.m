function [waveBeamVelocity,nonwaveBeamVelocity] = spectralFilter(beamVelocity,sampFreq,filterFunction,varargin)

waveBeamVelocity = nan(size(beamVelocity)); nonwaveBeamVelocity = waveBeamVelocity;
for depthCtr = 1:size(beamVelocity,2)
    [tempWSST,tempFreqVec] = wsst(beamVelocity(:,depthCtr),sampFreq);
%Call filter function here
%If the optional argument in is supplied (i.e., nargin > 3), then it should
%be a structure containing key parameters to be passed to the filter.
    if nargin > 3,
        [waveTransformComponent,nonwaveTransformComponent] = filterFunction(tempWSST,tempFreqVec,varargin{1});
    else
        [waveTransformComponent,nonwaveTransformComponent] = filterFunction(tempWSST,tempFreqVec);
    end
    waveBeamVelocity(:,depthCtr) = iwsst(waveTransformComponent);
    nonwaveBeamVelocity(:,depthCtr) = iwsst(nonwaveTransformComponent);
end

end