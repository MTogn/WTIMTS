function [filterStoppedBeamVelocity,filterPassedBeamVelocity] = spectralFilter(beamVelocity,sampFreq,filterFunction,varargin)

filterStoppedBeamVelocity = nan(size(beamVelocity)); filterPassedBeamVelocity = filterStoppedBeamVelocity;
for depthCtr = 1:size(beamVelocity,2)
%Transform each velocity record into its pseudo-periodogram
    [tempWSST,tempFreqVec] = wsst(beamVelocity(:,depthCtr),sampFreq);
%Call filter function here
%If the optional argument in is supplied (i.e., nargin > 3), then it should
%be a structure containing key parameters to be passed to the filter.
    if nargin > 3,
        [filterStoppedComponent,filterPassedComponent] = filterFunction(tempWSST,tempFreqVec,varargin{1});
    else
        [filterStoppedComponent,filterPassedComponent] = filterFunction(tempWSST,tempFreqVec);
    end
%After the filter has split the WSST representation into its two
%components, use the inverse WSST to transform back into the time domain.
    filterStoppedBeamVelocity(:,depthCtr) = iwsst(filterStoppedComponent);
    filterPassedBeamVelocity(:,depthCtr) = iwsst(filterPassedComponent);
end

end