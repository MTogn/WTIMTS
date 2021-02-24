function [waveTransformComponent,nonwaveTransformComponent] = ridgeTriangleNotch(velocityWSST,wsstFreqVec,filterParameters);

%Function may be called without specifying the filter parameters (frequency
%width, max filter intensity etc.)
if nargin < 3,
    %Specify default filter parameters
    filterParameters.halfWidth = 2.5*(max(wsstFreqVec) - min(wsstFreqVec))/length(wsstFreqVec);
    filterParameters.filterDepth = 0.9;
end

%Ridge finding
freqRidge = wsstridge(velocityWSST,wsstFreqVec);
for tCtr = 1:size(velocityWSST,2)
    filterMask = ones(length(wsstFreqVec),1);
    waveFreq = freqRidge(tCtr);
    filterMask(wsstFreqVec > (waveFreq - filterParameters.halfWidth) & testFreqVec <= waveFreq) = ...
            1 - filterDepth*(testFreqVec(testFreqVec > (waveFreq - filterHalfWidth) & testFreqVec <= waveFreq) - (waveFreq - filterHalfWidth))/filterHalfWidth;