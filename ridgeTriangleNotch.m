function [waveTransformComponent,nonwaveTransformComponent] = ridgeTriangleNotch(velocityWSST,wsstFreqVec,varargin);

%Function may be called without specifying the filter parameters (frequency
%width, max filter intensity etc.)
if nargin < 3,
    %Specify default filter parameters
    filterParameters.halfWidth = 2.5*(max(wsstFreqVec) - min(wsstFreqVec))/length(wsstFreqVec);
    filterParameters.filterDepth = 0.9;
else
    filterParameters = varargin{1};
end

nonwaveTransformComponent = velocityWSST;
waveTransformComponent = zeros(size(velocityWSST));
%Ridge finding
freqRidge = wsstridge(velocityWSST,wsstFreqVec);
for tCtr = 1:size(velocityWSST,2)
    filterMask = ones(length(wsstFreqVec),1);
    waveFreq = freqRidge(tCtr);
    filterMask(wsstFreqVec > (waveFreq - filterParameters.halfWidth) & wsstFreqVec <= waveFreq) = ...
            1 - filterParameters.filterDepth*(wsstFreqVec(wsstFreqVec > (waveFreq - filterParameters.halfWidth) & wsstFreqVec <= waveFreq) - (waveFreq - filterParameters.halfWidth))/filterParameters.halfWidth;
    filterMask(wsstFreqVec > waveFreq & wsstFreqVec < (waveFreq + filterParameters.halfWidth)) = ...
        1 - filterParameters.filterDepth*((waveFreq + filterParameters.halfWidth) - wsstFreqVec(wsstFreqVec > waveFreq & wsstFreqVec < (waveFreq + filterParameters.halfWidth)))/filterParameters.halfWidth;
    nonwaveTransformComponent(:,tCtr) = filterMask.*velocityWSST(:,tCtr);
    waveTransformComponent(:,tCtr) = velocityWSST(:,tCtr) - nonwaveTransformComponent(:,tCtr);
end