function [waveTransformComponent,nonwaveTransformComponent] = ridgeTriangleNotch(velocityWSST,wsstFreqVec,varargin);

%Function may be called without specifying the filter parameters (frequency
%width, max filter intensity etc.)
if nargin < 3,
    %Specify default filter parameters
    filterParameters.halfWidth = 5*(max(wsstFreqVec) - min(wsstFreqVec))/100;
    filterParameters.filterDepth = 1.2;
    filterParameters.maxSwellFreq = 3;
    filterParameters.wsstWaveThreshold = 0.01;
else
    filterParameters = varargin{1};
    filterParameters.halfWidth = filterParameters.halfWidthPercent*(max(wsstFreqVec) - min(wsstFreqVec))/100;
end

nonwaveTransformComponent = velocityWSST;
waveTransformComponent = zeros(size(velocityWSST));
%Ridge finding
freqRidge = wsstridge(velocityWSST,wsstFreqVec);
for tCtr = 1:size(velocityWSST,2)
    filterMask = ones(length(wsstFreqVec),1);
    waveFreq = freqRidge(tCtr);
%We will only apply the filter if the peak frequency corresponds to a
%period greater than 3s; any higher frequencies almost certainly do not
%correspond to a swell wave.
    if waveFreq < 1/filterParameters.maxSwellFreq & max(abs(velocityWSST(:,tCtr))) > filterParameters.wsstWaveThreshold;
        filterMask(wsstFreqVec > (waveFreq - filterParameters.halfWidth) & wsstFreqVec <= waveFreq) = ...
            1 - filterParameters.filterDepth*(wsstFreqVec(wsstFreqVec > (waveFreq - filterParameters.halfWidth) & wsstFreqVec <= waveFreq) - (waveFreq - filterParameters.halfWidth))/filterParameters.halfWidth;
        filterMask(wsstFreqVec > waveFreq & wsstFreqVec < (waveFreq + filterParameters.halfWidth)) = ...
            1 - filterParameters.filterDepth*((waveFreq + filterParameters.halfWidth) - wsstFreqVec(wsstFreqVec > waveFreq & wsstFreqVec < (waveFreq + filterParameters.halfWidth)))/filterParameters.halfWidth;
    end
    filterMask(filterMask < 0) = 0;
    nonwaveTransformComponent(:,tCtr) = filterMask.*velocityWSST(:,tCtr);
    waveTransformComponent(:,tCtr) = velocityWSST(:,tCtr) - nonwaveTransformComponent(:,tCtr);
end