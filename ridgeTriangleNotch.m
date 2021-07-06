function [filterStopComponent,filterPassComponent] = ridgeTriangleNotch(velocityWSST,wsstFreqVec,varargin);

%Function may be called without specifying the filter parameters (frequency
%width, max filter intensity etc.)
if nargin < 3,
    %Specify default filter parameters
    filterParameters.halfWidth = 5*(max(wsstFreqVec) - min(wsstFreqVec))/100;
    filterParameters.filterDepth = 1.2;
    filterParameters.maxSwellFreq = (1/3);
    filterParameters.wsstWaveThreshold = 0.01;
else
    filterParameters = varargin{1};
    filterParameters.halfWidth = filterParameters.halfWidthPercent*(max(wsstFreqVec) - min(wsstFreqVec))/100;
end

filterPassComponent = velocityWSST;
filterStopComponent = zeros(size(velocityWSST));
%Ridge finding
freqRidge = wsstridge(velocityWSST,wsstFreqVec);
for tCtr = 1:size(velocityWSST,2)
    filterMask = ones(length(wsstFreqVec),1);
%We will only apply the filter if the peak frequency corresponds to a
%period greater than 3s; any higher frequencies almost certainly do not
%correspond to a swell wave.
    if freqRidge(tCtr) < filterParameters.maxSwellFreq & max(abs(velocityWSST(:,tCtr))) > filterParameters.wsstWaveThreshold;
        filterMask(wsstFreqVec > (freqRidge(tCtr) - filterParameters.halfWidth) & wsstFreqVec <= freqRidge(tCtr)) = ...
            1 - filterParameters.filterDepth*(wsstFreqVec(wsstFreqVec > (freqRidge(tCtr) - filterParameters.halfWidth) & wsstFreqVec <= freqRidge(tCtr)) - (freqRidge(tCtr) - filterParameters.halfWidth))/filterParameters.halfWidth;
        filterMask(wsstFreqVec > freqRidge(tCtr) & wsstFreqVec < (freqRidge(tCtr) + filterParameters.halfWidth)) = ...
            1 - filterParameters.filterDepth*((freqRidge(tCtr) + filterParameters.halfWidth) - wsstFreqVec(wsstFreqVec > freqRidge(tCtr) & wsstFreqVec < (freqRidge(tCtr) + filterParameters.halfWidth)))/filterParameters.halfWidth;
    end
    filterMask(filterMask < 0) = 0;
    filterPassComponent(:,tCtr) = filterMask.*velocityWSST(:,tCtr);
    filterStopComponent(:,tCtr) = velocityWSST(:,tCtr) - filterPassComponent(:,tCtr);
end