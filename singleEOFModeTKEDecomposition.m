%This function takes a whole array of ADCP estimates of TKE and decomposes
%it into a part due to turbulence and a part due to waves based on a
%previously-performed EOF decomposition (1st EOF and corresponding
%expansion coefficients must be supplied as arguments).
%Inputs:
%TKE - the array of TKE values as estimated by the ADCP. Must be QCd and
%zeroed to the surface.
%EOF1 - the first EOF for data in the array TKE. Must have the same size as
%each map in the array TKE
%EC1 - the EC corresponding to EOF1. Must be the same length as the time
%series in the array TKE.
%plotFlag - 1 to make plots, 0 to switch plots off.
%plotParams - structure containing time and depth vectors to make the
%plots, if needed.

function [TKETurb,TKEWave,figHand] = singleEOFModeTKEDecomposition(TKE,EOF1,EC1,plotFlag,plotParams)

%EOF can capture the time-varying pseudo-TKE introduced by wave action, but
%it cannot capture any bias introduced to the mean TKE estimate. In order
%to estimate what part of the mean TKE estimate is due to wave action, we
%recalculate the mean using only bursts when waves are not very strong
%i.e., when the EC of the 1st EOF is <0. The function that performs the
%separation accepts the whole-deployment TKE records and the expansion
%coefficient of ONLY the first EOF as arguments.
[meanTKEWave,meanTKETurb] = separateWaveTurbMeanTKE(TKE,EC1);

%The time-mean portion of the wave TKE is constant throughout the whole
%data record and the time-varying portion approximated by the EOF is simply
%added on top of this.
TKEWave = repmat(meanTKEWave,size(TKE,1),1);
TKEWave = TKEWave + repmat(EOF1',size(TKE,1),1).*repmat(EC1,1,size(TKE,2));
TKETurb = TKE - TKEWave;

%Plot the results now
if nargin > 4,
    if plotFlag == 1,
        
        if nargout > 2,
            figure(figHand), clf
        else
            figure
        end
        
        [depthArr,timeArr] = meshgrid(plotParams.depthVec,plotParams.timeVec)
        subplot(3,1,1), contourf(timeArr',depthArr',log10(TKE'))
        tempTKE = TKEWave; tempTKE(tempTKE < 0) = 0;
        subplot(3,1,2), contourf(timeArr',depthArr',log10(tempTKE'))
        tempTKE = TKETurb; tempTKE(tempTKE < 0) = 0;
        subplot(3,1,3), contourf(timeArr',depthArr',log10(tempTKE'))
        
    end
end

end