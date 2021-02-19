%Anything that needs to happen before the main loop over our shorter
%averaging periods goes here.

%For the original WADZ data set, preprocessing tasks are gathered into the
%function WADZPreprocessing, which includes:
% - Loading and formatting tilt data
tiltDataLong = WADZPreprocessing();
%Initialise some key parameters
paramStruc.beamAngle = 25*pi/180; paramStruc.anisoParam = 0.1684;
paramStruc.blankDist = 1.89; paramStruc.binVertSize = 0.6;
%The maximum number of depth bins should be known ahead of time in order to
%preallocate array sizes for depth-varying data. If the user does not know,
%it is better to err on the side of too many bins than too few, as the
%function call will fail if the preallocated arrays are too small.
maxDepthBins = 91;

%Preallocate the whole-record variables based on the number of bursts being
%analysed
burstStartIndex = 5;
burstEndIndex = 1470;
wholeRecordEnsNos = nan(burstEndIndex,2);
wholeRecordDatenums = nan(burstEndIndex,2);
wholeRecordTKE = nan(burstEndIndex,maxDepthBins);

%Burst loop
for burstCtr = burstStartIndex:burstEndIndex
    %Import a burst into the workspace
    [burstEnsembleNos,burstDatenums,burstBeamVelocities] = importWADZBurst(burstCtr,tiltDataLong);
    wholeRecordEnsNos(burstCtr,:) = [burstEnsembleNos(1) burstEnsembleNos(end)];
    wholeRecordDatenums(burstCtr,:) = [burstDatenums(1) burstDatenums(end)];

    %Carry out WSST
    %Filter
    %IWSST
    
    %Calculate TKE for unfiltered burst velocities
    wholeRecordTKE(burstCtr,:) = calcBurst4BeamTKE(burstBeamVelocities,paramStruc);
    
    %CalculateTKE for filtered burst velocities
    
    if rem(burstCtr,10) == 0,
        fprintf("Burst # is %d \r",burstCtr)
    end
    %End burst loop
end

%Carry out statistical filter
%We assume the wave pseudo-TKE to be a function of depth from surface,
%rather than height above bed. If the EOF is to capture this, we need to
%re-zero the dataset to reflect this, such that each point corresponds to
%the same number of bins below the surface rather than the same number of
%bins above the bed.
%For the WADZ data set, we use the pre-processed depth data in order to
%calculate the mean depth for a given burst.
[burstDepths,burstMaxBins] = demozoneDepthPreprocessing(paramStruc,wholeRecordEnsNos,burstStartIndex,burstEndIndex);
surfRelativeTKE = wholeRecordBed2Surf(wholeRecordTKE,burstMaxBins,burstStartIndex,burstEndIndex);

%Once the TKE array has been zeroed to the surface, we can perform the EOF
%analysis using the function EOF. numEOFs restricts the number of basis
%functions that are used in the analysis, as higher modes are unlikely to
%contain useful or meaningful information.

numEOFs = 10;
[TKEEigenvals,TKEEigenvalsNormd,TKEEOFs,TKEExpanCoeffs,TKETruncnErr] = EOFWrapper(surfRelativeTKE(burstStartIndex:burstEndIndex,:),numEOFs);

%EOF can capture the time-varying pseudo-TKE introduced by wave action, but
%it cannot capture any bias introduced to the mean TKE estimate. In order
%to estimate what part of the mean TKE estimate is due to wave action, we
%recalculate the mean using only bursts when waves are not very strong
%i.e., when the EC of the 1st EOF is <0. The function that performs the
%separation accepts the whole-deployment TKE records and the expansion
%coefficient of ONLY the first EOF as arguments.
[meanTKEWave,meanTKETurb] = separateWaveTurbMeanTKE(surfRelativeTKE(burstStartIndex:burstEndIndex,:),TKEExpanCoeffs(:,1));