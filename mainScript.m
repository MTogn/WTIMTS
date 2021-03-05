%Anything that needs to happen before the main loop over our shorter
%averaging periods goes here.

%For the original WADZ data set, preprocessing tasks are gathered into the
%function WADZPreprocessing, which includes:
% - Loading and formatting tilt and depth data
[tiltDataLong,demozoneRawDepth] = WADZPreprocessing();
%Initialise some key parameters
paramStruc.beamAngle = 25*pi/180; paramStruc.anisoParam = 0.1684;
paramStruc.blankDist = 1.89; paramStruc.binVertSize = 0.6;
paramStruc.sampFreq = 2;

%Preallocate the whole-record variables based on the number of bursts being
%analysed
burstStartIndex = 5;
burstEndIndex = 50;
wholeRecordEnsNos = nan(burstEndIndex,2);
wholeRecordDatenums = nan(burstEndIndex,2);
%It is helpful to know the maximum possible number of bins in advance to
%aid with preprocessing; if you don't know err on the side of caution
%(i.e., a larger number of bins) to avoid assigning data to arrays that are
%too small.
maxBinNo = 91;

%For the WADZ data set, we use the pre-processed depth data in order to
%calculate the mean depth for a given burst. This allows us to ensure that
%the time-consuming WSST filter is not applied to bins in the sidelobe
%interference range with unreliable data, or even to bins above the
%surface.
burstDepths = nan(burstEndIndex,1); burstMaxBins = nan(burstEndIndex,1);
wholeRecordTKE = nan(burstEndIndex,maxBinNo);
specFilteredWaveTKE = nan(burstEndIndex,maxBinNo);
specFilteredNonwaveTKE = nan(burstEndIndex,maxBinNo);

%Burst loop
for burstCtr = burstStartIndex:burstEndIndex
    %Import a burst into the workspace
    [burstEnsembleNos,burstDatenums,burstBeamVelocities] = importWADZBurst(burstCtr,tiltDataLong);
    wholeRecordEnsNos(burstCtr,:) = [burstEnsembleNos(1) burstEnsembleNos(end)];
    wholeRecordDatenums(burstCtr,:) = [burstDatenums(1) burstDatenums(end)];
    [burstDepths(burstCtr),burstMaxBins(burstCtr)] = demozoneDepthPreprocessing(paramStruc,demozoneRawDepth,wholeRecordEnsNos(burstCtr,:),burstCtr);
    burstBeamVelocities = truncateBeamVelocityDepthRange(burstBeamVelocities,burstMaxBins(burstCtr));

    %Carry out spectral filter
    [specFilteredWaveVelocities.beam1,specFilteredNonwaveVelocities.beam1] = spectralFilter(burstBeamVelocities.beam1,paramStruc.sampFreq,@ridgeTriangleNotch);
    [specFilteredWaveVelocities.beam2,specFilteredNonwaveVelocities.beam2] = spectralFilter(burstBeamVelocities.beam2,paramStruc.sampFreq,@ridgeTriangleNotch);
    [specFilteredWaveVelocities.beam3,specFilteredNonwaveVelocities.beam3] = spectralFilter(burstBeamVelocities.beam3,paramStruc.sampFreq,@ridgeTriangleNotch);
    [specFilteredWaveVelocities.beam4,specFilteredNonwaveVelocities.beam4] = spectralFilter(burstBeamVelocities.beam4,paramStruc.sampFreq,@ridgeTriangleNotch);
%     [specFilteredWaveVelocities.beam5,specFilteredNonwaveVelocities.beam5] = spectralFilter(burstBeamVelocities.beam5,paramStruc.sampFreq,@ridgeTriangleNotch);
    
    %Calculate TKE for unfiltered burst velocities
    wholeRecordTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(burstBeamVelocities,paramStruc);
    
    %Calculate TKE for filtered burst velocities
    specFilteredWaveTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(specFilteredWaveVelocities,paramStruc);
    specFilteredNonwaveTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(specFilteredNonwaveVelocities,paramStruc);
    
    if rem(burstCtr,10) == 0,
        fprintf("Burst # is %d \r",burstCtr)
    end
    %End burst loop
end
%% 

%Carry out statistical filter
%We assume the wave pseudo-TKE to be a function of depth from surface,
%rather than height above bed. If the EOF is to capture this, we need to
%re-zero the dataset to reflect this, such that each point corresponds to
%the same number of bins below the surface rather than the same number of
%bins above the bed.
surfRelativeTKE = wholeRecordBed2Surf(wholeRecordTKE,burstMaxBins,burstStartIndex,burstEndIndex);

%Once the TKE array has been zeroed to the surface, we can perform the EOF
%analysis using the function EOFWrapper. numEOFs restricts the number of
%basis functions that are used in the analysis, as higher modes are
%unlikely to contain useful or meaningful information. EOFWrapper includes
%a call to the main function EOF that actually performs the mode
%decomposition, plus some data hygiene and normalisation. Small EOF numbers
%can lead to the first EOFs/ECs being the opposite sign to what is expected
%- reason for this is unclear. Stick to larger numbers for now.
numEOFs = size(surfRelativeTKE,2);
[TKEEigenvals,TKEEigenvalsNormd,TKEEOFs,TKEExpanCoeffs,TKETruncnErr] = EOFWrapper(surfRelativeTKE(burstStartIndex:burstEndIndex,:),numEOFs);

%Some preprocessing of date/depth arrays so we have a unified way of
%plotting data that depends on time and/or depth.
depthBins = size(surfRelativeTKE,2);
plotParams.timeVec = wholeRecordDatenums(burstStartIndex:burstEndIndex,1);
meanDepth = nanmean(burstDepths); sidelobeDepth = meanDepth*(1 - cos(paramStruc.beamAngle));
plotParams.depthVec = -paramStruc.binVertSize*((depthBins - 1):-1:0) - sidelobeDepth;

[surfRelTKEWave,surfRelTKETurb] = singleEOFModeTKEDecomposition(surfRelativeTKE(burstStartIndex:burstEndIndex,:),TKEEOFs(:,1),TKEExpanCoeffs(:,1),1,plotParams);
%%
%Now repeat the statistical filter for the TKE dataset that has already
%undergone spectral filtering.