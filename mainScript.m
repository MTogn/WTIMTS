%Anything that needs to happen before the main loop over our shorter
%averaging periods goes here.

%makePlots tells the code whether to plot results as they are calculated or
%not.
makePlots = true;
%Not all cases need mean flow data, or may have mean flow data already
%calculated; this section is therefore controlled with a flag.
calcMeanFlowFlag = true;
%calcErrorFlag tells the code whether to calculate the error between
%analytic wave pseudo-TKE estimates and the estimates from the filters.
calcErrorFlag = true;

%Preallocate the whole-record variables based on the number of bursts being
%analysed
burstStartIndex = 50;
burstEndIndex = 8958;
wholeRecordEnsNos = nan(burstEndIndex,2);
wholeRecordDatenums = nan(burstEndIndex,2);
burstMeanDepths = nan(burstEndIndex,1);
burstMaxBins = nan(burstEndIndex,1);
%It is helpful to know the maximum possible number of bins in advance to
%aid with preprocessing; if you don't know err on the side of caution
%(i.e., a larger number of bins) to avoid assigning data to arrays that are
%too small.
maxBinNo = 90;

%If your data requires any preprocessing, alter the script
%dataPreprocessing to contain all the relevant commands and initialisation
%of key parameters - in particular, the structure paramStruc should be
%defined inside dataPreprocessing. paramStruc should define the directory
%from which data is drawn as paramStruc.dataLocation
dataPreprocessing

%If data is to be saved automatically, switch the flag to true and indicate
%the save location
saveDataFlag = true;
paramStruc.saveDirectory = 'C:\Users\michael\Documents\ADCP\NWDZ_north\Results\';
paramStruc.saveFilename = 'completeWorkspace0p9Deep15Wide.mat';

%Preallocate whole-record variables whose size depends on both the number
%of bursts and the number of bins.
wholeRecordADCPTKE = nan(burstEndIndex,maxBinNo);
specFilterStoppedTKE = nan(burstEndIndex,maxBinNo);
specFilterPassedTKE = nan(burstEndIndex,maxBinNo);

%If there are any parameters to be set for the filter, they can be
%set immediately in advance of the main loop; this may be better moved to a
%general "initialise values" script, but otoh the fewer scripts that have
%to be modified by a user the better.
filterParameters.halfWidthPercent = 7.5;
filterParameters.filterDepth = 0.9;
filterParameters.maxSwellFreq = (1/3);
filterParameters.wsstWaveThreshold = 0.02;

%%

switch calcMeanFlowFlag
    case 1
        %velLat and velLong are the lateral and longitudinal velocities
        %relative to the ADCP i.e., the device-relative y- and
        %x-velocities.
        velLat = nan(maxBinNo,burstEndIndex);
        velLong = nan(maxBinNo,burstEndIndex);
        velMag = nan(maxBinNo,burstEndIndex);
        velDirn = nan(maxBinNo,burstEndIndex);
        for burstCtr = burstStartIndex:burstEndIndex
            burstLoading

            burstMeanFlow = calcMeanFlowProperties(burstBeamVelocities,burstHeadingRaw,paramStruc);
            velLong(1:burstMaxBins(burstCtr),burstCtr) = burstMeanFlow.uMean(1:burstMaxBins(burstCtr));
            velLat(1:burstMaxBins(burstCtr),burstCtr) = burstMeanFlow.vMean(1:burstMaxBins(burstCtr));
            velMag(1:burstMaxBins(burstCtr),burstCtr) = sqrt(velLong(1:burstMaxBins(burstCtr),burstCtr).^2 + velLat(1:burstMaxBins(burstCtr),burstCtr).^2);
            velDirn(1:burstMaxBins(burstCtr),burstCtr) = burstMeanFlow.flowHeadingMean;

            if rem(burstCtr,10) == 0,
                fprintf("Burst # is %d \r",burstCtr)
            end

        end

        save([paramStruc.dataLocation 'burstMeanVelsDirs'],'velLong','velLat','velMag','velDirn');
end
%%
%Burst loop
for burstCtr = burstStartIndex:burstEndIndex
    %Call the script burstLoading to load a burst into the workspace in the
    %correct format. This script is likely to need modification in order to
    %suit the format in which your ADCP data is stored. Check the user
    %guide, section 4.1.2 to see the format the the burst data must be in
    %before beginning the spectral filter.
    burstLoading

    %Carry out spectral filter
    [specFilterStoppedVelocities.beam1,specFilterPassedVelocities.beam1] = spectralFilter(burstBeamVelocities.beam1,paramStruc.sampFreq,@ridgeTriangleNotch,filterParameters);
    [specFilterStoppedVelocities.beam2,specFilterPassedVelocities.beam2] = spectralFilter(burstBeamVelocities.beam2,paramStruc.sampFreq,@ridgeTriangleNotch,filterParameters);
    [specFilterStoppedVelocities.beam3,specFilterPassedVelocities.beam3] = spectralFilter(burstBeamVelocities.beam3,paramStruc.sampFreq,@ridgeTriangleNotch,filterParameters);
    [specFilterStoppedVelocities.beam4,specFilterPassedVelocities.beam4] = spectralFilter(burstBeamVelocities.beam4,paramStruc.sampFreq,@ridgeTriangleNotch,filterParameters);
    if isfield(burstBeamVelocities,'beam5'),
        [specFilteredWaveVelocities.beam5,specFilteredNonwaveVelocities.beam5] = spectralFilter(burstBeamVelocities.beam5,paramStruc.sampFreq,@ridgeTriangleNotch);
    end
    
    %Calculate TKE for unfiltered burst velocities
    wholeRecordADCPTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(burstBeamVelocities,paramStruc);
    
    %Calculate TKE for filtered burst velocities
    specFilterStoppedTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(specFilterStoppedVelocities,paramStruc);
    specFilterPassedTKE(burstCtr,1:burstMaxBins(burstCtr)) = calcBurst4BeamTKE(specFilterPassedVelocities,paramStruc);
    
    if rem(burstCtr,10) == 0,
        fprintf("Burst # is %d \r",burstCtr)
    end
    %End burst loop
end
%% 
%Carry out statistical filter

%The first data set to which we apply the statistical filter is the TKE
%calculated 'naively' from the ADCP data, without any spectral filtering.
%We assume the wave pseudo-TKE to be a function of depth from surface,
%rather than height above bed. If the EOF is to capture this, we need to
%re-zero the dataset to reflect this, such that each point corresponds to
%the same number of bins below the surface rather than the same number of
%bins above the bed.
surfRelADCPTKE = wholeRecordBed2Surf(wholeRecordADCPTKE,burstMaxBins,burstStartIndex,burstEndIndex);

%Once the TKE array has been zeroed to the surface, we can perform the EOF
%analysis using the function EOFWrapper. numEOFs restricts the number of
%basis functions that are used in the analysis, as higher modes are
%unlikely to contain useful or meaningful information. EOFWrapper includes
%a call to the main function EOF that actually performs the mode
%decomposition, plus some data hygiene and normalisation. Small EOF numbers
%can lead to the first EOFs/ECs being the opposite sign to what is expected
%- reason for this is unclear. Stick to larger numbers for now.
numEOFs = size(surfRelADCPTKE,2);
[TKEEigenvals,TKEEigenvalsNormd,TKEEOFs,TKEExpanCoeffs,TKETruncnErr] = EOFWrapper(surfRelADCPTKE(burstStartIndex:burstEndIndex,:),numEOFs);
%EOF decomposition of the unfiltered ADCP estimate of TKE
[surfRelTKETurb,surfRelTKEWave] = singleEOFModeTKEDecomposition(surfRelADCPTKE(burstStartIndex:burstEndIndex,:),TKEEOFs(:,1),TKEExpanCoeffs(:,1));

%Now repeat the statistical filter for the TKE dataset that has already
%undergone spectral filtering to remove (part of) the wave pseudo-TKE.
surfRelFilterPassedTKE = wholeRecordBed2Surf(specFilterPassedTKE,burstMaxBins,burstStartIndex,burstEndIndex);
numEOFs = size(surfRelFilterPassedTKE,2);
[TKEEigenvals,TKEEigenvalsNormd,TKEEOFs,TKEExpanCoeffs,TKETruncnErr] = EOFWrapper(surfRelFilterPassedTKE(burstStartIndex:burstEndIndex,:),numEOFs);
%EOF decomposition of the spectrally-filtered ADCP estimate of TKE
[surfRelFilterPassedTKETurb,surfRelFilterPassedTKEWave] = singleEOFModeTKEDecomposition(surfRelFilterPassedTKE(burstStartIndex:burstEndIndex,:),TKEEOFs(:,1),TKEExpanCoeffs(:,1));

%This part of the routine takes the outputs of the EOF decomposition and
%returns it to a bed-zeroed format. Since some data was lost or excerpted
%to make the original inputs to the EOF routine, we have to pad the output
%with nans at the bed AND at the start, if there are any initial bursts
%that are excluded.
bedRelTKETurb = wholeRecordSurf2Bed([nan(burstStartIndex - 1,size(surfRelTKETurb,2)); surfRelTKETurb],burstMaxBins,burstStartIndex,burstEndIndex);
bedRelTKEWave = wholeRecordSurf2Bed([nan(burstStartIndex - 1,size(surfRelTKEWave,2)); surfRelTKEWave],burstMaxBins,burstStartIndex,burstEndIndex);
%We also bed-zero the EOF decomposition of the double-filtered TKE arrays.
bedRelDblFilteredTKETurb = wholeRecordSurf2Bed([nan(burstStartIndex - 1,size(surfRelFilterPassedTKETurb,2)); surfRelFilterPassedTKETurb],burstMaxBins,burstStartIndex,burstEndIndex);
bedRelDblFilteredTKEWave = wholeRecordSurf2Bed([nan(burstStartIndex - 1,size(surfRelFilterPassedTKEWave,2)); surfRelFilterPassedTKEWave],burstMaxBins,burstStartIndex,burstEndIndex);
%To get the double-filtered estimate of wave pseudo-TKE, we must add
%together the spectral-filter-stopped TKE component and the EOF estimate of
%the wave portion of the spectral-filter-passed TKE component.
bedRelDblFilteredTKEWave = bedRelDblFilteredTKEWave + specFilterStoppedTKE(:,1:max(burstMaxBins));

%%
%Some preprocessing of date/depth arrays so we have a unified way of
%plotting data that depends on time and/or depth.
minNumBinsWithData = size(surfRelADCPTKE,2);
maxNumBinsWithData = max(burstMaxBins(burstStartIndex:burstEndIndex));
%The datenum record must be a column vector
if size(wholeRecordDatenums,1) == 1, wholeRecordDatenums = wholeRecordDatenums'; end
plotParams.timeVec = wholeRecordDatenums(burstStartIndex:burstEndIndex,1);
%recordMeanDepth is the mean depth across all bursts; sidelobeDepth is the depth
%range near the surface of the water column from which we cannot obtain
%useful data due to sidelobe interference.
recordMeanDepth = mean(burstMeanDepths,'omitnan'); sidelobeDepth = recordMeanDepth*(1 - cos(paramStruc.beamAngle));
burstSurfBins = 1 + floor((burstMeanDepths - paramStruc.blankDist)/paramStruc.binVertSize);
%surfRelDepthVec is the depth vector relative to the surface for any burst
%individually; meanSurfRelDepthVec is the depth vector relative to the mean
%surface depth
plotParams.surfRelDepthVec = -paramStruc.binVertSize*((minNumBinsWithData - 1):-1:0) - sidelobeDepth;
plotParams.meanSurfRelDepthVec = paramStruc.binVertSize*(0:1:(maxNumBinsWithData - 1)) + (paramStruc.blankDist - recordMeanDepth);

%We can now plot all the TKE arrays in separate figures:
% - 'Naive' ADCP estimate of TKE
% - EOF-only estimate of true turbulent TKE
% - EOF-only estimate of wave pseudo-TKE
% - WSST-only estimate of true turbulent TKE
% - WSST-only estimate of wave pseudo-TKE
% - Double-filtered estimate of true turbulent TKE
% - Double-filtered estimate of wave pseudo-TKE
%including both wave and turbulent contributions
switch makePlots
    case 1
        [wholeRecordTKEFig,wholeRecordTKEAx,wholeRecordTKECont] = plotTKE(wholeRecordADCPTKE(burstStartIndex:burstEndIndex,1:maxNumBinsWithData),plotParams);
        title(wholeRecordTKEAx,'Unfiltered ADCP estimate of TKE')

        [EOFOnlyTurbFig,EOFOnlyTurbAx,EOFOnlyTurbCont] = plotTKE(bedRelTKETurb(burstStartIndex:burstEndIndex,:),plotParams);
        title(EOFOnlyTurbAx,'EOF-only estimate of turbulent k (log_{10}, J\cdotkg^{-1})')

        [EOFOnlyWaveFig,EOFOnlyWaveAx,EOFOnlyWaveCont] = plotTKE(bedRelTKEWave(burstStartIndex:burstEndIndex,:),plotParams);
        title(EOFOnlyWaveAx,'EOF-only estimate of wave pseudo-k (log_{10}, J\cdotkg^{-1})')

        [WSSTOnlyTurbFig,WSSTOnlyTurbAx,WSSTOnlyTurbCont] = plotTKE(specFilterPassedTKE(burstStartIndex:burstEndIndex,1:maxNumBinsWithData),plotParams);
        title(WSSTOnlyTurbAx,'WSST-only estimate of turbulent k (log_{10}, J\cdotkg^{-1})')

        [WSSTOnlyWaveFig,WSSTOnlyWaveAx,WSSTOnlyWaveCont] = plotTKE(specFilterStoppedTKE(burstStartIndex:burstEndIndex,1:maxNumBinsWithData),plotParams);
        title(WSSTOnlyWaveAx,'WSST-only estimate of wave pseudo-k (log_{10}, J\cdotkg^{-1})')

        [bothFilterTurbFig,bothFilterTurbAx,bothFilterTurbCont] = plotTKE(bedRelDblFilteredTKETurb(burstStartIndex:burstEndIndex,:),plotParams);
        title(bothFilterTurbAx,'Estimate of turbulent k with both filters (log_{10}, J\cdotkg^{-1})')

        [bothFilterWaveFig,bothFilterWaveAx,bothFilterWaveCont] = plotTKE(bedRelDblFilteredTKEWave(burstStartIndex:burstEndIndex,:),plotParams);
        title(bothFilterWaveAx,'Estimate of wave pseudo-k with both filters (log_{10}, J\cdotkg^{-1})')
end

%%
switch saveDataFlag
    case 1
        saveVarsList = mkSaveVarsList;
        save([paramStruc.saveDirectory paramStruc.saveFilename],saveVarsList{:});
end
%%
%This section of the code calculates or imports the expected wave
%pseudo-TKE if it exists, and calculates the error vs. the filtered TKE.
switch calcErrorFlag
    case 1
    
        anycWavePseudoTKE = calcApparentWaveTKE([burstStartIndex burstEndIndex],paramStruc);
        bedRelPseudoTKE = wholeRecordSurf2Bed(anycWavePseudoTKE,burstMaxBins,burstStartIndex,burstEndIndex);
        if makePlots == true,
            [wavePseudoTurbFig,wavePseudoTurbCont] = plotTKE(bedRelPseudoTKE(burstStartIndex:burstEndIndex,:),plotParams);
            title(get(wavePseudoTurbFig,'Children'),'Estimate of wave pseudo-k from AWT (log_{10}, J\cdotkg^{-1})')
        end

%First the error vs. the statistical-only filtered wave condition. As some
%initia bursts had to be excluded from the EOF analysis, the array of EOF
%estimates of wave pseudo-TKE (surfRelFilteredTKEWave) has to be padded to
%match the dimension of the calculated pseudo-TKE.
        [EOFOnlyRelError,EOFOnlyAbsError] = TKEArrayErrorCalc([nan(burstStartIndex - 1,size(surfRelTKEWave,2)); surfRelTKEWave],anycWavePseudoTKE,[burstStartIndex burstEndIndex]);

%To get the error calculation for the spectral-only filtered estimate of
%wave pseudo-TKE, we first need to zero to the surface.
        wsstOnlySurfZeroWaveTKE = wholeRecordBed2Surf(specFilterStoppedTKE,burstMaxBins,burstStartIndex,burstEndIndex);
        [WSSTOnlyRelError,WSSTOnlyAbsError] = TKEArrayErrorCalc(wsstOnlySurfZeroWaveTKE,anycWavePseudoTKE,[burstStartIndex burstEndIndex]);

%The double-filtered wave pseudo-TKE is obtained by taking the WSST
%filtered wave pseudo-TKE, and then adding the fraction of the remainder
%attributed to wave action by the subsequent EOF filter. As in the EOF-only
%case, this requires some initial padding to account for the bursts
%excluded in the EOF analysis.
        bothFilterSurfZeroWaveTKE = wsstOnlySurfZeroWaveTKE + [nan(burstStartIndex - 1,size(surfRelFilterPassedTKEWave,2)); surfRelFilterPassedTKEWave];
        [bothFilterRelError,bothFilterAbsError] = TKEArrayErrorCalc(bothFilterSurfZeroWaveTKE,anycWavePseudoTKE,[burstStartIndex burstEndIndex]);

end