%This function determines the apparent TKE due to wave action, given
%certain  wave parameters - specifically, wave frequency, wave amplitude,
%the surface current and the directions of the wave and currents. Its
%intended output is an array of wave pseudo-TKE of the same size as the
%filtered TKE arrays based on ADCP measurements.
function waveBuoyPseudoTKE = calcApparentWaveTKE(burstIndexRange,paramStruc)

%burstIndexRange is a two-element array that specifies the first and last
%burst number for which we must calculated the expected pseudo-TKED
burstStartIndex = burstIndexRange(1); burstEndIndex = burstIndexRange(2);

%First determine the frequency in radians per second. Note this is the
%observed frequency under combined wave-current conditions, not the
%intrinsic frequency.
load([paramStruc.dataLocation 'burstWaveModalPeriod.mat']);
burstObservedWaveFreq = 2*pi./burstWaveModalPeriodDespiked;
%Several variables relating to current and wave properties in this code are
%found to not be reliably in the correct row order. The first rank should
%be along the water column, and the second rank along the bursts. Thus,
%variables such as burstObservedWaveFreq (with one value per burst) should
%be row vectors. This is repeated for similar variables throughout.
if ~isrow(burstObservedWaveFreq), burstObservedWaveFreq = burstObservedWaveFreq'; end

%Wave amplitude is half of significant wave height. Waverider buoys record
%SWH in cm, so an additional factor of 100 is needed.
load([paramStruc.dataLocation 'burstWaveHeights.mat']);
burstWaveAmp = burstWaveHeightsDespiked/200;
if ~isrow(burstWaveAmp), burstWaveAmp = burstWaveAmp'; end

%Surface current speed is approximated as the magnitude of velocity at the
%bin closest to the surface for which there is data available
load([paramStruc.dataLocation 'burstDepthSurfBins.mat']);
%burstMeanVelsDirs must include at least:
%   - a 2D array velDirn of size (maxBinNos,numBursts) containing the mean
%flow direction for each height in the water column for each burst
%   - a 2D array velMag of the same size containing the mean flow magnitude
%for the same locations
load([paramStruc.dataLocation 'burstMeanVelsDirs.mat']);
%Note that some bursts with bad data may have a recorded max bin of 0 - 
%relabel these as being at the highest bin possible.
burstMaxBins_noZeros = burstMaxBins;
burstMaxBins_noZeros(burstMaxBins_noZeros == 0) = size(velMag,1);
if ~isrow(burstMaxBins_noZeros), burstMaxBins_noZeros = burstMaxBins_noZeros'; end
surfBinIndices = sub2ind(size(velDirn),burstMaxBins_noZeros,1:burstEndIndex);

%Some bursts may have current or wave data as nans due to bad data; this
%ensures that we get only the data corresponding to bursts with all surface
%and wave data.
surfCurrentSpeed = nan(size(surfBinIndices));
surfCurrentSpeed(~isnan(surfBinIndices)) = velMag(surfBinIndices(~isnan(surfBinIndices)));

if ~isrow(surfCurrentSpeed), surfCurrentSpeed = surfCurrentSpeed'; end
if ~isrow(burstMeanDepths), burstMeanDepths = burstMeanDepths'; end

%Determine the directions of the waves and the surface current, and the
%difference between them. This is an important input for calculating the
%wave wavenumber and intrinsic frequency.
load([paramStruc.dataLocation 'burstWaveHeadings.mat']), burstWaveHeadings = burstWaveHeadings';
surfCurrentDirn = nan(size(surfBinIndices));
surfCurrentDirn(~isnan(surfBinIndices)) = velDirn(surfBinIndices(~isnan(surfBinIndices)));
currWaveDirnDiff = wrapToPi((surfCurrentDirn - burstWaveHeadings)*pi/180);

%The loop below is the main calculation; for each burst, we calculate the
%wavenumber from the wave buoy observed frequency, depth and surface
%current speed and direction. From this we get the intrinsic frequency, the
%amplitude of the wave pseudo-TKE and, finally, the profile of wave
%pseudo-TKE itself.
burstWvno = nan(1,burstEndIndex);
burstIntrinsWaveFreq = nan(1,burstEndIndex);
kWaveAmplitude = nan(1,burstEndIndex);
%For the array of wave pseudo-TKE profiles, we need to know the number of
%bins in the profile; this is also used inside the loop to calculate the
%burst-wise depths.
minBinWithData = min(burstMaxBins_noZeros);
waveBuoyPseudoTKE = nan(minBinWithData,burstEndIndex);
for burstCtr = burstStartIndex:burstEndIndex
    burstWvno(burstCtr) = iterCalcWvno_WavesCurrents(burstObservedWaveFreq(burstCtr),...
                                        burstMeanDepths(burstCtr),...
                                        surfCurrentSpeed(burstCtr),...
                                        currWaveDirnDiff(burstCtr));
    burstIntrinsWaveFreq(burstCtr) = burstObservedWaveFreq(burstCtr) - burstWvno(burstCtr)*surfCurrentSpeed(burstCtr)*cos(currWaveDirnDiff(burstCtr));
    kWaveAmplitude(burstCtr) = (0.5*burstIntrinsWaveFreq(burstCtr)*burstWaveAmp(burstCtr)/sinh(burstWvno(burstCtr)*burstMeanDepths(burstCtr)))^2;
    if isnan(kWaveAmplitude(burstCtr)), kWaveAmplitude(burstCtr) = 0; end
%Depth range is slightly different for each burst, as the sidelobe depth
%depends on the mean water depth during that burst. This has to be
%simplified when plotting the whole-record data.
    sidelobeDepth = burstMeanDepths(burstCtr)*(1 - cos(paramStruc.beamAngle));
    depthFromSurf = -paramStruc.binVertSize*((minBinWithData - 1):-1:0) - sidelobeDepth;
    waveBuoyPseudoTKE(:,burstCtr) = kWaveAmplitude(burstCtr)*cosh(2*burstWvno(burstCtr)*(depthFromSurf + burstMeanDepths(burstCtr))); 
end

end