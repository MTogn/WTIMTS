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

%Next determine the wave amplitude - take this as half the significant wave
%height. burstWaveHeights records the SWH in cm, so an additional factor of
%100 is needed.
load([paramStruc.dataLocation 'burstWaveHeights.mat']);
burstWaveAmp = burstWaveHeightsDespiked'/200;

%Next determine the surface current speed. Take this to be the magnitude of
%the velocity at the bin closest to the surface for which there is data
%available (data in the top ~15% of the column cannot be measured due to
%sidelobe interference).
load([paramStruc.dataLocation 'burstDepthSurfBins.mat']);
%burstMeanVelsDirs must include at least:
%   - a 2D array velDirn of size (maxBinNos,numBursts) containing the mean
%flow direction for each height in the water column for each burst
%   - a 2D array velMag of the same size containing the mean flow magnitude
%for the same locations
load([paramStruc.dataLocation 'burstMeanVelsDirs.mat']);
%Note that some bursts with bad data have a recorded max bin of 0 - this
%doesn't play well with subscript assignment, so we relabel these as being
%at the highest bin possible.
burstMaxBins_noZeros = burstMaxBins;
burstMaxBins_noZeros(burstMaxBins_noZeros == 0) = size(velMag,1);
surfBinIndices = sub2ind(size(velDirn),burstMaxBins_noZeros,1:burstEndIndex);
surfCurrentSpeed = velMag(surfBinIndices)';
burstMeanDepths = burstMeanDepths';

%Determine the directions of the waves and the surface current, and the
%difference between them. This is an important input for calculating the
%wave wavenumber and intrinsic frequency.
load([paramStruc.dataLocation 'burstWaveHeadings.mat']), burstWaveHeadings = burstWaveHeadings';
surfCurrentDirn = velDirn(surfBinIndices);
currWaveDirnDiff = wrapToPi((surfCurrentDirn - burstWaveHeadings)*pi/180);

%The loop below is the main calculation; for each burst, we calculate the
%wavenumber from the wave buoy observed frequency, depth and surface
%current speed and direction. From this we get the intrinsic frequency, the
%amplitude of the wave pseudo-TKE and, finally, the profile of wave
%pseudo-TKE itself.
burstWvno = nan(burstEndIndex,1);
burstIntrinsWaveFreq = nan(burstEndIndex,1);
kWaveAmplitude = nan(burstEndIndex,1);
%For the array of wave pseudo-TKE profiles, we need to know the number of
%bins in the profile; this is also used inside the loop to calculate the
%burst-wise depths.
minBinWithData = min(burstMaxBins_noZeros);
waveBuoyPseudoTKE = nan(burstEndIndex,minBinWithData);
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
    waveBuoyPseudoTKE(burstCtr,:) = kWaveAmplitude(burstCtr)*cosh(2*burstWvno(burstCtr)*(depthFromSurf + burstMeanDepths(burstCtr))); 
end

end