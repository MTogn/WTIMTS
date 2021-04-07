beamVelocity = burstBeamVelocities.beam1(:,50); 
sampFreq = 2;
[tempWSST,tempFreqVec] = wsst(beamVelocity,sampFreq);
tempWSSTp1 = zeros(size(tempWSST)); tempWSSTp2 = tempWSSTp1;
tempWSSTp1(1:2:end,:) = tempWSST(1:2:end,:); tempWSSTp2(2:2:end,:) = tempWSST(2:2:end,:);

beamVelyP1 = iwsst(tempWSSTp1); beamVelyP2 = iwsst(tempWSSTp2);
reconstructedBeamVely = beamVelyP1 + beamVelyP2 + mean(beamVelocity);

%Check that the spectra look like spectra
dT = 0.5; sampFreq = 1/dT;
recordLength = length(beamVelocity);
recordLengthSec = recordLength*dT;
dF = 1/recordLengthSec;

freqVec1Sided = dF*(0:1:recordLength);
nyquistFreq = 1/(2*dT); nyquistFreqIndex = find(freqVec1Sided == nyquistFreq);

rawPSD = periodogram(beamVelocity,[],freqVec1Sided,1/dT);
p1PSD = periodogram(beamVelyP1,[],freqVec1Sided,1/dT);
p2PSD = periodogram(beamVelyP2,[],freqVec1Sided,1/dT);
figure, plot(freqVec1Sided(1:nyquistFreqIndex),rawPSD(1:nyquistFreqIndex),...
                freqVec1Sided(1:nyquistFreqIndex),p1PSD(1:nyquistFreqIndex),...
                freqVec1Sided(1:nyquistFreqIndex),p2PSD(1:nyquistFreqIndex))
title('Non-physical spectral decomposition'), legend('Total','Odds','Evens')
set(gca,'XScale','log','YScale','log')

[tempWSSTp3,tempWSSTp4] = ridgeTriangleNotch(tempWSST,tempFreqVec);
beamVelyP3 = iwsst(tempWSSTp3); beamVelyP4 = iwsst(tempWSSTp4);
reconstructedBeamVelyv2 = beamVelyP3 + beamVelyP4 + mean(beamVelocity);

rawPSD = periodogram(beamVelocity,[],freqVec1Sided,1/dT);
wavePSD = periodogram(beamVelyP3,[],freqVec1Sided,1/dT);
nonwavePSD = periodogram(beamVelyP4,[],freqVec1Sided,1/dT);
figure, plot(freqVec1Sided(1:nyquistFreqIndex),rawPSD(1:nyquistFreqIndex),...
                freqVec1Sided(1:nyquistFreqIndex),wavePSD(1:nyquistFreqIndex),...
                freqVec1Sided(1:nyquistFreqIndex),nonwavePSD(1:nyquistFreqIndex))
title('Physical spectral decomposition'), legend('Total','Wave','Non-wave')
set(gca,'XScale','log','YScale','log')

figure, subplot(2,1,1), plot(beamVelocity,'k--','LineWidth',2), hold on, plot(beamVelyP1), plot(beamVelyP2), plot(reconstructedBeamVely), set(gca,'XLim',[0 120])
subplot(2,1,2), plot(beamVelocity,'k--','LineWidth',2), hold on, plot(beamVelyP3), plot(beamVelyP4), plot(reconstructedBeamVelyv2), set(gca,'XLim',[0 120])

%%
%Common scaling for contour plots + corresponding colorbars
set(wholeRecordTKECont,'LevelList',uniformLevelList); caxis(get(wholeRecordTKECont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(EOFOnlyTurbCont,'LevelList',uniformLevelList); caxis(get(EOFOnlyTurbCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(EOFOnlyWaveCont,'LevelList',uniformLevelList); caxis(get(EOFOnlyWaveCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(WSSTOnlyTurbCont,'LevelList',uniformLevelList); caxis(get(WSSTOnlyTurbCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(WSSTOnlyWaveCont,'LevelList',uniformLevelList); caxis(get(WSSTOnlyWaveCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(bothFilterTurbCont,'LevelList',uniformLevelList); caxis(get(bothFilterTurbCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(bothFilterWaveCont,'LevelList',uniformLevelList); caxis(get(bothFilterWaveCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;
set(wavePseudoTurbCont,'LevelList',uniformLevelList); caxis(get(wavePseudoTurbCont,'Parent'),[uniformLevelList(1) uniformLevelList(end)]), colorbar;

%%
%Testing wsst magnitudes for a magnitude test on the filter
for ctr = 1:size(burstBeamVelocities.beam1,2)
    [tempWSSTArrays{ctr},tempWSSTFreqVec] = wsst(burstBeamVelocities.beam1(:,ctr),paramStruc.sampFreq);
    maxWSST(ctr) = max(max(abs(tempWSSTArrays{ctr})));
end