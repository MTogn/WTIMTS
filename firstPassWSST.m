dT = 0.5; sampFreq = 2;
testTimeSeries = b1Fluc(end,:);
recordLength = length(testTimeSeries);
testSampleTimes = dT*(1:recordLength);
[testWSST,testFreqVec] = wsst(testTimeSeries,sampFreq);

waveFreqRidge = wsstridge(testWSST,testFreqVec);

figure
pcolor(testSampleTimes,testFreqVec,abs(testWSST));
shading interp, caxis([0 0.05])
xlabel('Seconds'), ylabel('Frequency (Hz)');
hold on, plot(testSampleTimes,waveFreqRidge,'r.')

testWSST(testFreqVec > 0.09 & testFreqVec < 0.125,:) = 0.25*testWSST(testFreqVec > 0.09 & testFreqVec < 0.125,:);
figure
pcolor(testSampleTimes,testFreqVec,abs(testWSST));
shading interp, caxis([0 0.05])
xlabel('Seconds'), ylabel('Frequency (Hz)');

testTimeSeriesRecon = iwsst(testWSST);
testPSD = periodogram(testTimeSeries,[],freqVec1Sided,1/dT);
testPSDRecon = periodogram(testTimeSeriesRecon,[],freqVec1Sided,1/dT);

figure
subplot(2,1,1)
plot(dT:dT:recordLengthSec,testTimeSeries), hold on, plot(dT:dT:recordLengthSec,testTimeSeriesRecon)
set(get(gca,'XLabel'),'String','Time (s)'), set(get(gca,'YLabel'),'String','Beam 1 fluctuation velocity (ms^{-1})')
subplot(2,1,2)
plot(freqVec1Sided(1:nyquistFreqIndex),testPSD(1:nyquistFreqIndex))
hold on, plot(freqVec1Sided(1:nyquistFreqIndex),testPSDRecon(1:nyquistFreqIndex))
set(gca,'XScale','log','YScale','log')
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')

b1FlucDewaved = b1Fluc;
for binCtr = 1:numBins
    testTimeSeries = b1Fluc(binCtr,:);
    [testWSST,testFreqVec] = wsst(testTimeSeries,sampFreq);
    testWSST(testFreqVec > 0.09 & testFreqVec < 0.125,:) = 0.25*testWSST(testFreqVec > 0.09 & testFreqVec < 0.125,:);
    dewavedTimeSeries = iwsst(testWSST);
    b1FlucDewaved(binCtr,:) = dewavedTimeSeries;
end

freqVec1Sided = dF*(0:1:recordLength);
nyquistFreq = 1/(2*dT); nyquistFreqIndex = find(freqVec1Sided == nyquistFreq);

% testTimeSeries = b1Fluc';
% testPSD = periodogram(testTimeSeries,[],freqVec1Sided,1/dT);
% figure
% subplot(2,2,1), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(1,1):quartileIndices(1,2)),2))
% set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1e-1])
% set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
% subplot(2,2,2), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(2,1):quartileIndices(2,2)),2))
% set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1e-1])
% set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
% subplot(2,2,3), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(3,1):quartileIndices(3,2)),2))
% set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1e-1])
% set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
% subplot(2,2,4), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(4,1):quartileIndices(4,2)),2))
% set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1e-1])
% set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')

testTimeSeries = b1FlucDewaved';
testPSD = periodogram(testTimeSeries,[],freqVec1Sided,1/dT);
figure
subplot(2,2,1), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(1,1):quartileIndices(1,2)),2))
title('0-25% depth')
set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1])
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
subplot(2,2,2), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(2,1):quartileIndices(2,2)),2))
title('25-50% depth')
set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1])
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
subplot(2,2,3), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(3,1):quartileIndices(3,2)),2))
title('50-75% depth')
set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1])
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')
subplot(2,2,4), plot(freqVec1Sided(1:nyquistFreqIndex),mean(testPSD(1:nyquistFreqIndex,quartileIndices(4,1):quartileIndices(4,2)),2))
title('75-100% depth')
set(gca,'XScale','log','YScale','log','XLim',[1e-3 1],'YLim',[1e-4 1])
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')