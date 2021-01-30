demozonePath = 'C:\Users\michael\Documents\ADCP\DEMOZONE\';
dataPath = 'burstData\burstData';
load([demozonePath 'burstEnsembleNumbers'])
load([demozonePath 'burstDepthSurfBins'])
load([demozonePath 'burstExcessTiltRecords'])
burstCtr = 5;
load([demozonePath dataPath int2str(burstCtr)])

burstLength = size(burstData,1)/5;
numBins = burstMaxBin(burstCtr);

%Preallocating the size of all these variables speeds everything up by a
%fairly chunky margin
b1Raw = zeros(burstLength,numBins);
b2Raw = zeros(burstLength,numBins);
b3Raw = zeros(burstLength,numBins);
b4Raw = zeros(burstLength,numBins);
bVertRaw = zeros(burstLength,numBins);
headVec = zeros(burstLength);

b1Raw = burstData(1:5:(5*burstLength - 4),9:(numBins + 8));
b2Raw = burstData(2:5:(5*burstLength - 3),9:(numBins + 8));
b3Raw = burstData(3:5:(5*burstLength - 2),9:(numBins + 8));
b4Raw = burstData(4:5:(5*burstLength - 1),9:(numBins + 8));
bVertRaw = burstData(5:5:5*burstLength,9:(numBins + 8));

b1Raw = NaNFilterv2(b1Raw);, b1Raw = b1Raw'/1000;
b2Raw = NaNFilterv2(b2Raw);, b2Raw = b2Raw'/1000;
b3Raw = NaNFilterv2(b3Raw);, b3Raw = b3Raw'/1000;
b4Raw = NaNFilterv2(b4Raw);, b4Raw = b4Raw'/1000;
bVertRaw= NaNFilterv2(bVertRaw);, bVertRaw = bVertRaw'/1000;

excessTiltSamps = listSamplesExcessTilt(burstCtr).excessTiltSamps;
if isempty(excessTiltSamps) == 0,
    for badSampCtr = numSamplesExcessTilt(burstCtr):-1:1
        badSamp = excessTiltSamps(badSampCtr);
        b1Raw = b1Raw(:,[1:(badSamp - 1) (badSamp + 1):size(b1Raw,2)]);
        b2Raw = b2Raw(:,[1:(badSamp - 1) (badSamp + 1):size(b2Raw,2)]);
        b3Raw = b3Raw(:,[1:(badSamp - 1) (badSamp + 1):size(b3Raw,2)]);
        b4Raw = b4Raw(:,[1:(badSamp - 1) (badSamp + 1):size(b4Raw,2)]);
        bVertRaw = bVertRaw(:,[1:(badSamp - 1) (badSamp + 1):size(bVertRaw,2)]);
    end
end

b1Mean = mean(b1Raw,2); b2Mean = mean(b2Raw,2);
    b3Mean = mean(b3Raw,2); b4Mean = mean(b4Raw,2);
    bVertMean = mean(bVertRaw,2);
b1Fluc = b1Raw - repmat(b1Mean,1,size(b1Raw,2)); b2Fluc = b2Raw - repmat(b2Mean,1,size(b2Raw,2));
    b3Fluc = b3Raw - repmat(b3Mean,1,size(b3Raw,2)); b4Fluc = b4Raw - repmat(b4Mean,1,size(b4Raw,2));
    bVertFluc = bVertRaw - repmat(bVertMean,1,size(bVertRaw,2));

dT = 0.5; sampFreq = 1/dT;
recordLength = length(b1Fluc);
recordLengthSec = recordLength*dT;
dF = 1/recordLengthSec;

freqVec1Sided = dF*(0:1:recordLength);
nyquistFreq = 1/(2*dT); nyquistFreqIndex = find(freqVec1Sided == nyquistFreq);
testTimeSeries = b1Fluc';
testPSD = periodogram(testTimeSeries,[],freqVec1Sided,1/dT);

quartileIndices(1,1) = 1; quartileIndices(1,2) = floor(numBins/4);
quartileIndices(2,1) = quartileIndices(1,2) + 1; quartileIndices(2,2) = floor(numBins/2);
quartileIndices(3,1) = quartileIndices(2,2) + 1; quartileIndices(3,2) = floor(3*numBins/4);
quartileIndices(4,1) = quartileIndices(3,2) + 1; quartileIndices(4,2) = numBins;

figure
subplot(2,1,1)
plot(dT:dT:recordLengthSec,testTimeSeries(:,end))
set(get(gca,'XLabel'),'String','Time (s)'), set(get(gca,'YLabel'),'String','Beam 1 fluctuation velocity (ms^{-1})')
subplot(2,1,2)
plot(freqVec1Sided(1:nyquistFreqIndex),testPSD(1:nyquistFreqIndex,end));
set(gca,'XScale','log','YScale','log')
set(get(gca,'XLabel'),'String','Frequency (Hz)'), set(get(gca,'YLabel'),'String','Velocity PSD (m^{2}s^{-2}Hz^{-1})')

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