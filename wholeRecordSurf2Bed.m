%The opposite of wholeRecordBed2Surf, this takes a surface-zeroed array
%from an EOF analysis and reorders it so that each bin is positioned along
%the depth rank relative to its distance from the bed. This will mean that
%there are some near-bed bins for which no data is available; these are
%represented as nans.
function bedRelativeVariable = wholeRecordSurf2Bed(surfRelativeVariable,binsByBurst,burstStartIndex,burstEndIndex)

if size(surfRelativeVariable,1) == burstEndIndex & size(surfRelativeVariable,2) ~= burstEndIndex,
    surfRelativeVariable = surfRelativeVariable';
end

bedRelativeVariable = nan(max(binsByBurst(burstStartIndex:burstEndIndex)),burstEndIndex);
for recordCtr = burstStartIndex:burstEndIndex
    binsMoreThanMin = binsByBurst(recordCtr) - min(binsByBurst);
    bedRelativeVariable((1 + binsMoreThanMin):binsByBurst(recordCtr),recordCtr) = surfRelativeVariable(:,recordCtr);
end

end