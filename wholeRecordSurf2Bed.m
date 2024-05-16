%The opposite of wholeRecordBed2Surf, this takes a surface-zeroed array
%from an EOF analysis and reorders it so that each bin is positioned along
%the depth rank relative to its distance from the bed. This will mean that
%there are some near-bed bins for which no data is available; these are
%represented as nans.
function bedRelativeVariable = wholeRecordSurf2Bed(surfRelativeVariable,binsByBurst,recordStartBurst,recordEndBurst)

bedRelativeVariable = nan(max(binsByBurst(recordStartBurst:recordEndBurst)),recordEndBurst);
for recordCtr = recordStartBurst:recordEndBurst
    binsMoreThanMin = binsByBurst(recordCtr) - min(binsByBurst);
    bedRelativeVariable((1 + binsMoreThanMin):binsByBurst(recordCtr),recordCtr) = surfRelativeVariable(:,recordCtr);
end

end