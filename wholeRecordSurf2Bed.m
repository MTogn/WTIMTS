%The opposite of wholeRecordBed2Surf, this takes a surface-zeroed array
%from an EOF analysis and reorders it so that each bin is positioned along
%the depth rank relative to its distance from the bed. This will mean that
%there are some near-bed bins for which no data is available; these are
%represented as nans.
function bedRelativeVariable = wholeRecordSurf2Bed(surfRelativeVariable,binsByBurst,recordStartBurst,recordEndBurst)

bedRelativeVariable = nan(size(surfRelativeVariable,1),max(binsByBurst(recordStartBurst:recordEndBurst)));
for recordCtr = recordStartBurst:recordEndBurst
    binsMoreThanMin = binsByBurst(recordCtr) - min(binsByBurst);
    bedRelativeVariable(recordCtr,(1 + binsMoreThanMin):binsByBurst(recordCtr)) = surfRelativeVariable(recordCtr,:);
end

end