%The opposite of wholeRecordBed2Surf, this takes a surface-zeroed array
%from an EOF analysis and reorders it so that each bin is positioned along
%the depth rank relative to its distance from the bed. This will mean that
%there are some near-bed bins for which no data is available; these are
%represented as nans.
function bedRelativeVariable = wholeRecordSurf2Bed(surfRelativeVariable,binsByRecord,recordStartNum,recordEndNum)

bedRelativeVariable = nan(size(surfRelativeVariable,1),max(binsByRecord(recordStartNum:recordEndNum)));
for recordCtr = recordStartNum:recordEndNum
    binsMoreThanMin = binsByRecord(recordCtr) - min(binsByRecord);
    bedRelativeVariable(recordCtr,(1 + binsMoreThanMin):binsByRecord(recordCtr)) = surfRelativeVariable(recordCtr,:);
end

end