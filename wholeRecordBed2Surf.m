%This function acts on a whole-record variable in order to change it from
%bed-relative (such that each entry along the depth dimension corresponds
%to the same distance from the seabed) to surface relative (such that each
%entry along that dimension corresponds to the same distance below the
%surface). This discards some near-bed data from records pertaining to
%times when the water column is at its deepest - this is unavoidable, as
%EOF analysis requires the same number of sample points in the space
%discretisation at all times.
function surfRelativeVariable = wholeRecordBed2Surf(bedRelativeVariable,binsByRecord,recordStartNum,recordEndNum)

surfRelativeVariable = nan(size(bedRelativeVariable,1),min(binsByRecord));
for recordCtr = recordStartNum:recordEndNum
    binsMoreThanMin = binsByRecord(recordCtr) - min(binsByRecord);
    surfRelativeVariable(recordCtr,:) = bedRelativeVariable(recordCtr,(binsMoreThanMin + 1):binsByRecord(recordCtr));
end

end