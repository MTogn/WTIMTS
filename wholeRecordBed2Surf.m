%This function acts on a whole-record variable in order to change it from
%bed-relative (such that each entry along the depth dimension corresponds
%to the same distance from the seabed) to surface relative (such that each
%entry along that dimension corresponds to the same distance below the
%surface). This discards some near-bed data from records pertaining to
%times when the water column is at its deepest - this is unavoidable, as
%EOF analysis requires the same number of sample points in the space
%discretisation at all times.
function surfRelativeVariable = wholeRecordBed2Surf(bedRelativeVariable,binsByBurst,burstStartIndex,burstEndIndex)

if size(bedRelativeVariable,1) == burstEndIndex & size(bedRelativeVariable,2) ~= burstEndIndex,
    bedRelativeVariable = bedRelativeVariable';
end

surfRelativeVariable = nan(min(binsByBurst),size(bedRelativeVariable,2));
for burstCtr = burstStartIndex:burstEndIndex
    binsMoreThanMin = binsByBurst(burstCtr) - min(binsByBurst);
    surfRelativeVariable(:,burstCtr) = bedRelativeVariable((binsMoreThanMin + 1):binsByBurst(burstCtr),burstCtr);
end

end