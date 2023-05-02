%This function calculates the error between two arrays, nominally of TKE
%values - one estimate from measured data and an estimate from wave theory.
%The function returns the mean absolute error per burst (or more generally
%the mean absolute error along the second dimension) and the relative error
%per burst as a percentage, referenced to the first argument.
function [relError,absError] = TKEArrayErrorCalc(TKEArr1,TKEArr2,burstIndexLims);

burstStartIndex = burstIndexLims(1);
burstEndIndex = burstIndexLims(2);

if size(TKEArr1) ~= size(TKEArr2),
    fprintf("Array size mismatch in error calculation \r")
    return
end

relError = nan(burstEndIndex,1);
absError = nan(burstEndIndex,1);
for burstCtr = burstStartIndex:burstEndIndex
   absError(burstCtr) = mean(abs(TKEArr1(burstCtr,:) - TKEArr2(burstCtr,:)),'omitnan');
   relError(burstCtr) = 100*absError(burstCtr)/mean(abs(TKEArr1(burstCtr,:)),'omitnan');
end

end