function condErr = errorECConditional(rawErr,expanCoeff,burstIndexRange,ECCondition)

burstStartIndex = burstIndexRange(1); burstEndIndex = burstIndexRange(2);

%Ensure both the non-conditional error and the expansion coefficients are
%row vectors of length burstEndIndex; if either are not, they are front-
%and end-padded with nans.
if ~isrow(rawErr), rawErr = rawErr'; end
if ~isrow(expanCoeff), expanCoeff = expanCoeff'; end

if length(rawErr) ~= burstEndIndex,
    rawErr = [nan(1,burstStartIndex -1) rawErr nan(1,burstEndIndex - length(rawErr))];
end
if length(expanCoeff) ~= burstEndIndex,
    expanCoeff = [nan(1,burstStartIndex -1) expanCoeff nan(1,burstEndIndex - length(rawErr))];
end

condErr = nan(1,burstEndIndex);
condErr(expanCoeff > ECCondition) = rawErr(expanCoeff > ECCondition);
end