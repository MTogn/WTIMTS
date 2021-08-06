%NaN filter for ADCP data; v2 includes a second output that counts how many
%NaNs have been replaced for each filtered column
function [FilteredArr,NaNCount] = NaNFilterv2(UnfilteredArr)

FilteredArr = UnfilteredArr;
NaNCount = zeros(1,size(FilteredArr,2));

%Some of our datasets have set all the bad data to -32768, let's turn all these to NaNs.
FilteredArr(FilteredArr == -32678) = NaN;

%Firstly set any initial or final NaNs to 0
for ctr2 = 1:size(FilteredArr,2)
    if isnan(FilteredArr(1,ctr2)) == 1,
        FilteredArr(1,ctr2) = 0;
        NaNCount(ctr2) = NaNCount(ctr2) + 1;
    end
    if isnan(FilteredArr(size(FilteredArr,1),ctr2)) == 1,
        FilteredArr(size(FilteredArr,1),ctr2) = 0;
        NaNCount(ctr2) = NaNCount(ctr2) + 1;
    end
end

%This subroutine replaces NaNs with interpolated values (N.B., this does
%not try to determine whether a bin has returned good data or not, that is
%dealt with outside this function; here we interpolate data whether or not it 
%ends up being meaningful).
for ctr1 = 1:size(FilteredArr,1)
    for ctr2 = 1:size(FilteredArr,2)
        if isnan(FilteredArr(ctr1,ctr2)) == 1,
            StartPt = ctr1;
            EndPt = StartPt;
            NaNCount(ctr2) = NaNCount(ctr2) + 1;
            while isnan(FilteredArr(EndPt + 1,ctr2)) == 1,
                EndPt = EndPt + 1;
                NaNCount(ctr2) = NaNCount(ctr2) + 1;
            end
            if EndPt == StartPt,
                FilteredArr(ctr1,ctr2) = 0.5*(FilteredArr(ctr1 + 1,ctr2) + FilteredArr(ctr1 - 1,ctr2));
            else
                for ctr3 = StartPt:EndPt
                    FilteredArr(ctr3,ctr2) = (ctr3 - StartPt + 1)*(FilteredArr(StartPt - 1,ctr2) + FilteredArr(EndPt + 1,ctr2))/(EndPt - StartPt + 2);
                end
            end                    
        end
    end
end