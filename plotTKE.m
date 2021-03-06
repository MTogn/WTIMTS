%plotParams must contain at the fields plotParams.meanSurfRelDepthVec &
%plotParams.timeVec
%The TKE array is passed to the function with time in the first rank and
%depth in the second rank.

function [figHand,contourHand] = plotTKE(TKEArray,plotParams)

[depthArr,timeArr] = meshgrid(plotParams.meanSurfRelDepthVec,plotParams.timeVec);
tempTKE = TKEArray; tempTKE(tempTKE < 0) = 0;

if nargout == 0,
    figure
else
    figHand = figure;
end

if nargout == 1,
    contourf(timeArr',depthArr',log10(tempTKE'),'LineStyle','none')
else
    [dummy,contourHand] = contourf(timeArr',depthArr',log10(tempTKE'),'LineStyle','none');
end

datetick('x','mmm-dd','keepticks')
set(get(gca,'YLabel'),'String','Depth below mean surface level (m)')

end