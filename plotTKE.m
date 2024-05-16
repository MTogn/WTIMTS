%plotParams must contain at least the fields plotParams.meanSurfRelDepthVec
%and plotParams.timeVec
%The TKE array is passed to the function with time in the first rank and
%depth in the second rank.

function [figHand,axHand,contourHand] = plotTKE(TKEArray,plotParams)

[depthArr,timeArr] = meshgrid(plotParams.meanSurfRelDepthVec,plotParams.timeVec);
tempTKE = TKEArray; tempTKE(tempTKE < 0) = 0;

switch nargout
    case 0
        figure;
        contourf(timeArr',depthArr',log10(tempTKE'),'LineStyle','none');
    case 1
        figHand = figure;
        contourf(timeArr',depthArr',log10(tempTKE'),'LineStyle','none');
    otherwise
        figHand = figure;
        [dummy,contourHand] = contourf(timeArr',depthArr',log10(tempTKE),'LineStyle','none');
        axHand = gca;
end
colorbar

datetick('x','mmm-dd','keepticks','keeplimits')
set(get(gca,'YLabel'),'String','Depth below mean surface level (m)')

end