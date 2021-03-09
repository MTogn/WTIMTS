%plotParams must contain at the fields plotParams.meanSurfRelDepthVec &
%plotParams.timeVec
%The TKE array is passed to the function with time in the first rank and
%depth in the second rank.

function figHand = plotTKE(TKEArray,plotParams)

[depthArr,timeArr] = meshgrid(plotParams.meanSurfRelDepthVec,plotParams.timeVec);
tempTKE = TKEArray; tempTKE(tempTKE < 0) = 0;
figure, contourf(timeArr',depthArr',log10(tempTKE'),'LineStyle','none')
datetick('x','mmm-dd','keepticks')
set(get(gca,'YLabel'),'String','Depth below mean surface level (m)')

end