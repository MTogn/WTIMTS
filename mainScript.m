%Anything that needs to happen before the main loop over our shorter
%averaging periods goes here.
tiltDataLong = WADZPreprocessing();

burstStartIndex = 10;
burstEndIndex = 20;

%Burst loop
for burstCtr = burstStartIndex:burstEndIndex
    tic
    %Import a burst into the workspace
    [burstDatenums,burstBeamVelocities] = importWADZBurst(burstCtr,tiltDataLong);

    %Calculate TKE for burst
    
    %Carry out WSST
    %Filter
    %IWSST

    %Calc burst turbulence
    
    toc
%End burst loop
end

%Carry out statistical filter