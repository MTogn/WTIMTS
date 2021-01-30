burstStartIndex = 5;
burstEndIndex = 1470;

%Burst loop
for burstCtr = burstStartIndex:burstEndIndex
    
    %Import a burst into the workspace
    [burstDatenums,burstBeamVelocities] = importWADZBurst(burstCtr)

    %Calculate TKE for burst
    
    %Carry out WSST
    %Filter
    %IWSST

    %Calc burst turbulence
    
%End burst loop
end

%Carry out statistical filter