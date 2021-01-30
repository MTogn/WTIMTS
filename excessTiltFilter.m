%This function takes a tilt record and velocity record corresponding to the
%same time range, removes any velocity data from times of excess tilt
%(defined vs. a threshold supplied as an argument to the function) and
%interpolates from good data to replace the measurements from times of
%excess tilt. Velocity, pitch and roll arrays should have their first
%dimension corresponding to time.
function tiltCorrectedBurstVelocity = excessTiltFilter(burstVelocityData,pitchData,rollData,tiltThreshold)

%Default tilt threshold is 5 degrees.
if nargin < 4,
    tiltThreshold = 5;
end

%For the RDI Workhorse, the pitch value that is supplied is not exactly
%what we need, so this needs to undergo a simple transformation.
pitchData = atan(tan(pitchData*pi/180).*cos(rollData*pi/180))*180/pi;
excessTiltEnsembles = find(pitchData > tiltThreshold | rollData > tiltThreshold);

%Firstly set initial or final velocities to zero if tilt is too great.
if (pitchData(1) > tiltThreshold | rollData(1) > tiltThreshold),
    burstVelocityData(1,:) = 0;
elseif (pitchData(end) > tiltThreshold | rollData(end) > tiltThreshold),
    burstVelocityData(end,:) = 0;
end

%Then step through the burst, interpolating along the first dimension if
%there is excess tilt.
ensCtr = 2;
while ensCtr < size(burstVelocityData(1))
    if (pitchData(ensCtr) > tiltThreshold | rollData(ensCtr) > tiltThreshold),
        interpStartEns = ensCtr;
        interpEndEns = interpStartEns;
        while (pitchData(interpEndEns + 1) > tiltThreshold | rollData(interpEndEns + 1) > tiltThreshold),
            interpEndEns = interpEndEns + 1;
        end
        burstVelocityData(interpStartEns:interpEndEns) = interp1([(interpStartEns - 1) (interpEndEns + 1)],...
                                                                    [burstVelocityData(interpStartEns - 1) burstVelocityData(interpEndEns + 1)],...
                                                                    interpStartEns:interpEndEns);
        ensCtr = interpEndEns + 1;
    else
        ensCtr = interpEndEns + 1;
    end
end

tiltCorrectedBurstVelocity = burstVelocityData;

%Function end
end
