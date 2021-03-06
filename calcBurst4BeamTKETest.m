%Inputs
%burstVelocities - structure containing (at least) time records of
%along-beam velocity for four off-verical beams of an ADCP. Each beam
%record must be the same duration and same number of bins. The required
%entries in the structure are burstVelocities.beam1, burstVelocities.beam2,
%burstVelocities.beam3 and burstVelocities.beam4.
%calcParams - structure containing (at least) the angle of the off-vertical
%beams (calcParams.beamAngle) and the parameter defining the proportion of
%TKE assumed to be due to vertical fluctuations (calcParams.anisoParam). A
%default value is supplied in the function for the anisotropy parameter
%below but not for the beam angle.
function [beamVars,burst4BeamTKE] = calcBurst4BeamTKETest(burstVelocities,calcParams)

if ~exist('calcParams.anisoParam'),
    calcParams.anisoParam = 0.1684;
end

%For the beam velocity data, time should be the first dimension and depth
%the second.
if size(burstVelocities.beam1,1) < size(burstVelocities.beam1,2)
    burstVelocities.beam1 = burstVelocities.beam1';
end
if size(burstVelocities.beam2,1) < size(burstVelocities.beam2,2)
    burstVelocities.beam2 = burstVelocities.beam2';
end
if size(burstVelocities.beam3,1) < size(burstVelocities.beam3,2)
    burstVelocities.beam3 = burstVelocities.beam3';
end
if size(burstVelocities.beam4,1) < size(burstVelocities.beam4,2)
    burstVelocities.beam4 = burstVelocities.beam4';
end

b1Var = var(burstVelocities.beam1,1,1);
b2Var = var(burstVelocities.beam2,1,1);
b3Var = var(burstVelocities.beam3,1,1);
b4Var = var(burstVelocities.beam4,1,1);
beamVars.b1Var = b1Var; beamVars.b2Var = b2Var; beamVars.b3Var = b3Var; beamVars.b4Var = b4Var;

burst4BeamTKE = b1Var + b2Var + b3Var + b4Var;
burst4BeamTKE = burst4BeamTKE/(4*(sin(calcParams.beamAngle)^2)*(1 + (2*(cot(calcParams.beamAngle)^2) - 1)*calcParams.anisoParam));

end