%This function will calculate the wavenumber for waves observed from a
%point observation of waves in significant currents. The user must supply a
%number of parameters:
%observedFreq - the observed frequency (the function assumes you do not
%   have the intrinsic frequency) in rad/s
%waterDepth - in m
%currSpeed - the speed of the current in m/s
%waveCurrAngle - the angle between the wave direction and the current
%   direction in rad

%There are also two optional parameters:
%convergeMargin - the convergence criterion, expressed as a marginal
%   change from one iteration to the next, in percent
%maxIters - maximum number of iterations before declaring an unconverged
%   result

%The function sets a default number of iterations as 100 and a convergence
%criterion of 0.1% if the user does not specify otherwise
function kappa = iterCalcWvno_WavesCurrents(observedFreq,...
                                            waterDepth,...
                                            currSpeed,...
                                            waveCurrAngle,...
                                            convergeMargin,...
                                            maxIters)

gravAccel = 9.81;
%Set default values of optional parameters
if nargin < 6,
    maxIters = 100;
end
if nargin < 5,
    convergeMargin = 0.1;
end
%Check that all necessary inputs are supplied
if nargin < 4,
    error('Insufficient data supplied - wavenumber calculation requires the user to specify the observed frequency, the water depth, the mean current speed and the angle between waves and current.')
end

%First specify an initial guess for the wavenumber, based on the deep water
%approximation for the wavenumber.
kappaOld = (gravAccel + 2*observedFreq*currSpeed*cos(waveCurrAngle)) - sqrt(gravAccel*(gravAccel + 4*observedFreq*currSpeed*cos(waveCurrAngle)));
kappaOld = kappaOld/(2*(currSpeed^2)*(cos(waveCurrAngle)^2));
%This should not be complex, so set it to a typical order-of-magnitude
%guess if it is.
if imag(kappaOld) ~= 0, kappaOld = 0.1; end

%This is the main iteration loop
iterCtr = 1;
while iterCtr < maxIters
    kappaNew = observedFreq^2 + (currSpeed*kappaOld*cos(waveCurrAngle))^2;
    kappaNew = kappaNew/(2*observedFreq*currSpeed*cos(waveCurrAngle) + gravAccel*tanh(kappaOld*waterDepth));
    marginPercent = abs((kappaNew - kappaOld)/kappaOld);
    if marginPercent <= convergeMargin
        break
    end
    kappaOld = kappaNew;
    iterCtr = iterCtr + 1;
end
kappa = kappaNew;

if iterCtr == maxIters
    'Wavenumber unconverged'
end

end