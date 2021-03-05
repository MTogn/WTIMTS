%This function accepts a structure containing arrays of along-beam
%velocities in the as beamVelocityStructure.beamx for x = 1:4 or 1:5, with
%time along the first rank and depth along the second. It then truncates
%all beam records according to the length of the maximum useful bin -
%usually the highest bin before sidelobe interference makes measurements
%unreliable.
function beamVelocityStructure = truncateBeamVelocityDepthRange(beamVelocityStructure,maxUsefulBin)

beamVelocityStructure.beam1 = beamVelocityStructure.beam1(:,1:maxUsefulBin);
beamVelocityStructure.beam2 = beamVelocityStructure.beam2(:,1:maxUsefulBin);
beamVelocityStructure.beam3 = beamVelocityStructure.beam3(:,1:maxUsefulBin);
beamVelocityStructure.beam4 = beamVelocityStructure.beam4(:,1:maxUsefulBin);
if isfield(beamVelocityStructure,'beam5')
    beamVelocityStructure.beam5 = beamVelocityStructure.beam5(:,1:maxUsefulBin);
end

end