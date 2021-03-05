figure
subplot(2,2,1), contourf(specFilteredWaveVelocities.beam1(:,1:60)','LineStyle','none')
subplot(2,2,2), contourf(specFilteredWaveVelocities.beam2(:,1:60)','LineStyle','none')
subplot(2,2,3), contourf(specFilteredWaveVelocities.beam3(:,1:60)','LineStyle','none')
subplot(2,2,4), contourf(specFilteredWaveVelocities.beam4(:,1:60)','LineStyle','none')

figure
subplot(2,2,1), contourf(specFilteredNonwaveVelocities.beam1(:,1:60)','LineStyle','none')
subplot(2,2,2), contourf(specFilteredNonwaveVelocities.beam2(:,1:60)','LineStyle','none')
subplot(2,2,3), contourf(specFilteredNonwaveVelocities.beam3(:,1:60)','LineStyle','none')
subplot(2,2,4), contourf(specFilteredNonwaveVelocities.beam4(:,1:60)','LineStyle','none')