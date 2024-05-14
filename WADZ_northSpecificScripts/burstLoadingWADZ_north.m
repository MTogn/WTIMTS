[burstEnsembleNos,burstDatenums,burstBeamVelocities] = importWADZBurst_north(burstCtr,tiltDataLong);
wholeRecordEnsNos(burstCtr,:) = [burstEnsembleNos(1) burstEnsembleNos(end)];
wholeRecordDatenums(burstCtr,:) = [burstDatenums(1) burstDatenums(end)];
%paramStruc is defined in the main script, demozoneRawDepth is loaded
%during the call to dataPreprocessingWADZ_north (aliased to
%dataPreprocessing) in mainScript
[burstMeanDepths(burstCtr),burstMaxBins(burstCtr)] = demozoneDepthPreprocessing_north(paramStruc,demozoneRawDepth,wholeRecordEnsNos(burstCtr,:),burstCtr);
burstBeamVelocities = truncateBeamVelocityDepthRange(burstBeamVelocities,burstMaxBins(burstCtr));