[burstEnsembleNos,burstDatenums,burstBeamVelocities] = importWADZBurst(burstCtr,tiltDataLong);
wholeRecordEnsNos(:,burstCtr) = [burstEnsembleNos(1); burstEnsembleNos(end)];
wholeRecordDatenums(:,burstCtr) = [burstDatenums(1); burstDatenums(end)];
[burstDepths(burstCtr),burstMaxBins(burstCtr)] = demozoneDepthPreprocessing(paramStruc,demozoneRawDepth,wholeRecordEnsNos(:,burstCtr),burstCtr);
burstBeamVelocities = truncateBeamVelocityDepthRange(burstBeamVelocities,burstMaxBins(burstCtr));