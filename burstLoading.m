[burstEnsembleNos,burstDatenums,burstBeamVelocities,burstDepth] = importPartracBurst(paramStruc.dataLocation,burstCtr);
wholeRecordEnsNos(burstCtr,:) = [burstEnsembleNos(1) burstEnsembleNos(end)];
wholeRecordDatenums(burstCtr,:) = [burstDatenums(1) burstDatenums(end)];
[burstMeanDepths(burstCtr),burstMaxBins(burstCtr)] = partracDepthPreprocessing(paramStruc,burstDepth);
burstBeamVelocities = truncateBeamVelocityDepthRange(burstBeamVelocities,burstMaxBins(burstCtr));