import FWCore.ParameterSet.Config as cms

L1TkElectrons = cms.EDProducer("L1TkElectronTrackProducer",
        #label = cms.string("ElecTrk"),
	label = cms.string("EG"),	# labels the collection of L1TkEmParticleProducer that is produced.
                                        # e.g. EG or IsoEG if all objects are kept, or
                                        # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                        # objects that pass a cut RelIso < RelIsoCut are written
                                        # into the new collection.
        L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),     # input EGamma collection
					# When the standard sequences are used :
                                                #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                #     or ("l1extraParticles","Isolated")
                                                #   - for the "old stage-2" algo (2x2 clustering), use 
                                                #     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
        ETmin = cms.double( -1.0 ),             # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
        L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
        # Quality cuts on Track and Track L1EG matching criteria                                
        TrackChi2           = cms.double(100.0), # minimum Chi2 to select tracks
        TrackMinPt          = cms.double(12.0), # minimum Pt to select tracks                                     
        TrackEGammaDeltaPhi = cms.vdouble(0.046,0.5,-0.25),  # functional Delta Phi cut parameters to match Track with L1EG objects
        TrackEGammaDeltaR   = cms.vdouble(0.066,0.27,-0.15), # functional Delta R cut parameters to match Track with L1EG objects
        TrackEGammaDeltaEta = cms.double(0.08), # Delta Eta cutoff to match Track with L1EG objects
                                                # are considered. 
	RelativeIsolation = cms.bool( True ),	# default = True. The isolation variable is relative if True,
						# else absolute.
        IsoCut = cms.double( -0.1 ), 		# Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
						# When RelativeIsolation = False, IsoCut is in GeV.
        # Determination of the isolation w.r.t. L1Tracks :
        PTMINTRA = cms.double( 2. ),	# in GeV
	DRmin = cms.double( 0.03),
	DRmax = cms.double( 0.2 ),
	DeltaZ = cms.double( 0.6 )    # in cm. Used for tracks to be used isolation calculation
)
