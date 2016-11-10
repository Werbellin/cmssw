import FWCore.ParameterSet.Config as cms

##########################################################################
## TREE CONTENT
#########################################################################
    
ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    abseta = cms.string("abs(eta)"),
    pt  = cms.string("pt"),
    mass  = cms.string("mass"),
    )   

SCProbeVariablesToStore = cms.PSet(
    probe_sc_eta    = cms.string("eta"),
    probe_sc_abseta = cms.string("abs(eta)"),
    probe_sc_pt     = cms.string("pt"),
    probe_sc_et     = cms.string("et"),
    probe_sc_e      = cms.string("energy"),

    probe_sc_tkIso  = cms.InputTag("recoEcalCandidateHelper:scTkIso"),
    )

EleProbeVariablesToStore = cms.PSet(
    probe_Ele_eta    = cms.string("eta"),
    probe_Ele_abseta = cms.string("abs(eta)"),
    probe_Ele_phi    = cms.string("phi"),
    probe_Ele_pt     = cms.string("pt"),
    probe_Ele_et     = cms.string("et"),
    probe_Ele_e      = cms.string("energy"),
    probe_Ele_q      = cms.string("charge"),
    
    ## super cluster quantities
    probe_sc_energy     = cms.string("superCluster.energy"),
    probe_sc_rawE       = cms.string("superCluster.rawEnergy"),
    probe_sc_preshowerE = cms.string("superCluster.preshowerEnergy"),
    probe_sc_et         = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_eta        = cms.string("-log(tan(superCluster.position.theta/2))"),
    probe_sc_abseta     = cms.string("abs(-log(tan(superCluster.position.theta/2)))"),
    probe_sc_phi        = cms.string("superCluster.phi"),    
    probe_sc_phiW       = cms.string("superCluster.phiWidth"),
    probe_sc_etaW       = cms.string("superCluster.etaWidth"),

    #id based
#    probe_Ele_dEtaSeeOut       = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    probe_Ele_dEtaIn        = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_Ele_dPhiIn        = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_Ele_dEtaSeed      = cms.string("deltaEtaSuperClusterTrackAtVtx+log(tan(superCluster.position.theta/2))-log(tan(superCluster.seed.position.theta/2))"),
    probe_Ele_sieie         = cms.string("sigmaIetaIeta"),
    probe_Ele_e1x5          = cms.string("e1x5"),
    probe_Ele_e2x5          = cms.string("e2x5Max"),
    probe_Ele_e5x5          = cms.string("e5x5"),
    probe_Ele_r9            = cms.string("r9"),
    probe_Ele_r9_5x5        = cms.string("full5x5_r9"),
    probe_Ele_sieie_5x5     = cms.string("full5x5_sigmaIetaIeta"),
    probe_Ele_hoe           = cms.string("hadronicOverEm"),
    probe_Ele_ooemoop       = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),

    probe_Ele_chisq         = cms.InputTag("eleVarHelper:chisq"),
    probe_Ele_mHits         = cms.InputTag("eleVarHelper:missinghits"),
    probe_Ele_dz            = cms.InputTag("eleVarHelper:dz"),
    probe_Ele_dxy           = cms.InputTag("eleVarHelper:dxy"),
     
    #isolation
    probe_Ele_ecalIso       = cms.string("ecalPFClusterIso"),
    probe_Ele_hcalIso       = cms.string("hcalPFClusterIso"),

    # tracker
    probe_Ele_trkPt         = cms.string("gsfTrack().ptMode"), 
    
    probe_ele_Spring15_NonTrig = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    probe_ele_Spring16 = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16V1Values"),
    probe_ele_Spring16_HZZ_WP = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-V1-wpLoose"),
    probe_ele_Spring16_90_WP = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-V1-wp90"),
    probe_ele_Spring16_80_WP = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-V1-wp80"),

    #probe_ele_Spring15_Trig = cms.InputTag("X"),
    #probe_ele_Spring15_Trig_Classic = cms.InputTag("X"),


    probe_ele_deltaPhiSeedClusterTrackAtCalo = cms.string("deltaPhiSeedClusterTrackAtCalo"),
    probe_ele_deltaEtaSeedClusterTrackAtCalo = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    probe_ele_deltaPhiSuperClusterTrackAtVtx = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_ele_deltaEtaSuperClusterTrackAtVtx = cms.string("deltaEtaSuperClusterTrackAtVtx"),


    probe_ele_eSuperClusterOverP = cms.string("eSuperClusterOverP"),
    probe_ele_eEleClusterOverPout = cms.string("eEleClusterOverPout"),
    probe_ele_eSeedClusterOverP = cms.string("eSeedClusterOverP"),
    probe_ele_eSeedClusterOverPout = cms.string("eSeedClusterOverPout"),

    probe_ele_trackMomentumAtVtx = cms.string("trackMomentumAtVtx.R"),
    probe_ele_trackMomentumAtVtx_pT = cms.string("trackMomentumAtVtx.Rho"),
    probe_ele_trackMomentumOut = cms.string("trackMomentumOut.R"),
    probe_ele_trackMomentumOut_pT = cms.string("trackMomentumOut.Rho"),


    probe_ele_full_5x5_sigmaIEtaIEta = cms.string("full5x5_sigmaIetaIeta"),
    probe_ele_full_5x5_sigmaIphiIphi = cms.string("full5x5_sigmaIphiIphi"),
    probe_ele_phoESEffSigmaRR = cms.InputTag("photonIDValueMapProducer:phoESEffSigmaRR"),
    probe_ele_full5x5_e1x5 = cms.string("full5x5_e1x5"),
    probe_ele_full5x5_e5x5 = cms.string("full5x5_e5x5"),
    probe_ele_full5x5_r9 = cms.string("full5x5_r9"),
    probe_ele_full5x5_e2x5Max = cms.string("full5x5_e2x5Max"),
    probe_ele_full5x5_hcalOverEcal = cms.string("full5x5_hcalOverEcal"),
    probe_ele_full5x5_hcalOverEcalBc = cms.string("full5x5_hcalOverEcalBc"),


#    probe_scl_eta = cms.string("superCluster.eta"),
    probe_scl_phi = cms.string("superCluster.phi"),
    probe_scl_energy = cms.string("superCluster.energy"),
    probe_scl_rawEnergy = cms.string("superCluster.rawEnergy"),
    probe_scl_phiWidth = cms.string("superCluster.phiWidth"),
    probe_scl_etaWidth = cms.string("superCluster.etaWidth"),
    probe_scl_clustersSize = cms.string("superCluster.clustersSize"),
    probe_scl_preshowerEnergy = cms.string("superCluster.preshowerEnergy"),
    probe_ele_psEoverEraw = cms.string("superCluster.preshowerEnergy/superCluster.rawEnergy"),

    probe_ele_fbrem = cms.string("fbrem"),
    probe_ele_superClusterFbrem = cms.string("superClusterFbrem"),
    probe_ele_numberOfBrems = cms.string("numberOfBrems"),
    probe_ele_classification = cms.string("classification"),
    probe_ele_gsf_normalizedChi2 = cms.string("gsfTrack.normalizedChi2"),  
    probe_ele_gsf_trackerLayersWithMeasurement = cms.string("gsfTrack.hitPattern.trackerLayersWithMeasurement"),  
    probe_ele_gsf_found = cms.string("gsfTrack.found"),
    probe_ele_gsf_lost = cms.string("gsfTrack.lost"),
    probe_ele_gsf_stripLayersWithMeasurement = cms.string("gsfTrack.hitPattern.stripLayersWithMeasurement"),  
    probe_ele_gsf_pixelLayersWithMeasurement = cms.string("gsfTrack.hitPattern.pixelLayersWithMeasurement"),  


    probe_ele_kf_normalizedChi2 = cms.InputTag("eleVarHelper:ele-kf-normalizedChi2"),  
    probe_ele_kf_trackerLayersWithMeasurement = cms.InputTag("eleVarHelper:ele-kf-trackerLayersWithMeasurement"),  

    probe_ele_ecalEnergyError = cms.string("ecalEnergyError"),
    probe_ele_trackMomentumError = cms.string("trackMomentumError"),
    probe_ele_p4Error_combination = cms.InputTag("eleVarHelper:ele-p4Error-combination"),
    probe_ele_p4Error_PFcombination = cms.InputTag("eleVarHelper:ele-p4Error-PFcombination"),
    

    probe_ele_convDist = cms.string("convDist"),
    probe_ele_convDcot = cms.string("convDcot"),
    probe_ele_convRadius = cms.string("convRadius"),
    probe_ele_missing_inner_hits = cms.InputTag("eleVarHelper:ele-missing-inner-hits"),  
    probe_ele_conv_vertex_fit_prob = cms.InputTag("eleVarHelper:ele-conv-vertex-fit-prob"),

    #probe_ele_electronEcalPFClusterIso = cms.InputTag("electronEcalPFClusterIsolationProducer"),
    #probe_ele_electronHcalPFClusterIso = cms.InputTag("electronHcalPFClusterIsolationProducer"),
    probe_ele_dr03EcalRecHitSumEt = cms.string("dr03EcalRecHitSumEt"),
    probe_ele_dr03HcalTowerSumEt = cms.string("dr03HcalTowerSumEt"),
    probe_ele_dr03TkSumPt = cms.string("dr03TkSumPt"),


    probe_ele_pfIso_sumChargedHadronPt = cms.string("pfIsolationVariables.sumChargedHadronPt"),
    probe_ele_pfIso_sumChargedParticlePt = cms.string("pfIsolationVariables.sumChargedParticlePt"),
    probe_ele_pfIso_sumNeutralHadronEt = cms.string("pfIsolationVariables.sumNeutralHadronEt"),
    probe_ele_pfIso_sumPhotonEt = cms.string("pfIsolationVariables.sumPhotonEt"),
    probe_ele_pfIso_sumPUPt = cms.string("pfIsolationVariables.sumPUPt"),


    probe_ele_isEB = cms.string("isEB"),
    probe_ele_isEE = cms.string("isEE"),
    probe_ele_isGap = cms.string("isGap"),
    probe_ele_isEBEEGap = cms.string("isEBEEGap"),
    probe_ele_isEBEtaGap = cms.string("isEBEtaGap"),
    probe_ele_isEBPhiGap = cms.string("isEBPhiGap"),
    probe_ele_isEEDeeGap = cms.string("isEEDeeGap"),
    probe_ele_isEERingGap = cms.string("isEERingGap"),
    probe_ele_ecalDrivenSeed = cms.string("ecalDrivenSeed"),
    probe_ele_trackerDrivenSeed = cms.string("trackerDrivenSeed"),

    probe_ele_fixedGridRhoFastjetAll = cms.InputTag("eleVarHelper:ele-fixedGridRhoFastjetAll"),
    probe_ele_fixedGridRhoFastjetCentralNeutral = cms.InputTag("eleVarHelper:ele-fixedGridRhoFastjetCentralNeutral"),

    #sip
    probe_Ele_ip            = cms.InputTag("eleVarHelper:ip"),
    probe_Ele_iperror       = cms.InputTag("eleVarHelper:iperror"),
    probe_Ele_sip           = cms.InputTag("eleVarHelper:sip"),
                                 
    #isolation
    probe_Ele_chIso         = cms.string("pfIsolationVariables().sumChargedHadronPt"),
    probe_Ele_phoIso        = cms.string("pfIsolationVariables().sumPhotonEt"),
    probe_Ele_neuIso        = cms.string("pfIsolationVariables().sumNeutralHadronEt"),
    probe_Ele_effarea       = cms.InputTag("eleVarHelper:effarea"),
    probe_Ele_rho           = cms.InputTag("eleVarHelper:rho"),
    probe_Ele_ed            = cms.InputTag("eleVarHelper:enedens"),
    probe_Ele_iso           = cms.InputTag("eleVarHelper:iso"),
    probe_Ele_correctediso  = cms.InputTag("eleVarHelper:correctediso"),
    probe_Ele_trkIso        = cms.string("dr03IsolationVariables().tkSumPt"),

    #FSR photon
    probe_FSRcorr                = cms.InputTag("eleVarHelper:FSRcorr"),
    probe_FSRphoton_px           = cms.InputTag("eleVarHelper:FSRpx"),
    probe_FSRphoton_py           = cms.InputTag("eleVarHelper:FSRpy"),
    probe_FSRphoton_pz           = cms.InputTag("eleVarHelper:FSRpz"),
    probe_FSRphoton_e            = cms.InputTag("eleVarHelper:FSRe"),
    )

PhoProbeVariablesToStore = cms.PSet(
    probe_Pho_eta    = cms.string("eta"),
    probe_Pho_abseta = cms.string("abs(eta)"),
    probe_Pho_et     = cms.string("et"),
    probe_Pho_e      = cms.string("energy"),

## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et     = cms.string("superCluster.energy*sin(superCluster.position.theta)"),    
    probe_sc_eta    = cms.string("-log(tan(superCluster.position.theta/2))"),
    probe_sc_abseta = cms.string("abs(-log(tan(superCluster.position.theta/2)))"),


#id based
    probe_Pho_full5x5x_r9   = cms.string("full5x5_r9"),
    probe_Pho_r9            = cms.string("r9"),
    probe_Pho_sieie         = cms.string("full5x5_sigmaIetaIeta"),
    probe_Pho_sieip         = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIPhi"),
    probe_Pho_ESsigma       = cms.InputTag("photonIDValueMapProducer:phoESEffSigmaRR"),
    probe_Pho_hoe           = cms.string("hadronicOverEm"),

#iso
    probe_Pho_chIso    = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    probe_Pho_neuIso   = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    probe_Pho_phoIso   = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    probe_Pho_chWorIso = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"), 

#pho mva
    probe_mva          = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values"),
)




TagVariablesToStore = cms.PSet(
    Ele_eta    = cms.string("eta"),
    Ele_phi    = cms.string("phi"),
    Ele_abseta = cms.string("abs(eta)"),
    Ele_pt     = cms.string("pt"),
    Ele_et     = cms.string("et"),
    Ele_e      = cms.string("energy"),
    Ele_q      = cms.string("charge"),
    
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    sc_eta    = cms.string("-log(tan(superClusterPosition.theta/2))"),
    sc_abseta = cms.string("abs(-log(tan(superCluster.position.theta/2)))"),
 #   
    Ele_mHits         = cms.InputTag("eleVarHelper:missinghits"),
    Ele_dz            = cms.InputTag("eleVarHelper:dz"),
    Ele_dxy           = cms.InputTag("eleVarHelper:dxy"),
    Ele_nonTrigMVA    = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    Ele_trigMVA       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),

    Ele_dEtaIn        = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    Ele_dPhiIn        = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    Ele_dEtaSeed      = cms.string("deltaEtaSuperClusterTrackAtVtx+log(tan(superCluster.position.theta/2))-log(tan(superCluster.seed.position.theta/2))"),
    Ele_sieie         = cms.string("sigmaIetaIeta"),
    Ele_e1x5          = cms.string("e1x5"),
    Ele_e2x5          = cms.string("e2x5Max"),
    Ele_e5x5          = cms.string("e5x5"),
    Ele_r9            = cms.string("r9"),
    Ele_r9_5x5        = cms.string("full5x5_r9"),
    Ele_sieie_5x5     = cms.string("full5x5_sigmaIetaIeta"),
    Ele_hoe           = cms.string("hadronicOverEm"),
    Ele_ooemoop       = cms.string("(1.0/ecalEnergy - eSuperClusterOverP/ecalEnergy)"),



    #FSR photon
    FSRphoton_px           = cms.InputTag("eleVarHelper:FSRpx"),
    FSRphoton_py           = cms.InputTag("eleVarHelper:FSRpy"),
    FSRphoton_pz           = cms.InputTag("eleVarHelper:FSRpz"),
    FSRphoton_e            = cms.InputTag("eleVarHelper:FSRe"),

    )

CommonStuffForGsfElectronProbe = cms.PSet(
    addEventVariablesInfo   =  cms.bool(True),

    variables        = cms.PSet(EleProbeVariablesToStore),
    pairVariables    =  cms.PSet(ZVariablesToStore),
    tagVariables   =  cms.PSet(TagVariablesToStore),

    ignoreExceptions =  cms.bool (True),
    addRunLumiInfo   =  cms.bool (True),
    pileupInfoTag    = cms.InputTag("slimmedAddPileupInfo"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot         = cms.InputTag("offlineBeamSpot"),
    pfMet            = cms.InputTag("slimmedMETsPuppi"),
    rho              = cms.InputTag("fixedGridRhoFastjetAll"),
    #    pfMet = cms.InputTag("slimmedMETsNoHF"),
    #    rho = cms.InputTag("fixedGridRhoAll"),

    pairFlags     =  cms.PSet(
        mass60to120 = cms.string("60 < mass < 120")
        ),
    tagFlags       =  cms.PSet(),    

    )

CommonStuffForPhotonProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForPhotonProbe.variables = cms.PSet(PhoProbeVariablesToStore)

CommonStuffForSuperClusterProbe = CommonStuffForGsfElectronProbe.clone()
CommonStuffForSuperClusterProbe.variables = cms.PSet(SCProbeVariablesToStore)

mcTruthCommonStuff = cms.PSet(
    isMC = cms.bool(True),
    #tagMatches = cms.InputTag("McMatchTag"),
    motherPdgId = cms.vint32(),
    #motherPdgId = cms.vint32(22,23),
    #motherPdgId = cms.vint32(443), # JPsi
    #motherPdgId = cms.vint32(553), # Yupsilon
    #makeMCUnbiasTree = cms.bool(False),
    #checkMotherInUnbiasEff = cms.bool(False),
    genParticles = cms.InputTag("prunedGenParticles"),
    useTauDecays = cms.bool(False),
    checkCharge = cms.bool(False),
    pdgId = cms.int32(11),
    mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_abseta = cms.string("abs(eta)"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        ),
    mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),      
    )



def setupTnPVariablesForAOD():
    mcTruthCommonStuff. genParticles                 = cms.InputTag("genParticles")

    CommonStuffForSuperClusterProbe.pileupInfoTag    = cms.InputTag("addPileupInfo")
    CommonStuffForSuperClusterProbe.vertexCollection = cms.InputTag("offlinePrimaryVerticesWithBS")
    CommonStuffForSuperClusterProbe.pfMet            = cms.InputTag("pfMet")

    CommonStuffForGsfElectronProbe.pileupInfoTag     = cms.InputTag("addPileupInfo")
    CommonStuffForGsfElectronProbe.vertexCollection  = cms.InputTag("offlinePrimaryVerticesWithBS")
    CommonStuffForGsfElectronProbe.pfMet             = cms.InputTag("pfMet")

    CommonStuffForPhotonProbe.pileupInfoTag          = cms.InputTag("addPileupInfo")
    CommonStuffForPhotonProbe.vertexCollection       = cms.InputTag("offlinePrimaryVerticesWithBS")
    CommonStuffForPhotonProbe.pfMet                  = cms.InputTag("pfMet")
