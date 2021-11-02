import FWCore.ParameterSet.Config as cms

process = cms.Process("NTuple")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/10000/05134D30-E44E-B649-890A-492D4E85AF3C.root')
)
process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.CondDBTauConnection = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06
    )
)

process.METSignificanceParams = cms.PSet(
    dRMatch = cms.double(0.4),
    jetThreshold = cms.double(15),
    jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
    jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
    pjpar = cms.vdouble(-0.2586, 0.6173)
)

process.METSignificance_params = cms.PSet(
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EB_PhiResPar = cms.vdouble(0.00502),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EE_PhiResPar = cms.vdouble(0.02511),
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    HB_PhiResPar = cms.vdouble(0.02511),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    HE_PhiResPar = cms.vdouble(0.02511),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    HF_PhiResPar = cms.vdouble(0.05022),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    HO_PhiResPar = cms.vdouble(0.02511),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType4 = cms.vdouble(0.042, 0.1, 0.0),
    PF_EtResType5 = cms.vdouble(0.41, 0.52, 0.25),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    PF_PhiResType4 = cms.vdouble(0.0028, 0.0, 0.0022),
    PF_PhiResType5 = cms.vdouble(0.1, 0.1, 0.13),
    PF_PhiResType6 = cms.vdouble(0.02511),
    PF_PhiResType7 = cms.vdouble(0.02511),
    jdphi0 = cms.vdouble(
        0.034, 0.034, 0.034, 0.034, 0.032, 
        0.031, 0.028, 0.027, 0.027, 0.027
    ),
    jdphi1 = cms.vdouble(
        0.034, 0.035, 0.035, 0.035, 0.035, 
        0.034, 0.031, 0.03, 0.029, 0.027
    ),
    jdphi2 = cms.vdouble(
        0.04, 0.04, 0.04, 0.04, 0.04, 
        0.038, 0.036, 0.035, 0.034, 0.033
    ),
    jdphi3 = cms.vdouble(
        0.042, 0.043, 0.044, 0.043, 0.041, 
        0.039, 0.039, 0.036, 0.034, 0.031
    ),
    jdphi4 = cms.vdouble(
        0.042, 0.042, 0.043, 0.042, 0.038, 
        0.036, 0.036, 0.033, 0.031, 0.031
    ),
    jdphi5 = cms.vdouble(
        0.069, 0.069, 0.064, 0.058, 0.053, 
        0.049, 0.049, 0.043, 0.039, 0.04
    ),
    jdphi6 = cms.vdouble(
        0.084, 0.08, 0.072, 0.065, 0.066, 
        0.06, 0.051, 0.049, 0.045, 0.045
    ),
    jdphi7 = cms.vdouble(
        0.077, 0.072, 0.059, 0.05, 0.045, 
        0.042, 0.039, 0.039, 0.037, 0.031
    ),
    jdphi8 = cms.vdouble(
        0.059, 0.057, 0.051, 0.044, 0.038, 
        0.035, 0.037, 0.032, 0.028, 0.028
    ),
    jdphi9 = cms.vdouble(
        0.062, 0.059, 0.053, 0.047, 0.042, 
        0.045, 0.036, 0.032, 0.034, 0.044
    ),
    jdpt0 = cms.vdouble(
        0.749, 0.829, 1.099, 1.355, 1.584, 
        1.807, 2.035, 2.217, 2.378, 2.591
    ),
    jdpt1 = cms.vdouble(
        0.718, 0.813, 1.133, 1.384, 1.588, 
        1.841, 2.115, 2.379, 2.508, 2.772
    ),
    jdpt2 = cms.vdouble(
        0.841, 0.937, 1.316, 1.605, 1.919, 
        2.295, 2.562, 2.722, 2.943, 3.293
    ),
    jdpt3 = cms.vdouble(
        0.929, 1.04, 1.46, 1.74, 2.042, 
        2.289, 2.639, 2.837, 2.946, 2.971
    ),
    jdpt4 = cms.vdouble(
        0.85, 0.961, 1.337, 1.593, 1.854, 
        2.005, 2.209, 2.533, 2.812, 3.047
    ),
    jdpt5 = cms.vdouble(
        1.049, 1.149, 1.607, 1.869, 2.012, 
        2.219, 2.289, 2.412, 2.695, 2.865
    ),
    jdpt6 = cms.vdouble(
        1.213, 1.298, 1.716, 2.015, 2.191, 
        2.612, 2.863, 2.879, 2.925, 2.902
    ),
    jdpt7 = cms.vdouble(
        1.094, 1.139, 1.436, 1.672, 1.831, 
        2.05, 2.267, 2.549, 2.785, 2.86
    ),
    jdpt8 = cms.vdouble(
        0.889, 0.939, 1.166, 1.365, 1.553, 
        1.805, 2.06, 2.22, 2.268, 2.247
    ),
    jdpt9 = cms.vdouble(
        0.843, 0.885, 1.245, 1.665, 1.944, 
        1.981, 1.972, 2.875, 3.923, 7.51
    ),
    ptresolthreshold = cms.double(10.0),
    resolutionsAlgo = cms.string('AK5PF'),
    resolutionsEra = cms.string('Spring10')
)

process.calibratedEgammaPatSettings = cms.PSet(
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True)
)

process.calibratedEgammaSettings = cms.PSet(
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True)
)

process.ecalTrkCombinationRegression = cms.PSet(
    ecalTrkRegressionConfig = cms.PSet(
        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        lowEtHighEtBoundary = cms.double(50.0),
        rangeMaxHighEt = cms.double(3.0),
        rangeMaxLowEt = cms.double(3.0),
        rangeMinHighEt = cms.double(-1.0),
        rangeMinLowEt = cms.double(-1.0)
    ),
    ecalTrkRegressionUncertConfig = cms.PSet(
        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        lowEtHighEtBoundary = cms.double(50.0),
        rangeMaxHighEt = cms.double(0.5),
        rangeMaxLowEt = cms.double(0.5),
        rangeMinHighEt = cms.double(0.0002),
        rangeMinLowEt = cms.double(0.0002)
    ),
    maxEPDiffInSigmaForComb = cms.double(15.0),
    maxEcalEnergyForComb = cms.double(200.0),
    maxRelTrkMomErrForComb = cms.double(10.0),
    minEOverPForComb = cms.double(0.025)
)

process.isolationSumsCalculator = cms.PSet(
    ComponentName = cms.string('isolationSumsCalculator'),
    EcalRecHitEtaSliceA_Barrel = cms.double(2.5),
    EcalRecHitEtaSliceA_Endcap = cms.double(2.5),
    EcalRecHitEtaSliceB_Barrel = cms.double(2.5),
    EcalRecHitEtaSliceB_Endcap = cms.double(2.5),
    EcalRecHitInnerRadiusA_Barrel = cms.double(3.5),
    EcalRecHitInnerRadiusA_Endcap = cms.double(3.5),
    EcalRecHitInnerRadiusB_Barrel = cms.double(3.5),
    EcalRecHitInnerRadiusB_Endcap = cms.double(3.5),
    EcalRecHitOuterRadiusA_Barrel = cms.double(0.4),
    EcalRecHitOuterRadiusA_Endcap = cms.double(0.4),
    EcalRecHitOuterRadiusB_Barrel = cms.double(0.3),
    EcalRecHitOuterRadiusB_Endcap = cms.double(0.3),
    EcalRecHitThreshEA_Barrel = cms.double(0.095),
    EcalRecHitThreshEA_Endcap = cms.double(0.0),
    EcalRecHitThreshEB_Barrel = cms.double(0.095),
    EcalRecHitThreshEB_Endcap = cms.double(0.0),
    EcalRecHitThreshEtA_Barrel = cms.double(0.0),
    EcalRecHitThreshEtA_Endcap = cms.double(0.11),
    EcalRecHitThreshEtB_Barrel = cms.double(0.0),
    EcalRecHitThreshEtB_Endcap = cms.double(0.11),
    HcalDepth1TowerInnerRadiusA_Barrel = cms.double(0.15),
    HcalDepth1TowerInnerRadiusA_Endcap = cms.double(0.15),
    HcalDepth1TowerInnerRadiusB_Barrel = cms.double(0.15),
    HcalDepth1TowerInnerRadiusB_Endcap = cms.double(0.15),
    HcalDepth1TowerOuterRadiusA_Barrel = cms.double(0.4),
    HcalDepth1TowerOuterRadiusA_Endcap = cms.double(0.4),
    HcalDepth1TowerOuterRadiusB_Barrel = cms.double(0.3),
    HcalDepth1TowerOuterRadiusB_Endcap = cms.double(0.3),
    HcalDepth1TowerThreshEA_Barrel = cms.double(0.0),
    HcalDepth1TowerThreshEA_Endcap = cms.double(0.0),
    HcalDepth1TowerThreshEB_Barrel = cms.double(0.0),
    HcalDepth1TowerThreshEB_Endcap = cms.double(0.0),
    HcalDepth2TowerInnerRadiusA_Barrel = cms.double(0.15),
    HcalDepth2TowerInnerRadiusA_Endcap = cms.double(0.15),
    HcalDepth2TowerInnerRadiusB_Barrel = cms.double(0.15),
    HcalDepth2TowerInnerRadiusB_Endcap = cms.double(0.15),
    HcalDepth2TowerOuterRadiusA_Barrel = cms.double(0.4),
    HcalDepth2TowerOuterRadiusA_Endcap = cms.double(0.4),
    HcalDepth2TowerOuterRadiusB_Barrel = cms.double(0.3),
    HcalDepth2TowerOuterRadiusB_Endcap = cms.double(0.3),
    HcalDepth2TowerThreshEA_Barrel = cms.double(0.0),
    HcalDepth2TowerThreshEA_Endcap = cms.double(0.0),
    HcalDepth2TowerThreshEB_Barrel = cms.double(0.0),
    HcalDepth2TowerThreshEB_Endcap = cms.double(0.0),
    HcalRecHitCollection = cms.InputTag("towerMaker"),
    HcalTowerInnerRadiusA_Barrel = cms.double(0.15),
    HcalTowerInnerRadiusA_Endcap = cms.double(0.15),
    HcalTowerInnerRadiusB_Barrel = cms.double(0.15),
    HcalTowerInnerRadiusB_Endcap = cms.double(0.15),
    HcalTowerOuterRadiusA_Barrel = cms.double(0.4),
    HcalTowerOuterRadiusA_Endcap = cms.double(0.4),
    HcalTowerOuterRadiusB_Barrel = cms.double(0.3),
    HcalTowerOuterRadiusB_Endcap = cms.double(0.3),
    HcalTowerThreshEA_Barrel = cms.double(0.0),
    HcalTowerThreshEA_Endcap = cms.double(0.0),
    HcalTowerThreshEB_Barrel = cms.double(0.0),
    HcalTowerThreshEB_Endcap = cms.double(0.0),
    TrackConeInnerRadiusA_Barrel = cms.double(0.04),
    TrackConeInnerRadiusA_Endcap = cms.double(0.04),
    TrackConeInnerRadiusB_Barrel = cms.double(0.04),
    TrackConeInnerRadiusB_Endcap = cms.double(0.04),
    TrackConeOuterRadiusA_Barrel = cms.double(0.4),
    TrackConeOuterRadiusA_Endcap = cms.double(0.4),
    TrackConeOuterRadiusB_Barrel = cms.double(0.3),
    TrackConeOuterRadiusB_Endcap = cms.double(0.3),
    barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    beamSpotProducer = cms.InputTag("offlineBeamSpot"),
    endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    isolationtrackEtaSliceA_Barrel = cms.double(0.015),
    isolationtrackEtaSliceA_Endcap = cms.double(0.015),
    isolationtrackEtaSliceB_Barrel = cms.double(0.015),
    isolationtrackEtaSliceB_Endcap = cms.double(0.015),
    isolationtrackThresholdA_Barrel = cms.double(0.0),
    isolationtrackThresholdA_Endcap = cms.double(0.0),
    isolationtrackThresholdB_Barrel = cms.double(0.0),
    isolationtrackThresholdB_Endcap = cms.double(0.0),
    longImpactParameterA_Barrel = cms.double(0.2),
    longImpactParameterA_Endcap = cms.double(0.2),
    longImpactParameterB_Barrel = cms.double(0.2),
    longImpactParameterB_Endcap = cms.double(0.2),
    moduleEtaBoundary = cms.vdouble(
        0.0, 0.02, 0.43, 0.46, 0.78, 
        0.81, 1.13, 1.15, 1.45, 1.58
    ),
    modulePhiBoundary = cms.double(0.0087),
    trackProducer = cms.InputTag("generalTracks"),
    transImpactParameterA_Barrel = cms.double(0.1),
    transImpactParameterA_Endcap = cms.double(0.1),
    transImpactParameterB_Barrel = cms.double(0.1),
    transImpactParameterB_Endcap = cms.double(0.1),
    useNumCrystals = cms.bool(True),
    vetoClustered = cms.bool(False)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mvaPhoID_RunIIFall17_v1p1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'abs(superCluster.eta) <  1.479', 
        'abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('PhotonMVAEstimator'),
    mvaTag = cms.string('RunIIFall17v1p1'),
    nCategories = cms.int32(2),
    variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V1.weights.xml.gz', 
        'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V1.weights.xml.gz'
    )
)

process.mvaPhoID_RunIIFall17_v1p1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Categories"),
        mvaCuts = cms.vdouble(0.67, 0.54),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1p1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_RunIIFall17_v1p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Categories"),
        mvaCuts = cms.vdouble(0.27, 0.14),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v1p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v1p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_RunIIFall17_v2_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'abs(superCluster.eta) <  1.479', 
        'abs(superCluster.eta) >= 1.479'
    ),
    mvaName = cms.string('PhotonMVAEstimator'),
    mvaTag = cms.string('RunIIFall17v2'),
    nCategories = cms.int32(2),
    variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V2.weights.xml.gz', 
        'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V2.weights.xml.gz'
    )
)

process.mvaPhoID_RunIIFall17_v2_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v2Categories"),
        mvaCuts = cms.vdouble(0.42, 0.14),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v2Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v2-wp80'),
    isPOGApproved = cms.bool(True)
)

process.mvaPhoID_RunIIFall17_v2_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v2Categories"),
        mvaCuts = cms.vdouble(-0.02, -0.26),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRunIIFall17v2Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-RunIIFall17-v2-wp90'),
    isPOGApproved = cms.bool(True)
)

process.mvaPhoID_Spring16_nonTrig_V1_producer_config = cms.PSet(
    categoryCuts = cms.vstring(
        'abs(superCluster.eta) <  1.479', 
        'abs(superCluster.eta) >= 1.479'
    ),
    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
    mvaName = cms.string('PhotonMVAEstimator'),
    mvaTag = cms.string('Run2Spring16NonTrigV1'),
    nCategories = cms.int32(2),
    phoIsoCutoff = cms.double(2.5),
    phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
    variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesSpring16.txt'),
    weightFileNames = cms.vstring(
        'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EB_V1.weights.xml.gz', 
        'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EE_V1.weights.xml.gz'
    )
)

process.options = cms.untracked.PSet(
    fileMode = cms.untracked.string('NOMERGE'),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(4)
)

process.pset = cms.PSet(
    etaMax = cms.double(-2.901376),
    etaMin = cms.double(-5.2),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('egammaHFMinus'),
    px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
    py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
    type = cms.int32(7),
    varType = cms.int32(0)
)

process.trkIsol03CfgV2 = cms.PSet(
    barrelCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.3),
        maxDZ = cms.double(0.1),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    ),
    endcapCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.3),
        maxDZ = cms.double(0.5),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    )
)

process.trkIsol04CfgV2 = cms.PSet(
    barrelCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.4),
        maxDZ = cms.double(0.1),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    ),
    endcapCuts = cms.PSet(
        algosToReject = cms.vstring(),
        allowedQualities = cms.vstring(),
        maxDPtPt = cms.double(0.1),
        maxDR = cms.double(0.4),
        maxDZ = cms.double(0.5),
        minDEta = cms.double(0.005),
        minDR = cms.double(0.0),
        minHits = cms.int32(8),
        minPixelHits = cms.int32(1),
        minPt = cms.double(1.0)
    )
)

process.variableJTAPars = cms.PSet(
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5)
)

process.vertexCutsBlock = cms.PSet(
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    )
)

process.vertexRecoBlock = cms.PSet(
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    )
)

process.vertexSelectionBlock = cms.PSet(
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)

process.vertexTrackSelectionBlock = cms.PSet(
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    )
)

process.mvaConfigsForPhoProducer = cms.VPSet(
    cms.PSet(
        categoryCuts = cms.vstring(
            'abs(superCluster.eta) <  1.479', 
            'abs(superCluster.eta) >= 1.479'
        ),
        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
        mvaName = cms.string('PhotonMVAEstimator'),
        mvaTag = cms.string('Run2Spring16NonTrigV1'),
        nCategories = cms.int32(2),
        phoIsoCutoff = cms.double(2.5),
        phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
        variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesSpring16.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EB_V1.weights.xml.gz', 
            'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EE_V1.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'abs(superCluster.eta) <  1.479', 
            'abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('PhotonMVAEstimator'),
        mvaTag = cms.string('RunIIFall17v1p1'),
        nCategories = cms.int32(2),
        variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V1.weights.xml.gz', 
            'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V1.weights.xml.gz'
        )
    ), 
    cms.PSet(
        categoryCuts = cms.vstring(
            'abs(superCluster.eta) <  1.479', 
            'abs(superCluster.eta) >= 1.479'
        ),
        mvaName = cms.string('PhotonMVAEstimator'),
        mvaTag = cms.string('RunIIFall17v2'),
        nCategories = cms.int32(2),
        variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
        weightFileNames = cms.vstring(
            'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V2.weights.xml.gz', 
            'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V2.weights.xml.gz'
        )
    )
)

process.patMultPhiCorrParams_T0pcT1SmearTxy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1SmearTxy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1T2Txy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1T2Txy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1Txy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcT1Txy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcTxy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T0pcTxy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1SmearTxy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1SmearTxy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1T2SmearTxy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1T2SmearTxy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1T2Txy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1T2Txy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1Txy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_T1Txy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_Txy_25ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.patMultPhiCorrParams_Txy_50ns = cms.VPSet(
    cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
        py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    )
)

process.puppiCentral = cms.VPSet(cms.PSet(
    algoId = cms.int32(5),
    applyLowPUCorr = cms.bool(True),
    combOpt = cms.int32(0),
    cone = cms.double(0.4),
    rmsPtMin = cms.double(0.1),
    rmsScaleFactor = cms.double(1.0),
    useCharged = cms.bool(True)
))

process.puppiForward = cms.VPSet(cms.PSet(
    algoId = cms.int32(5),
    applyLowPUCorr = cms.bool(True),
    combOpt = cms.int32(0),
    cone = cms.double(0.4),
    rmsPtMin = cms.double(0.5),
    rmsScaleFactor = cms.double(1.0),
    useCharged = cms.bool(False)
))

process.ak10PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1FastL2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1FastL2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1L2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1L2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak1PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1FastL2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1FastL2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1L2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1L2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak2PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1FastL2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1FastL2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1L2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1L2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak3PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1FastL2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1FastL2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1L2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1L2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak4CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1FastL2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1FastL2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1L2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1L2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4CaloL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTJetsL1FastL2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1FastL2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1FastL2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1FastL2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1L2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1L2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1L2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1L2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4L1JPTFastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4L1JPTFastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTFastjetCorrector = cms.EDProducer("L1JPTOffsetCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.InputTag("ak4CaloL1FastjetCorrector")
)


process.ak4L1JPTOffsetCorrector = cms.EDProducer("L1JPTOffsetCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.InputTag("ak4CaloL1OffsetCorrector")
)


process.ak4PFCHSL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFCHSL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1FastL2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1FastL2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1L2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1L2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFPuppiL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFPuppiL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFPuppiL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2Relative')
)


process.ak4PFPuppiL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L3Absolute')
)


process.ak4PFPuppiResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2L3Residual')
)


process.ak4PFResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackJetsL2L3 = cms.EDProducer("TrackJetCorrectionProducer",
    correctors = cms.vstring('ak4TrackL2L3'),
    src = cms.InputTag("ak4TrackJets")
)


process.ak4TrackL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4TrackL2RelativeCorrector", "ak4TrackL3AbsoluteCorrector")
)


process.ak4TrackL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1FastL2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1FastL2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1L2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1L2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak6PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1FastL2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1FastL2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1L2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1L2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak7CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1FastL2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1FastL2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1L2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1L2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1FastL2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1FastL2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1L2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1L2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak8PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1FastL2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1FastL2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1L2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1L2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak9PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1FastL2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1FastL2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1L2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1L2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.basicJetsForMet = cms.EDProducer("PATJetCleanerForType1MET",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("patJetsReapplyJEC"),
    type1JetPtThreshold = cms.double(15.0)
)


process.basicJetsForMetPuppi = cms.EDProducer("PATJetCleanerForType1MET",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("patJetsReapplyJECPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.calibratedElectrons = cms.EDProducer("CalibratedElectronProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    epCombConfig = cms.PSet(
        ecalTrkRegressionConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMaxHighEt = cms.double(3.0),
            rangeMaxLowEt = cms.double(3.0),
            rangeMinHighEt = cms.double(-1.0),
            rangeMinLowEt = cms.double(-1.0)
        ),
        ecalTrkRegressionUncertConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMaxHighEt = cms.double(0.5),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMinLowEt = cms.double(0.0002)
        ),
        maxEPDiffInSigmaForComb = cms.double(15.0),
        maxEcalEnergyForComb = cms.double(200.0),
        maxRelTrkMomErrForComb = cms.double(10.0),
        minEOverPForComb = cms.double(0.025)
    ),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("gedGsfElectrons")
)


process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain_v2'),
    epCombConfig = cms.PSet(
        ecalTrkRegressionConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMaxHighEt = cms.double(3.0),
            rangeMaxLowEt = cms.double(3.0),
            rangeMinHighEt = cms.double(-1.0),
            rangeMinLowEt = cms.double(-1.0)
        ),
        ecalTrkRegressionUncertConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMaxHighEt = cms.double(0.5),
            rangeMaxLowEt = cms.double(0.5),
            rangeMinHighEt = cms.double(0.0002),
            rangeMinLowEt = cms.double(0.0002)
        ),
        maxEPDiffInSigmaForComb = cms.double(15.0),
        maxEcalEnergyForComb = cms.double(200.0),
        maxRelTrkMomErrForComb = cms.double(10.0),
        minEOverPForComb = cms.double(0.025)
    ),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(False),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("updatedElectrons")
)


process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_Step2Closure_CoarseEtaR9Gain_v2'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(False),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("updatedPhotons")
)


process.calibratedPhotons = cms.EDProducer("CalibratedPhotonProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("gedPhotons")
)


process.caloMetT1 = cms.EDProducer("CorrectedCaloMETProducer",
    src = cms.InputTag("caloMetM"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrCaloMetType1","type1"))
)


process.caloMetT1T2 = cms.EDProducer("CorrectedCaloMETProducer",
    src = cms.InputTag("caloMetM"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrCaloMetType1","type1"), cms.InputTag("corrCaloMetType2"))
)


process.cleanPatElectrons = cms.EDProducer("PATElectronCleaner",
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatMuons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("selectedPatElectrons")
)


process.cleanPatJets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatMuons")
        ),
        photons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatPhotons")
        ),
        taus = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatTaus")
        ),
        tkIsoElectrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            preselection = cms.string('pt > 10 && trackIso < 3'),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatElectrons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("selectedPatJets")
)


process.cleanPatMuons = cms.EDProducer("PATMuonCleaner",
    checkOverlaps = cms.PSet(

    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("selectedPatMuons")
)


process.cleanPatPhotons = cms.EDProducer("PATPhotonCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('bySuperClusterSeed'),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatElectrons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("selectedPatPhotons")
)


process.cleanPatTaus = cms.EDProducer("PATTauCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.3),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("cleanPatMuons")
        )
    ),
    finalCut = cms.string('pt > 18. & abs(eta) < 2.3'),
    preselection = cms.string('tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 & tauID("againstMuonTight3") > 0.5 & tauID("againstElectronVLooseMVA6") > 0.5'),
    src = cms.InputTag("selectedPatTaus")
)


process.cleanedPatJets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedMuons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("jetSelectorForMet")
)


process.cleanedPatJetsPuppi = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedMuons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("jetSelectorForMetPuppi")
)


process.corrCaloMetType1 = cms.EDProducer("CaloJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4CaloL2L3Corrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    src = cms.InputTag("ak4CaloJets"),
    srcMET = cms.InputTag("caloMetM"),
    type1JetPtThreshold = cms.double(20.0)
)


process.corrCaloMetType2 = cms.EDProducer("Type2CorrectionProducer",
    srcUnclEnergySums = cms.VInputTag(cms.InputTag("corrCaloMetType1","type2"), cms.InputTag("muCaloMetCorr")),
    type2CorrFormula = cms.string('A + B*TMath::Exp(-C*x)'),
    type2CorrParameter = cms.PSet(
        A = cms.double(2.0),
        B = cms.double(1.3),
        C = cms.double(0.1)
    )
)


process.corrPfMetType1 = cms.EDProducer("PFJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelRes = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    offsetCorrLabel = cms.InputTag("ak4PFCHSL1FastjetCorrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("ak4PFJetsCHS"),
    type1JetPtThreshold = cms.double(15.0)
)


process.corrPfMetType1Puppi = cms.EDProducer("PFJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4PFPuppiL1FastL2L3Corrector"),
    jetCorrLabelRes = cms.InputTag("ak4PFPuppiL1FastL2L3ResidualCorrector"),
    offsetCorrLabel = cms.InputTag("ak4PFPuppiL1FastjetCorrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("ak4PFJetsPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.corrPfMetType2 = cms.EDProducer("Type2CorrectionProducer",
    srcUnclEnergySums = cms.VInputTag(cms.InputTag("corrPfMetType1","type2"), cms.InputTag("corrPfMetType1","offset"), cms.InputTag("pfCandMETcorr")),
    type2CorrFormula = cms.string('A'),
    type2CorrParameter = cms.PSet(
        A = cms.double(1.4)
    )
)


process.deepTau2017v2p1 = cms.EDProducer("DeepTauId",
    VSeWP = cms.PSet(
        Loose = cms.string('0.6815435'),
        Medium = cms.string('0.8847544'),
        Tight = cms.string('0.9675541'),
        VLoose = cms.string('0.362813'),
        VTight = cms.string('0.9859251'),
        VVLoose = cms.string('0.1686942'),
        VVTight = cms.string('0.9928449'),
        VVVLoose = cms.string('0.0630386')
    ),
    VSjetWP = cms.PSet(
        Loose = cms.string('0.7848675'),
        Medium = cms.string('0.8834768'),
        Tight = cms.string('0.9308689'),
        VLoose = cms.string('0.5983682'),
        VTight = cms.string('0.9573137'),
        VVLoose = cms.string('0.4249705'),
        VVTight = cms.string('0.9733927'),
        VVVLoose = cms.string('0.2599605')
    ),
    VSmuWP = cms.PSet(
        Loose = cms.string('0.2158633'),
        Medium = cms.string('0.5551894'),
        Tight = cms.string('0.8754835'),
        VLoose = cms.string('0.1058354')
    ),
    debug_level = cms.int32(0),
    disable_dxy_pca = cms.bool(True),
    electrons = cms.InputTag("slimmedElectrons"),
    graph_file = cms.vstring(
        'core:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_core.pb', 
        'inner:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_inner.pb', 
        'outer:RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v2p6_e6_outer.pb'
    ),
    mem_mapped = cms.bool(True),
    muons = cms.InputTag("slimmedMuons"),
    pfcands = cms.InputTag("packedPFCandidates"),
    rho = cms.InputTag("fixedGridRhoAll"),
    taus = cms.InputTag("slimmedTaus"),
    version = cms.uint32(2),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.egmPhotonIDs = cms.EDProducer("VersionedPhotonIdProducer",
    physicsObjectIDs = cms.VPSet(),
    physicsObjectSrc = cms.InputTag("slimmedPhotons")
)


process.egmPhotonIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
    isolationConeDefinitions = cms.VPSet(
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h+'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h0'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('gamma'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
            vertexIndex = cms.int32(0)
        )
    ),
    srcForIsolationCone = cms.InputTag("packedPFCandidates"),
    srcToIsolate = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
)


process.elPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(True),
        SCMatch_Veto = cms.bool(False),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoValueCharged03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.electronMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(0.5),
    maxDeltaR = cms.double(0.5),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons")
)


process.genMetExtractor = cms.EDProducer("GenMETExtractor",
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.genMetExtractorPuppi = cms.EDProducer("GenMETExtractor",
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.heepIDVarValueMaps = cms.EDProducer("ElectronHEEPIDValueMapProducer",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    candVetosAOD = cms.vstring(
        'ELES', 
        'NONE', 
        'NONELES'
    ),
    candVetosMiniAOD = cms.vstring(
        'ELES', 
        'NONE', 
        'NONELES'
    ),
    candsAOD = cms.VInputTag("packedCandsForTkIso", "lostTracksForTkIso", "lostTracksForTkIso:eleTracks"),
    candsMiniAOD = cms.VInputTag("packedPFCandidates", "lostTracks", "lostTracks:eleTracks"),
    dataFormat = cms.int32(2),
    ebRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    ebRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    eeRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    elesAOD = cms.InputTag("gedGsfElectrons"),
    elesMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
    makeTrkIso04 = cms.bool(True),
    trkIso04Config = cms.PSet(
        barrelCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.4),
            maxDZ = cms.double(0.1),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        ),
        endcapCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.4),
            maxDZ = cms.double(0.5),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        )
    ),
    trkIsoConfig = cms.PSet(
        barrelCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.1),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        ),
        endcapCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.5),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        )
    )
)


process.hpsPFTauChargedIsoPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByLoosePileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByMediumPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(1.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(True),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByRawPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByTightPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    customOuterCone = cms.double(-1.0),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(0.8),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauFootprintCorrection = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(True),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauNeutralIsoPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauNeutralIsoPtSumWeight = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauPUcorrPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(True),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(True),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauPhotonPtSumOutsideSignalCone = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.bool(True),
    deltaBetaPUTrackPtCutOverride_val = cms.double(0.5),
    footprintCorrections = cms.VPSet(
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 0')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )
    ),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(1.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(False),
    storeRawOccupancy = cms.bool(False),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(True),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.ic5CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1FastL2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1FastL2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1L2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1L2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1FastL2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1FastL2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1L2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1L2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.islandPhotonCore = cms.EDProducer("PhotonCoreProducer",
    conversionProducer = cms.InputTag(""),
    minSCEt = cms.double(8.0),
    photonCoreCollection = cms.string(''),
    pixelSeedProducer = cms.InputTag("electronMergedSeeds"),
    risolveConversionAmbiguity = cms.bool(True),
    scHybridBarrelProducer = cms.InputTag("correctedIslandBarrelSuperClusters"),
    scIslandEndcapProducer = cms.InputTag("correctedIslandEndcapSuperClusters")
)


process.kt4CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1FastL2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1FastL2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1L2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1L2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1FastL2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1FastL2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1L2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1L2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt6CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1FastL2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1FastL2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1L2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1L2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1FastL2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1FastL2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1L2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1L2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.lowPtElectronMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(0.5),
    maxDeltaR = cms.double(0.5),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("lowPtGsfElectrons")
)


process.metrawCalo = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('rawCalo'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.metrawCaloPuppi = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('rawCalo'),
    metSource = cms.InputTag("slimmedMETsPuppi","","@skipCurrentProcess")
)


process.muCaloMetCorr = cms.EDProducer("MuonMETcorrInputProducer",
    src = cms.InputTag("muons"),
    srcMuonCorrections = cms.InputTag("muonMETValueMapProducer","muCorrData")
)


process.muPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(
            '0.0001', 
            'Threshold(0.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(1.0)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring(
            '0.01', 
            'Threshold(0.5)'
        ),
        weight = cms.string('1')
    ))
)


process.muonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(0.5),
    maxDeltaR = cms.double(0.5),
    mcPdgId = cms.vint32(13),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("muons")
)


process.ootPhotonCore = cms.EDProducer("PhotonCoreProducer",
    conversionProducer = cms.InputTag("conversions"),
    minSCEt = cms.double(10.0),
    photonCoreCollection = cms.string(''),
    pixelSeedProducer = cms.InputTag("electronMergedSeeds"),
    risolveConversionAmbiguity = cms.bool(True),
    scHybridBarrelProducer = cms.InputTag("particleFlowSuperClusterOOTECAL","particleFlowSuperClusterOOTECALBarrel"),
    scIslandEndcapProducer = cms.InputTag("particleFlowSuperClusterOOTECAL","particleFlowSuperClusterOOTECALEndcapWithPreshower")
)


process.ootPhotonEcalPFClusterIsolationProducer = cms.EDProducer("PhotonEcalPFClusterIsolationProducer",
    candidateProducer = cms.InputTag("ootPhotonsTmp"),
    drMax = cms.double(0.3),
    drVetoBarrel = cms.double(0),
    drVetoEndcap = cms.double(0),
    energyBarrel = cms.double(0),
    energyEndcap = cms.double(0),
    etaStripBarrel = cms.double(0),
    etaStripEndcap = cms.double(0),
    pfClusterProducer = cms.InputTag("particleFlowClusterOOTECAL")
)


process.ootPhotonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(1.0),
    maxDeltaR = cms.double(0.2),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("ootPhotons")
)


process.ootPhotons = cms.EDProducer("GEDPhotonProducer",
    PFIsolationCalculatorSet = cms.PSet(
        ComponentName = cms.string('pfIsolationCalculator'),
        applyDzDxyVeto = cms.bool(True),
        applyMissHitPhVeto = cms.bool(False),
        applyPFPUVeto = cms.bool(True),
        applyVeto = cms.bool(True),
        checkClosestZVertex = cms.bool(True),
        coneDR = cms.double(0.3),
        deltaRVetoBarrel = cms.bool(True),
        deltaRVetoBarrelCharged = cms.double(0.02),
        deltaRVetoBarrelNeutrals = cms.double(-1.0),
        deltaRVetoBarrelPhotons = cms.double(-1.0),
        deltaRVetoEndcap = cms.bool(True),
        deltaRVetoEndcapCharged = cms.double(0.02),
        deltaRVetoEndcapNeutrals = cms.double(-1.0),
        deltaRVetoEndcapPhotons = cms.double(0.07),
        numberOfCrystalEndcapPhotons = cms.double(4.0),
        numberOfRings = cms.int32(1),
        particleType = cms.int32(1),
        rectangleDeltaEtaVetoBarrelCharged = cms.double(-1),
        rectangleDeltaEtaVetoBarrelNeutrals = cms.double(-1),
        rectangleDeltaEtaVetoBarrelPhotons = cms.double(0.015),
        rectangleDeltaEtaVetoEndcapCharged = cms.double(-1),
        rectangleDeltaEtaVetoEndcapNeutrals = cms.double(-1),
        rectangleDeltaEtaVetoEndcapPhotons = cms.double(-1),
        rectangleDeltaPhiVetoBarrelCharged = cms.double(-1),
        rectangleDeltaPhiVetoBarrelNeutrals = cms.double(-1),
        rectangleDeltaPhiVetoBarrelPhotons = cms.double(1.0),
        rectangleDeltaPhiVetoEndcapCharged = cms.double(-1),
        rectangleDeltaPhiVetoEndcapNeutrals = cms.double(-1),
        rectangleDeltaPhiVetoEndcapPhotons = cms.double(-1),
        rectangleVetoBarrel = cms.bool(True),
        rectangleVetoEndcap = cms.bool(False),
        ringSize = cms.double(0.3),
        useCrystalSize = cms.bool(True)
    ),
    RecHitFlagToBeExcludedEB = cms.vstring(
        'kFaultyHardware', 
        'kTowerRecovered', 
        'kDead'
    ),
    RecHitFlagToBeExcludedEE = cms.vstring(
        'kFaultyHardware', 
        'kNeighboursRecovered', 
        'kTowerRecovered', 
        'kDead', 
        'kWeird'
    ),
    RecHitSeverityToBeExcludedEB = cms.vstring(
        'kWeird', 
        'kBad', 
        'kTime'
    ),
    RecHitSeverityToBeExcludedEE = cms.vstring(
        'kWeird', 
        'kBad', 
        'kTime'
    ),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    candidateP4type = cms.string('fromEcalEnergy'),
    checkHcalStatus = cms.bool(True),
    ecalRecHitSumEtOffsetBarrel = cms.double(999999999),
    ecalRecHitSumEtOffsetEndcap = cms.double(999999999),
    ecalRecHitSumEtSlopeBarrel = cms.double(0.0),
    ecalRecHitSumEtSlopeEndcap = cms.double(0.0),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    energyRegressionWeightsDBLocation = cms.string('wgbrph'),
    energyRegressionWeightsFileLocation = cms.string('/afs/cern.ch/user/b/bendavid/cmspublic/regweights/gbrph.root'),
    hOverEConeSize = cms.double(0.15),
    hbheInstance = cms.string(''),
    hbheModule = cms.string('hbhereco'),
    hcalTowerSumEtOffsetBarrel = cms.double(999999999),
    hcalTowerSumEtOffsetEndcap = cms.double(999999999),
    hcalTowerSumEtSlopeBarrel = cms.double(0.0),
    hcalTowerSumEtSlopeEndcap = cms.double(0.0),
    hcalTowers = cms.InputTag("towerMaker"),
    highEt = cms.double(100.0),
    isolationSumsCalculatorSet = cms.PSet(
        ComponentName = cms.string('isolationSumsCalculator'),
        EcalRecHitEtaSliceA_Barrel = cms.double(2.5),
        EcalRecHitEtaSliceA_Endcap = cms.double(2.5),
        EcalRecHitEtaSliceB_Barrel = cms.double(2.5),
        EcalRecHitEtaSliceB_Endcap = cms.double(2.5),
        EcalRecHitInnerRadiusA_Barrel = cms.double(3.5),
        EcalRecHitInnerRadiusA_Endcap = cms.double(3.5),
        EcalRecHitInnerRadiusB_Barrel = cms.double(3.5),
        EcalRecHitInnerRadiusB_Endcap = cms.double(3.5),
        EcalRecHitOuterRadiusA_Barrel = cms.double(0.4),
        EcalRecHitOuterRadiusA_Endcap = cms.double(0.4),
        EcalRecHitOuterRadiusB_Barrel = cms.double(0.3),
        EcalRecHitOuterRadiusB_Endcap = cms.double(0.3),
        EcalRecHitThreshEA_Barrel = cms.double(0.095),
        EcalRecHitThreshEA_Endcap = cms.double(0.0),
        EcalRecHitThreshEB_Barrel = cms.double(0.095),
        EcalRecHitThreshEB_Endcap = cms.double(0.0),
        EcalRecHitThreshEtA_Barrel = cms.double(0.0),
        EcalRecHitThreshEtA_Endcap = cms.double(0.11),
        EcalRecHitThreshEtB_Barrel = cms.double(0.0),
        EcalRecHitThreshEtB_Endcap = cms.double(0.11),
        HcalDepth1TowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalDepth1TowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalDepth1TowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalDepth1TowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalDepth1TowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalDepth1TowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalDepth1TowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalDepth1TowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalDepth1TowerThreshEA_Barrel = cms.double(0.0),
        HcalDepth1TowerThreshEA_Endcap = cms.double(0.0),
        HcalDepth1TowerThreshEB_Barrel = cms.double(0.0),
        HcalDepth1TowerThreshEB_Endcap = cms.double(0.0),
        HcalDepth2TowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalDepth2TowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalDepth2TowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalDepth2TowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalDepth2TowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalDepth2TowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalDepth2TowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalDepth2TowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalDepth2TowerThreshEA_Barrel = cms.double(0.0),
        HcalDepth2TowerThreshEA_Endcap = cms.double(0.0),
        HcalDepth2TowerThreshEB_Barrel = cms.double(0.0),
        HcalDepth2TowerThreshEB_Endcap = cms.double(0.0),
        HcalRecHitCollection = cms.InputTag("towerMaker"),
        HcalTowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalTowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalTowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalTowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalTowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalTowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalTowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalTowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalTowerThreshEA_Barrel = cms.double(0.0),
        HcalTowerThreshEA_Endcap = cms.double(0.0),
        HcalTowerThreshEB_Barrel = cms.double(0.0),
        HcalTowerThreshEB_Endcap = cms.double(0.0),
        TrackConeInnerRadiusA_Barrel = cms.double(0.04),
        TrackConeInnerRadiusA_Endcap = cms.double(0.04),
        TrackConeInnerRadiusB_Barrel = cms.double(0.04),
        TrackConeInnerRadiusB_Endcap = cms.double(0.04),
        TrackConeOuterRadiusA_Barrel = cms.double(0.4),
        TrackConeOuterRadiusA_Endcap = cms.double(0.4),
        TrackConeOuterRadiusB_Barrel = cms.double(0.3),
        TrackConeOuterRadiusB_Endcap = cms.double(0.3),
        barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        beamSpotProducer = cms.InputTag("offlineBeamSpot"),
        endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        isolationtrackEtaSliceA_Barrel = cms.double(0.015),
        isolationtrackEtaSliceA_Endcap = cms.double(0.015),
        isolationtrackEtaSliceB_Barrel = cms.double(0.015),
        isolationtrackEtaSliceB_Endcap = cms.double(0.015),
        isolationtrackThresholdA_Barrel = cms.double(0.0),
        isolationtrackThresholdA_Endcap = cms.double(0.0),
        isolationtrackThresholdB_Barrel = cms.double(0.0),
        isolationtrackThresholdB_Endcap = cms.double(0.0),
        longImpactParameterA_Barrel = cms.double(0.2),
        longImpactParameterA_Endcap = cms.double(0.2),
        longImpactParameterB_Barrel = cms.double(0.2),
        longImpactParameterB_Endcap = cms.double(0.2),
        moduleEtaBoundary = cms.vdouble(
            0.0, 0.02, 0.43, 0.46, 0.78, 
            0.81, 1.13, 1.15, 1.45, 1.58
        ),
        modulePhiBoundary = cms.double(0.0087),
        trackProducer = cms.InputTag("generalTracks"),
        transImpactParameterA_Barrel = cms.double(0.1),
        transImpactParameterA_Endcap = cms.double(0.1),
        transImpactParameterB_Barrel = cms.double(0.1),
        transImpactParameterB_Endcap = cms.double(0.1),
        useNumCrystals = cms.bool(True),
        vetoClustered = cms.bool(False)
    ),
    maxHoverEBarrel = cms.double(0.5),
    maxHoverEEndcap = cms.double(0.5),
    minR9Barrel = cms.double(0.94),
    minR9Endcap = cms.double(0.95),
    minSCEtBarrel = cms.double(10.0),
    minSCEtEndcap = cms.double(10.0),
    mipVariableSet = cms.PSet(
        ComponentName = cms.string('mipVariable'),
        HaloDiscThreshold = cms.double(70.0),
        ResidualWidth = cms.double(0.23),
        XRangeFit = cms.double(180.0),
        YRangeFit = cms.double(7.0),
        barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE")
    ),
    nTrackHollowConeBarrel = cms.double(999999999),
    nTrackHollowConeEndcap = cms.double(999999999),
    nTrackSolidConeBarrel = cms.double(999999999),
    nTrackSolidConeEndcap = cms.double(999999999),
    outputPhotonCollection = cms.string(''),
    pfCandidates = cms.InputTag("particleFlowTmp"),
    pfECALClusIsolation = cms.InputTag("ootPhotonEcalPFClusterIsolationProducer"),
    pfEgammaCandidates = cms.InputTag(""),
    pfHCALClusIsolation = cms.InputTag("ootPhotonHcalPFClusterIsolationProducer"),
    pfIsolCfg = cms.PSet(
        chargedHadronIso = cms.InputTag(""),
        chargedHadronPFPVIso = cms.InputTag(""),
        chargedHadronWorstVtxGeomVetoIso = cms.InputTag(""),
        chargedHadronWorstVtxIso = cms.InputTag(""),
        neutralHadronIso = cms.InputTag(""),
        photonIso = cms.InputTag("")
    ),
    photonEcalEnergyCorrFunction = cms.string('EcalClusterEnergyCorrectionObjectSpecific'),
    photonProducer = cms.InputTag("ootPhotonsTmp"),
    posCalcParameters = cms.PSet(
        LogWeighted = cms.bool(True),
        T0_barl = cms.double(7.4),
        T0_endc = cms.double(6.3),
        T0_endcPresh = cms.double(3.6),
        W0 = cms.double(4.2),
        X0 = cms.double(0.89)
    ),
    posCalc_logweight = cms.bool(True),
    posCalc_t0_barl = cms.double(7.7),
    posCalc_t0_endc = cms.double(6.3),
    posCalc_t0_endcPresh = cms.double(3.6),
    posCalc_w0 = cms.double(4.2),
    posCalc_x0 = cms.double(0.89),
    preshowerHits = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    primaryVertexProducer = cms.InputTag("offlinePrimaryVerticesWithBS"),
    reconstructionStep = cms.string('ootfinal'),
    regressionWeightsFromDB = cms.bool(True),
    runMIPTagger = cms.bool(True),
    sigmaIetaIetaCutBarrel = cms.double(999999999),
    sigmaIetaIetaCutEndcap = cms.double(999999999),
    superClusterCrackEnergyCorrFunction = cms.string('EcalClusterCrackCorrection'),
    superClusterEnergyCorrFunction = cms.string('EcalClusterEnergyCorrection'),
    superClusterEnergyErrorFunction = cms.string('EcalClusterEnergyUncertainty'),
    trackPtSumHollowConeBarrel = cms.double(999999999),
    trackPtSumHollowConeEndcap = cms.double(999999999),
    trackPtSumSolidConeBarrel = cms.double(999999999),
    trackPtSumSolidConeEndcap = cms.double(999999999),
    usePrimaryVertex = cms.bool(True),
    useRegression = cms.bool(True),
    valueMapPhotons = cms.string('')
)


process.ootPhotonsTmp = cms.EDProducer("GEDPhotonProducer",
    PFIsolationCalculatorSet = cms.PSet(
        ComponentName = cms.string('pfIsolationCalculator'),
        applyDzDxyVeto = cms.bool(True),
        applyMissHitPhVeto = cms.bool(False),
        applyPFPUVeto = cms.bool(True),
        applyVeto = cms.bool(True),
        checkClosestZVertex = cms.bool(True),
        coneDR = cms.double(0.3),
        deltaRVetoBarrel = cms.bool(True),
        deltaRVetoBarrelCharged = cms.double(0.02),
        deltaRVetoBarrelNeutrals = cms.double(-1.0),
        deltaRVetoBarrelPhotons = cms.double(-1.0),
        deltaRVetoEndcap = cms.bool(True),
        deltaRVetoEndcapCharged = cms.double(0.02),
        deltaRVetoEndcapNeutrals = cms.double(-1.0),
        deltaRVetoEndcapPhotons = cms.double(0.07),
        numberOfCrystalEndcapPhotons = cms.double(4.0),
        numberOfRings = cms.int32(1),
        particleType = cms.int32(1),
        rectangleDeltaEtaVetoBarrelCharged = cms.double(-1),
        rectangleDeltaEtaVetoBarrelNeutrals = cms.double(-1),
        rectangleDeltaEtaVetoBarrelPhotons = cms.double(0.015),
        rectangleDeltaEtaVetoEndcapCharged = cms.double(-1),
        rectangleDeltaEtaVetoEndcapNeutrals = cms.double(-1),
        rectangleDeltaEtaVetoEndcapPhotons = cms.double(-1),
        rectangleDeltaPhiVetoBarrelCharged = cms.double(-1),
        rectangleDeltaPhiVetoBarrelNeutrals = cms.double(-1),
        rectangleDeltaPhiVetoBarrelPhotons = cms.double(1.0),
        rectangleDeltaPhiVetoEndcapCharged = cms.double(-1),
        rectangleDeltaPhiVetoEndcapNeutrals = cms.double(-1),
        rectangleDeltaPhiVetoEndcapPhotons = cms.double(-1),
        rectangleVetoBarrel = cms.bool(True),
        rectangleVetoEndcap = cms.bool(False),
        ringSize = cms.double(0.3),
        useCrystalSize = cms.bool(True)
    ),
    RecHitFlagToBeExcludedEB = cms.vstring(
        'kFaultyHardware', 
        'kTowerRecovered', 
        'kDead'
    ),
    RecHitFlagToBeExcludedEE = cms.vstring(
        'kFaultyHardware', 
        'kNeighboursRecovered', 
        'kTowerRecovered', 
        'kDead', 
        'kWeird'
    ),
    RecHitSeverityToBeExcludedEB = cms.vstring(
        'kWeird', 
        'kBad', 
        'kTime'
    ),
    RecHitSeverityToBeExcludedEE = cms.vstring(
        'kWeird', 
        'kBad', 
        'kTime'
    ),
    barrelEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    candidateP4type = cms.string('fromEcalEnergy'),
    checkHcalStatus = cms.bool(True),
    ecalRecHitSumEtOffsetBarrel = cms.double(999999999),
    ecalRecHitSumEtOffsetEndcap = cms.double(999999999),
    ecalRecHitSumEtSlopeBarrel = cms.double(0.0),
    ecalRecHitSumEtSlopeEndcap = cms.double(0.0),
    endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    energyRegressionWeightsDBLocation = cms.string('wgbrph'),
    energyRegressionWeightsFileLocation = cms.string('/afs/cern.ch/user/b/bendavid/cmspublic/regweights/gbrph.root'),
    hOverEConeSize = cms.double(0.15),
    hbheInstance = cms.string(''),
    hbheModule = cms.string('hbhereco'),
    hcalTowerSumEtOffsetBarrel = cms.double(999999999),
    hcalTowerSumEtOffsetEndcap = cms.double(999999999),
    hcalTowerSumEtSlopeBarrel = cms.double(0.0),
    hcalTowerSumEtSlopeEndcap = cms.double(0.0),
    hcalTowers = cms.InputTag("towerMaker"),
    highEt = cms.double(100.0),
    isolationSumsCalculatorSet = cms.PSet(
        ComponentName = cms.string('isolationSumsCalculator'),
        EcalRecHitEtaSliceA_Barrel = cms.double(2.5),
        EcalRecHitEtaSliceA_Endcap = cms.double(2.5),
        EcalRecHitEtaSliceB_Barrel = cms.double(2.5),
        EcalRecHitEtaSliceB_Endcap = cms.double(2.5),
        EcalRecHitInnerRadiusA_Barrel = cms.double(3.5),
        EcalRecHitInnerRadiusA_Endcap = cms.double(3.5),
        EcalRecHitInnerRadiusB_Barrel = cms.double(3.5),
        EcalRecHitInnerRadiusB_Endcap = cms.double(3.5),
        EcalRecHitOuterRadiusA_Barrel = cms.double(0.4),
        EcalRecHitOuterRadiusA_Endcap = cms.double(0.4),
        EcalRecHitOuterRadiusB_Barrel = cms.double(0.3),
        EcalRecHitOuterRadiusB_Endcap = cms.double(0.3),
        EcalRecHitThreshEA_Barrel = cms.double(0.095),
        EcalRecHitThreshEA_Endcap = cms.double(0.0),
        EcalRecHitThreshEB_Barrel = cms.double(0.095),
        EcalRecHitThreshEB_Endcap = cms.double(0.0),
        EcalRecHitThreshEtA_Barrel = cms.double(0.0),
        EcalRecHitThreshEtA_Endcap = cms.double(0.11),
        EcalRecHitThreshEtB_Barrel = cms.double(0.0),
        EcalRecHitThreshEtB_Endcap = cms.double(0.11),
        HcalDepth1TowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalDepth1TowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalDepth1TowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalDepth1TowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalDepth1TowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalDepth1TowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalDepth1TowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalDepth1TowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalDepth1TowerThreshEA_Barrel = cms.double(0.0),
        HcalDepth1TowerThreshEA_Endcap = cms.double(0.0),
        HcalDepth1TowerThreshEB_Barrel = cms.double(0.0),
        HcalDepth1TowerThreshEB_Endcap = cms.double(0.0),
        HcalDepth2TowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalDepth2TowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalDepth2TowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalDepth2TowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalDepth2TowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalDepth2TowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalDepth2TowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalDepth2TowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalDepth2TowerThreshEA_Barrel = cms.double(0.0),
        HcalDepth2TowerThreshEA_Endcap = cms.double(0.0),
        HcalDepth2TowerThreshEB_Barrel = cms.double(0.0),
        HcalDepth2TowerThreshEB_Endcap = cms.double(0.0),
        HcalRecHitCollection = cms.InputTag("towerMaker"),
        HcalTowerInnerRadiusA_Barrel = cms.double(0.15),
        HcalTowerInnerRadiusA_Endcap = cms.double(0.15),
        HcalTowerInnerRadiusB_Barrel = cms.double(0.15),
        HcalTowerInnerRadiusB_Endcap = cms.double(0.15),
        HcalTowerOuterRadiusA_Barrel = cms.double(0.4),
        HcalTowerOuterRadiusA_Endcap = cms.double(0.4),
        HcalTowerOuterRadiusB_Barrel = cms.double(0.3),
        HcalTowerOuterRadiusB_Endcap = cms.double(0.3),
        HcalTowerThreshEA_Barrel = cms.double(0.0),
        HcalTowerThreshEA_Endcap = cms.double(0.0),
        HcalTowerThreshEB_Barrel = cms.double(0.0),
        HcalTowerThreshEB_Endcap = cms.double(0.0),
        TrackConeInnerRadiusA_Barrel = cms.double(0.04),
        TrackConeInnerRadiusA_Endcap = cms.double(0.04),
        TrackConeInnerRadiusB_Barrel = cms.double(0.04),
        TrackConeInnerRadiusB_Endcap = cms.double(0.04),
        TrackConeOuterRadiusA_Barrel = cms.double(0.4),
        TrackConeOuterRadiusA_Endcap = cms.double(0.4),
        TrackConeOuterRadiusB_Barrel = cms.double(0.3),
        TrackConeOuterRadiusB_Endcap = cms.double(0.3),
        barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        beamSpotProducer = cms.InputTag("offlineBeamSpot"),
        endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        isolationtrackEtaSliceA_Barrel = cms.double(0.015),
        isolationtrackEtaSliceA_Endcap = cms.double(0.015),
        isolationtrackEtaSliceB_Barrel = cms.double(0.015),
        isolationtrackEtaSliceB_Endcap = cms.double(0.015),
        isolationtrackThresholdA_Barrel = cms.double(0.0),
        isolationtrackThresholdA_Endcap = cms.double(0.0),
        isolationtrackThresholdB_Barrel = cms.double(0.0),
        isolationtrackThresholdB_Endcap = cms.double(0.0),
        longImpactParameterA_Barrel = cms.double(0.2),
        longImpactParameterA_Endcap = cms.double(0.2),
        longImpactParameterB_Barrel = cms.double(0.2),
        longImpactParameterB_Endcap = cms.double(0.2),
        moduleEtaBoundary = cms.vdouble(
            0.0, 0.02, 0.43, 0.46, 0.78, 
            0.81, 1.13, 1.15, 1.45, 1.58
        ),
        modulePhiBoundary = cms.double(0.0087),
        trackProducer = cms.InputTag("generalTracks"),
        transImpactParameterA_Barrel = cms.double(0.1),
        transImpactParameterA_Endcap = cms.double(0.1),
        transImpactParameterB_Barrel = cms.double(0.1),
        transImpactParameterB_Endcap = cms.double(0.1),
        useNumCrystals = cms.bool(True),
        vetoClustered = cms.bool(False)
    ),
    maxHoverEBarrel = cms.double(0.5),
    maxHoverEEndcap = cms.double(0.5),
    minR9Barrel = cms.double(0.94),
    minR9Endcap = cms.double(0.95),
    minSCEtBarrel = cms.double(10.0),
    minSCEtEndcap = cms.double(10.0),
    mipVariableSet = cms.PSet(
        ComponentName = cms.string('mipVariable'),
        HaloDiscThreshold = cms.double(70.0),
        ResidualWidth = cms.double(0.23),
        XRangeFit = cms.double(180.0),
        YRangeFit = cms.double(7.0),
        barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE")
    ),
    nTrackHollowConeBarrel = cms.double(999999999),
    nTrackHollowConeEndcap = cms.double(999999999),
    nTrackSolidConeBarrel = cms.double(999999999),
    nTrackSolidConeEndcap = cms.double(999999999),
    outputPhotonCollection = cms.string(''),
    pfCandidates = cms.InputTag("particleFlowTmp"),
    pfEgammaCandidates = cms.InputTag(""),
    photonEcalEnergyCorrFunction = cms.string('EcalClusterEnergyCorrectionObjectSpecific'),
    photonProducer = cms.InputTag("ootPhotonCore"),
    posCalcParameters = cms.PSet(
        LogWeighted = cms.bool(True),
        T0_barl = cms.double(7.4),
        T0_endc = cms.double(6.3),
        T0_endcPresh = cms.double(3.6),
        W0 = cms.double(4.2),
        X0 = cms.double(0.89)
    ),
    posCalc_logweight = cms.bool(True),
    posCalc_t0_barl = cms.double(7.7),
    posCalc_t0_endc = cms.double(6.3),
    posCalc_t0_endcPresh = cms.double(3.6),
    posCalc_w0 = cms.double(4.2),
    posCalc_x0 = cms.double(0.89),
    preshowerHits = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
    primaryVertexProducer = cms.InputTag("offlinePrimaryVerticesWithBS"),
    reconstructionStep = cms.string('oot'),
    regressionWeightsFromDB = cms.bool(True),
    runMIPTagger = cms.bool(True),
    sigmaIetaIetaCutBarrel = cms.double(999999999),
    sigmaIetaIetaCutEndcap = cms.double(999999999),
    superClusterCrackEnergyCorrFunction = cms.string('EcalClusterCrackCorrection'),
    superClusterEnergyCorrFunction = cms.string('EcalClusterEnergyCorrection'),
    superClusterEnergyErrorFunction = cms.string('EcalClusterEnergyUncertainty'),
    trackPtSumHollowConeBarrel = cms.double(999999999),
    trackPtSumHollowConeEndcap = cms.double(999999999),
    trackPtSumSolidConeBarrel = cms.double(999999999),
    trackPtSumSolidConeEndcap = cms.double(999999999),
    usePrimaryVertex = cms.bool(True),
    useRegression = cms.bool(True),
    valueMapPhotons = cms.string('')
)


process.particleFlowClusterECAL = cms.EDProducer("CorrectedECALPFClusterProducer",
    energyCorrector = cms.PSet(
        algoName = cms.string('PFClusterEMEnergyCorrector'),
        applyCrackCorrections = cms.bool(False),
        applyMVACorrections = cms.bool(True),
        autoDetectBunchSpacing = cms.bool(True),
        bunchSpacing = cms.int32(25),
        ebSrFlagLabel = cms.InputTag("ecalDigis"),
        eeSrFlagLabel = cms.InputTag("ecalDigis"),
        maxPtForMVAEvaluation = cms.double(90.0),
        recHitsEBLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        recHitsEELabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        setEnergyUncertainty = cms.bool(False),
        srfAwareCorrection = cms.bool(False),
        verticesLabel = cms.InputTag("offlinePrimaryVertices")
    ),
    inputECAL = cms.InputTag("particleFlowClusterECALUncorrected"),
    inputPS = cms.InputTag("particleFlowClusterPS"),
    minimumPSEnergy = cms.double(0)
)


process.particleFlowClusterECALUncorrected = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                gatheringThreshold = cms.double(0.08),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                gatheringThreshold = cms.double(0.3),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(True)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(9),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        positionCalcForConvergence = cms.PSet(
            T0_EB = cms.double(7.4),
            T0_EE = cms.double(3.1),
            T0_ES = cms.double(1.2),
            W0 = cms.double(4.2),
            X0 = cms.double(0.89),
            algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
            minAllowedNormalization = cms.double(0.0),
            minFractionInCalc = cms.double(0.0)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                recHitEnergyNorm = cms.double(0.08)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                recHitEnergyNorm = cms.double(0.3)
            )
        ),
        showerSigma = cms.double(1.5),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(
        T0_EB = cms.double(7.4),
        T0_EE = cms.double(3.1),
        T0_ES = cms.double(1.2),
        W0 = cms.double(4.2),
        X0 = cms.double(0.89),
        algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
        minAllowedNormalization = cms.double(0.0),
        minFractionInCalc = cms.double(0.0)
    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("particleFlowRecHitECAL"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(8),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                seedingThreshold = cms.double(0.6),
                seedingThresholdPt = cms.double(0.15)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                seedingThreshold = cms.double(0.23),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.particleFlowClusterOOTECAL = cms.EDProducer("CorrectedECALPFClusterProducer",
    energyCorrector = cms.PSet(
        algoName = cms.string('PFClusterEMEnergyCorrector'),
        applyCrackCorrections = cms.bool(False),
        applyMVACorrections = cms.bool(True),
        autoDetectBunchSpacing = cms.bool(True),
        bunchSpacing = cms.int32(25),
        ebSrFlagLabel = cms.InputTag("ecalDigis"),
        eeSrFlagLabel = cms.InputTag("ecalDigis"),
        maxPtForMVAEvaluation = cms.double(90.0),
        recHitsEBLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        recHitsEELabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        setEnergyUncertainty = cms.bool(False),
        srfAwareCorrection = cms.bool(False),
        verticesLabel = cms.InputTag("offlinePrimaryVertices")
    ),
    inputECAL = cms.InputTag("particleFlowClusterOOTECALUncorrected"),
    inputPS = cms.InputTag("particleFlowClusterPS"),
    minimumPSEnergy = cms.double(0)
)


process.particleFlowClusterOOTECALUncorrected = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                gatheringThreshold = cms.double(0.08),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                gatheringThreshold = cms.double(0.3),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(True)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        allCellsPositionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(0.08),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(9),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(0.428192),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0510871),
                noiseTerm = cms.double(1.10889),
                noiseTermLowE = cms.double(1.31883),
                threshHighE = cms.double(5.0),
                threshLowE = cms.double(0.5)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(0.0),
                constantTermLowE = cms.double(0.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(5.72489999999),
                noiseTermLowE = cms.double(6.92683000001),
                threshHighE = cms.double(10.0),
                threshLowE = cms.double(1.0)
            )
        ),
        positionCalcForConvergence = cms.PSet(
            T0_EB = cms.double(7.4),
            T0_EE = cms.double(3.1),
            T0_ES = cms.double(1.2),
            W0 = cms.double(4.2),
            X0 = cms.double(0.89),
            algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
            minAllowedNormalization = cms.double(0.0),
            minFractionInCalc = cms.double(0.0)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                recHitEnergyNorm = cms.double(0.08)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                recHitEnergyNorm = cms.double(0.3)
            )
        ),
        showerSigma = cms.double(1.5),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(
        T0_EB = cms.double(7.4),
        T0_EE = cms.double(3.1),
        T0_ES = cms.double(1.2),
        W0 = cms.double(4.2),
        X0 = cms.double(0.89),
        algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
        minAllowedNormalization = cms.double(0.0),
        minFractionInCalc = cms.double(0.0)
    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("particleFlowRecHitOOTECAL"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(8),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('ECAL_ENDCAP'),
                seedingThreshold = cms.double(0.6),
                seedingThresholdPt = cms.double(0.15)
            ), 
            cms.PSet(
                detector = cms.string('ECAL_BARREL'),
                seedingThreshold = cms.double(0.23),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.particleFlowClusterPS = cms.EDProducer("PFClusterProducer",
    energyCorrector = cms.PSet(

    ),
    initialClusteringStep = cms.PSet(
        algoName = cms.string('Basic2DGenericTopoClusterizer'),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                gatheringThreshold = cms.double(6e-05),
                gatheringThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                gatheringThreshold = cms.double(6e-05),
                gatheringThresholdPt = cms.double(0.0)
            )
        ),
        useCornerCells = cms.bool(False)
    ),
    pfClusterBuilder = cms.PSet(
        algoName = cms.string('Basic2DGenericPFlowClusterizer'),
        excludeOtherSeeds = cms.bool(True),
        maxIterations = cms.uint32(50),
        minFracTot = cms.double(1e-20),
        minFractionToKeep = cms.double(1e-07),
        positionCalc = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
            logWeightDenominator = cms.double(6e-05),
            minAllowedNormalization = cms.double(1e-09),
            minFractionInCalc = cms.double(1e-09),
            posCalcNCrystals = cms.int32(-1)
        ),
        recHitEnergyNorms = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                recHitEnergyNorm = cms.double(6e-05)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                recHitEnergyNorm = cms.double(6e-05)
            )
        ),
        showerSigma = cms.double(0.3),
        stoppingTolerance = cms.double(1e-08)
    ),
    positionReCalc = cms.PSet(

    ),
    recHitCleaners = cms.VPSet(),
    recHitsSource = cms.InputTag("particleFlowRecHitPS"),
    seedFinder = cms.PSet(
        algoName = cms.string('LocalMaximumSeedFinder'),
        nNeighbours = cms.int32(4),
        thresholdsByDetector = cms.VPSet(
            cms.PSet(
                detector = cms.string('PS1'),
                seedingThreshold = cms.double(0.00012),
                seedingThresholdPt = cms.double(0.0)
            ), 
            cms.PSet(
                detector = cms.string('PS2'),
                seedingThreshold = cms.double(0.00012),
                seedingThresholdPt = cms.double(0.0)
            )
        )
    )
)


process.particleFlowDisplacedVertex = cms.EDProducer("PFDisplacedVertexProducer",
    avfParameters = cms.PSet(
        Tini = cms.double(256.0),
        ratio = cms.double(0.25),
        sigmacut = cms.double(6.0)
    ),
    debug = cms.untracked.bool(False),
    longSize = cms.double(5),
    mainVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    minAdaptWeight = cms.double(0.5),
    offlineBeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    primaryVertexCut = cms.double(1.8),
    switchOff2TrackVertex = cms.untracked.bool(True),
    tecCut = cms.double(220),
    tobCut = cms.double(100),
    tracksSelectorParameters = cms.PSet(
        bSelectTracks = cms.bool(True),
        dxy_min = cms.double(0.2),
        nChi2_max = cms.double(5.0),
        nChi2_min = cms.double(0.5),
        nHits_min = cms.int32(6),
        nOuterHits_max = cms.int32(9),
        pt_min = cms.double(0.2),
        quality = cms.string('HighPurity')
    ),
    transvSize = cms.double(1.0),
    verbose = cms.untracked.bool(False),
    vertexCandidatesLabel = cms.InputTag("particleFlowDisplacedVertexCandidate"),
    vertexIdentifierParameters = cms.PSet(
        angles = cms.vdouble(15, 15),
        bIdentifyVertices = cms.bool(True),
        logPrimSec_min = cms.double(0.0),
        looper_eta_max = cms.double(0.1),
        masses = cms.vdouble(
            0.05, 0.485, 0.515, 0.48, 0.52, 
            1.107, 1.125, 0.2
        ),
        pt_kink_min = cms.double(3.0),
        pt_min = cms.double(0.5)
    )
)


process.particleFlowPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
    src = cms.InputTag("particleFlow")
)


process.particleFlowRecHitECAL = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        barrel = cms.PSet(

        ),
        endcap = cms.PSet(

        ),
        name = cms.string('PFRecHitECALNavigator')
    ),
    producers = cms.VPSet(
        cms.PSet(
            name = cms.string('PFEBRecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    applySelectionsToAllCrystals = cms.bool(True),
                    name = cms.string('PFRecHitQTestDBThreshold')
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(True),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag(""),
            src = cms.InputTag("ecalRecHit","EcalRecHitsEB")
        ), 
        cms.PSet(
            name = cms.string('PFEERecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    applySelectionsToAllCrystals = cms.bool(True),
                    name = cms.string('PFRecHitQTestDBThreshold')
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(True),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag(""),
            src = cms.InputTag("ecalRecHit","EcalRecHitsEE")
        )
    )
)


process.particleFlowRecHitOOTECAL = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        barrel = cms.PSet(

        ),
        endcap = cms.PSet(

        ),
        name = cms.string('PFRecHitECALNavigator')
    ),
    producers = cms.VPSet(
        cms.PSet(
            name = cms.string('PFEBRecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    applySelectionsToAllCrystals = cms.bool(True),
                    name = cms.string('PFRecHitQTestDBThreshold')
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(False),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag(""),
            src = cms.InputTag("ecalRecHit","EcalRecHitsEB")
        ), 
        cms.PSet(
            name = cms.string('PFEERecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    applySelectionsToAllCrystals = cms.bool(True),
                    name = cms.string('PFRecHitQTestDBThreshold')
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(2.0),
                    name = cms.string('PFRecHitQTestECAL'),
                    skipTTRecoveredHits = cms.bool(True),
                    timingCleaning = cms.bool(False),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            srFlags = cms.InputTag(""),
            src = cms.InputTag("ecalRecHit","EcalRecHitsEE")
        )
    )
)


process.particleFlowRecHitPS = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
        name = cms.string('PFRecHitPreshowerNavigator')
    ),
    producers = cms.VPSet(cms.PSet(
        name = cms.string('PFPSRecHitCreator'),
        qualityTests = cms.VPSet(
            cms.PSet(
                name = cms.string('PFRecHitQTestThreshold'),
                threshold = cms.double(0.0)
            ), 
            cms.PSet(
                cleaningThreshold = cms.double(0.0),
                name = cms.string('PFRecHitQTestES'),
                topologicalCleaning = cms.bool(True)
            )
        ),
        src = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES")
    ))
)


process.particleFlowSuperClusterECAL = cms.EDProducer("PFECALSuperClusterProducer",
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ClusteringType = cms.string('Mustache'),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    EnergyWeight = cms.string('Raw'),
    PFBasicClusterCollectionBarrel = cms.string('particleFlowBasicClusterECALBarrel'),
    PFBasicClusterCollectionEndcap = cms.string('particleFlowBasicClusterECALEndcap'),
    PFBasicClusterCollectionPreshower = cms.string('particleFlowBasicClusterECALPreshower'),
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    PFSuperClusterCollectionBarrel = cms.string('particleFlowSuperClusterECALBarrel'),
    PFSuperClusterCollectionEndcap = cms.string('particleFlowSuperClusterECALEndcap'),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string('particleFlowSuperClusterECALEndcapWithPreshower'),
    applyCrackCorrections = cms.bool(False),
    doSatelliteClusterMerge = cms.bool(False),
    dropUnseedable = cms.bool(False),
    etawidth_SuperClusterBarrel = cms.double(0.04),
    etawidth_SuperClusterEndcap = cms.double(0.04),
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    phiwidth_SuperClusterEndcap = cms.double(0.6),
    regressionConfig = cms.PSet(
        applySigmaIetaIphiBug = cms.bool(False),
        ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        ecalRecHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
        regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
        uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
        uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
        vertexCollection = cms.InputTag("offlinePrimaryVertices")
    ),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    seedThresholdIsET = cms.bool(True),
    thresh_PFClusterBarrel = cms.double(0.0),
    thresh_PFClusterES = cms.double(0.0),
    thresh_PFClusterEndcap = cms.double(0.0),
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_SCEt = cms.double(4.0),
    useDynamicDPhiWindow = cms.bool(True),
    useRegression = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.particleFlowSuperClusterECALBox = cms.EDProducer("PFECALSuperClusterProducer",
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ClusteringType = cms.string('Box'),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    EnergyWeight = cms.string('Raw'),
    PFBasicClusterCollectionBarrel = cms.string('particleFlowBasicClusterECALBarrel'),
    PFBasicClusterCollectionEndcap = cms.string('particleFlowBasicClusterECALEndcap'),
    PFBasicClusterCollectionPreshower = cms.string('particleFlowBasicClusterECALPreshower'),
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    PFSuperClusterCollectionBarrel = cms.string('particleFlowSuperClusterECALBarrel'),
    PFSuperClusterCollectionEndcap = cms.string('particleFlowSuperClusterECALEndcap'),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string('particleFlowSuperClusterECALEndcapWithPreshower'),
    applyCrackCorrections = cms.bool(False),
    doSatelliteClusterMerge = cms.bool(False),
    dropUnseedable = cms.bool(False),
    etawidth_SuperClusterBarrel = cms.double(0.04),
    etawidth_SuperClusterEndcap = cms.double(0.04),
    phiwidth_SuperClusterBarrel = cms.double(0.28),
    phiwidth_SuperClusterEndcap = cms.double(0.28),
    regressionConfig = cms.PSet(
        applySigmaIetaIphiBug = cms.bool(False),
        ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        ecalRecHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
        regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
        uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
        uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
        vertexCollection = cms.InputTag("offlinePrimaryVertices")
    ),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    seedThresholdIsET = cms.bool(True),
    thresh_PFClusterBarrel = cms.double(0.5),
    thresh_PFClusterEndcap = cms.double(0.5),
    thresh_PFClusterSeedBarrel = cms.double(3.0),
    thresh_PFClusterSeedEndcap = cms.double(5.0),
    thresh_SCEt = cms.double(4.0),
    useDynamicDPhiWindow = cms.bool(False),
    useRegression = cms.bool(False),
    use_preshower = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.particleFlowSuperClusterECALMustache = cms.EDProducer("PFECALSuperClusterProducer",
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ClusteringType = cms.string('Mustache'),
    ESAssociation = cms.InputTag("particleFlowClusterECAL"),
    EnergyWeight = cms.string('Raw'),
    PFBasicClusterCollectionBarrel = cms.string('particleFlowBasicClusterECALBarrel'),
    PFBasicClusterCollectionEndcap = cms.string('particleFlowBasicClusterECALEndcap'),
    PFBasicClusterCollectionPreshower = cms.string('particleFlowBasicClusterECALPreshower'),
    PFClusters = cms.InputTag("particleFlowClusterECAL"),
    PFSuperClusterCollectionBarrel = cms.string('particleFlowSuperClusterECALBarrel'),
    PFSuperClusterCollectionEndcap = cms.string('particleFlowSuperClusterECALEndcap'),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string('particleFlowSuperClusterECALEndcapWithPreshower'),
    applyCrackCorrections = cms.bool(False),
    doSatelliteClusterMerge = cms.bool(False),
    dropUnseedable = cms.bool(False),
    etawidth_SuperClusterBarrel = cms.double(0.04),
    etawidth_SuperClusterEndcap = cms.double(0.04),
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    phiwidth_SuperClusterEndcap = cms.double(0.6),
    regressionConfig = cms.PSet(
        applySigmaIetaIphiBug = cms.bool(False),
        ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        ecalRecHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
        regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
        uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
        uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
        vertexCollection = cms.InputTag("offlinePrimaryVertices")
    ),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    seedThresholdIsET = cms.bool(True),
    thresh_PFClusterBarrel = cms.double(0.0),
    thresh_PFClusterES = cms.double(0.0),
    thresh_PFClusterEndcap = cms.double(0.0),
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_SCEt = cms.double(4.0),
    useDynamicDPhiWindow = cms.bool(True),
    useRegression = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.particleFlowSuperClusterOOTECAL = cms.EDProducer("PFECALSuperClusterProducer",
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ClusteringType = cms.string('Mustache'),
    ESAssociation = cms.InputTag("particleFlowClusterOOTECAL"),
    EnergyWeight = cms.string('Raw'),
    PFBasicClusterCollectionBarrel = cms.string('particleFlowBasicClusterOOTECALBarrel'),
    PFBasicClusterCollectionEndcap = cms.string('particleFlowBasicClusterOOTECALEndcap'),
    PFBasicClusterCollectionPreshower = cms.string('particleFlowBasicClusterOOTECALPreshower'),
    PFClusters = cms.InputTag("particleFlowClusterOOTECAL"),
    PFSuperClusterCollectionBarrel = cms.string('particleFlowSuperClusterOOTECALBarrel'),
    PFSuperClusterCollectionEndcap = cms.string('particleFlowSuperClusterOOTECALEndcap'),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string('particleFlowSuperClusterOOTECALEndcapWithPreshower'),
    applyCrackCorrections = cms.bool(False),
    barrelRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    doSatelliteClusterMerge = cms.bool(False),
    dropUnseedable = cms.bool(False),
    endcapRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    etawidth_SuperClusterBarrel = cms.double(0.04),
    etawidth_SuperClusterEndcap = cms.double(0.04),
    isOOTCollection = cms.bool(True),
    phiwidth_SuperClusterBarrel = cms.double(0.6),
    phiwidth_SuperClusterEndcap = cms.double(0.6),
    regressionConfig = cms.PSet(
        applySigmaIetaIphiBug = cms.bool(False),
        ecalRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        ecalRecHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        regressionKeyEB = cms.string('pfscecal_EBCorrection_offline_v2'),
        regressionKeyEE = cms.string('pfscecal_EECorrection_offline_v2'),
        uncertaintyKeyEB = cms.string('pfscecal_EBUncertainty_offline_v2'),
        uncertaintyKeyEE = cms.string('pfscecal_EEUncertainty_offline_v2'),
        vertexCollection = cms.InputTag("offlinePrimaryVertices")
    ),
    satelliteClusterSeedThreshold = cms.double(50.0),
    satelliteMajorityFraction = cms.double(0.5),
    seedThresholdIsET = cms.bool(True),
    thresh_PFClusterBarrel = cms.double(0.0),
    thresh_PFClusterES = cms.double(0.0),
    thresh_PFClusterEndcap = cms.double(0.0),
    thresh_PFClusterSeedBarrel = cms.double(1.0),
    thresh_PFClusterSeedEndcap = cms.double(1.0),
    thresh_SCEt = cms.double(4.0),
    useDynamicDPhiWindow = cms.bool(True),
    useRegression = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.patCHSMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    computeMETSignificant = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetCHS"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patCaloMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("metrawCalo"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patElectrons = cms.EDProducer("PATElectronProducer",
    addEfficiencies = cms.bool(False),
    addElectronID = cms.bool(False),
    addGenMatch = cms.bool(True),
    addMVAVariables = cms.bool(True),
    addPFClusterIso = cms.bool(False),
    addPuppiIsolation = cms.bool(False),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    computeMiniIso = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    electronIDSources = cms.PSet(

    ),
    electronSource = cms.InputTag("gedGsfElectrons"),
    embedBasicClusters = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedGsfElectronCore = cms.bool(True),
    embedGsfTrack = cms.bool(False),
    embedHighLevelSelection = cms.bool(True),
    embedPFCandidate = cms.bool(True),
    embedPflowBasicClusters = cms.bool(True),
    embedPflowPreshowerClusters = cms.bool(True),
    embedPflowSuperCluster = cms.bool(True),
    embedPreshowerClusters = cms.bool(True),
    embedRecHits = cms.bool(True),
    embedSeedCluster = cms.bool(True),
    embedSuperCluster = cms.bool(True),
    embedTrack = cms.bool(True),
    genParticleMatch = cms.InputTag("electronMatch"),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("elPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04PFIdPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04PFIdPAT"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04PFIdPAT")
    ),
    isolationValuesNoPFId = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04NoPFIdPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04NoPFIdPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04NoPFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04NoPFIdPAT"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04NoPFIdPAT")
    ),
    miniIsoParamsB = cms.vdouble(
        0.05, 0.2, 10.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0
    ),
    miniIsoParamsE = cms.vdouble(
        0.05, 0.2, 10.0, 0.0, 0.015, 
        0.015, 0.08, 0.0, 0.0
    ),
    pfCandidateMap = cms.InputTag("particleFlow","electrons"),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    pfElectronSource = cms.InputTag("particleFlow"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    useParticleFlow = cms.bool(False),
    usePfCandidateMultiMap = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patHemispheres = cms.EDProducer("PATHemisphereProducer",
    combinationMethod = cms.int32(3),
    maxElectronEta = cms.double(5),
    maxJetEta = cms.double(5),
    maxMuonEta = cms.double(5),
    maxPhotonEta = cms.double(5),
    maxTauEta = cms.double(-1),
    minElectronEt = cms.double(7),
    minJetEt = cms.double(30),
    minMuonEt = cms.double(7),
    minPhotonEt = cms.double(200000),
    minTauEt = cms.double(1000000),
    patElectrons = cms.InputTag("cleanLayer1Electrons"),
    patJets = cms.InputTag("cleanLayer1Jets"),
    patMets = cms.InputTag("layer1METs"),
    patMuons = cms.InputTag("cleanLayer1Muons"),
    patPhotons = cms.InputTag("cleanLayer1Photons"),
    patTaus = cms.InputTag("cleanLayer1Taus"),
    seedMethod = cms.int32(3)
)


process.patJetCharge = cms.EDProducer("JetChargeProducer",
    exp = cms.double(1.0),
    src = cms.InputTag("ak4JetTracksAssociatorAtVertexPF"),
    var = cms.string('Pt')
)


process.patJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(
        'L1FastJet', 
        'L2Relative', 
        'L3Absolute'
    ),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("ak4PFJetsCHS"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsReapplyJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(
        'L1FastJet', 
        'L2Relative', 
        'L3Absolute'
    ),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsReapplyJECPuppi = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(
        'L1FastJet', 
        'L2Relative', 
        'L3Absolute'
    ),
    payload = cms.string('AK4PFPuppi'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJetsPuppi"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsUpdatedJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(
        'L1FastJet', 
        'L2Relative', 
        'L3Absolute', 
        'L2L3Residual'
    ),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    jets = cms.InputTag("ak4PFJetsCHS"),
    leptons = cms.InputTag("patJetPartons","leptons"),
    partons = cms.InputTag("patJetPartons","physicsPartons"),
    rParam = cms.double(0.4)
)


process.patJetFlavourAssociationLegacy = cms.EDProducer("JetFlavourIdentifier",
    physicsDefinition = cms.bool(False),
    srcByReference = cms.InputTag("patJetPartonAssociationLegacy")
)


process.patJetGenJetMatch = cms.EDProducer("GenJetMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("ak4GenJets"),
    maxDeltaR = cms.double(0.4),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("ak4PFJetsCHS")
)


process.patJetPartonAssociationLegacy = cms.EDProducer("JetPartonMatcher",
    coneSizeToAssociate = cms.double(0.3),
    jets = cms.InputTag("ak4PFJetsCHS"),
    partons = cms.InputTag("patJetPartonsLegacy")
)


process.patJetPartonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(3.0),
    maxDeltaR = cms.double(0.4),
    mcPdgId = cms.vint32(
        1, 2, 3, 4, 5, 
        21
    ),
    mcStatus = cms.vint32(3, 23),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("ak4PFJetsCHS")
)


process.patJetPartons = cms.EDProducer("HadronAndPartonSelector",
    fullChainPhysPartons = cms.bool(True),
    particles = cms.InputTag("genParticles"),
    partonMode = cms.string('Auto'),
    src = cms.InputTag("generator")
)


process.patJetPartonsLegacy = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)


process.patJets = cms.EDProducer("PATJetProducer",
    JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociation"),
    JetPartonMapSource = cms.InputTag("patJetFlavourAssociationLegacy"),
    addAssociatedTracks = cms.bool(True),
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenPartonMatch = cms.bool(True),
    addJetCharge = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addJetFlavourInfo = cms.bool(True),
    addJetID = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(
        cms.InputTag("pfJetBProbabilityBJetTags"), cms.InputTag("pfJetProbabilityBJetTags"), cms.InputTag("pfTrackCountingHighEffBJetTags"), cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"), cms.InputTag("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"), 
        cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"), cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"), cms.InputTag("softPFMuonBJetTags"), cms.InputTag("softPFElectronBJetTags"), cms.InputTag("pfCombinedMVAV2BJetTags"), 
        cms.InputTag("pfCombinedCvsLJetTags"), cms.InputTag("pfCombinedCvsBJetTags"), cms.InputTag("pfDeepCSVJetTags","probb"), cms.InputTag("pfDeepCSVJetTags","probc"), cms.InputTag("pfDeepCSVJetTags","probudsg"), 
        cms.InputTag("pfDeepCSVJetTags","probbb")
    ),
    efficiencies = cms.PSet(

    ),
    embedGenJetMatch = cms.bool(True),
    embedGenPartonMatch = cms.bool(True),
    embedPFCandidates = cms.bool(False),
    genJetMatch = cms.InputTag("patJetGenJetMatch"),
    genPartonMatch = cms.InputTag("patJetPartonMatch"),
    getJetMCFlavour = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetCharge"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors")),
    jetIDMap = cms.InputTag("ak4JetID"),
    jetSource = cms.InputTag("ak4PFJetsCHS"),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    resolutions = cms.PSet(

    ),
    tagInfoSources = cms.VInputTag(),
    trackAssociationSource = cms.InputTag("ak4JetTracksAssociatorAtVertexPF"),
    useLegacyJetMCFlavour = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patJetsReapplyJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    printWarning = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patJetsReapplyJECPuppi = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECPuppi")),
    jetSource = cms.InputTag("slimmedJetsPuppi"),
    printWarning = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patLowPtElectrons = cms.EDProducer("PATElectronProducer",
    addEfficiencies = cms.bool(False),
    addElectronID = cms.bool(False),
    addGenMatch = cms.bool(True),
    addMVAVariables = cms.bool(False),
    addPFClusterIso = cms.bool(False),
    addPuppiIsolation = cms.bool(False),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    computeMiniIso = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    electronIDSources = cms.PSet(

    ),
    electronSource = cms.InputTag("lowPtGsfElectrons"),
    embedBasicClusters = cms.bool(False),
    embedGenMatch = cms.bool(False),
    embedGsfElectronCore = cms.bool(False),
    embedGsfTrack = cms.bool(False),
    embedHighLevelSelection = cms.bool(False),
    embedPFCandidate = cms.bool(False),
    embedPflowBasicClusters = cms.bool(False),
    embedPflowPreshowerClusters = cms.bool(False),
    embedPflowSuperCluster = cms.bool(False),
    embedPreshowerClusters = cms.bool(False),
    embedRecHits = cms.bool(False),
    embedSeedCluster = cms.bool(False),
    embedSuperCluster = cms.bool(False),
    embedTrack = cms.bool(True),
    genParticleMatch = cms.InputTag("lowPtElectronMatch"),
    isoDeposits = cms.PSet(

    ),
    isolationValues = cms.PSet(

    ),
    isolationValuesNoPFId = cms.PSet(

    ),
    miniIsoParamsB = cms.vdouble(),
    miniIsoParamsE = cms.vdouble(),
    pfCandidateMap = cms.InputTag("particleFlow","electrons"),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    pfElectronSource = cms.InputTag("particleFlow"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    useParticleFlow = cms.bool(False),
    usePfCandidateMultiMap = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patMETs = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetT1"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patMuons = cms.EDProducer("PATMuonProducer",
    addEfficiencies = cms.bool(False),
    addGenMatch = cms.bool(True),
    addPuppiIsolation = cms.bool(False),
    addResolutions = cms.bool(False),
    addTriggerMatching = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    caloMETMuonCorrs = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    computeMiniIso = cms.bool(False),
    computeMuonMVA = cms.bool(False),
    computePuppiCombinedIso = cms.bool(False),
    computeSoftMuonMVA = cms.bool(False),
    effectiveAreaVec = cms.vdouble(0.0566, 0.0562, 0.0363, 0.0119, 0.0064),
    efficiencies = cms.PSet(

    ),
    embedCaloMETMuonCorrs = cms.bool(True),
    embedCombinedMuon = cms.bool(True),
    embedDytMuon = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedHighLevelSelection = cms.bool(True),
    embedMuonBestTrack = cms.bool(True),
    embedPFCandidate = cms.bool(True),
    embedPfEcalEnergy = cms.bool(True),
    embedPickyMuon = cms.bool(True),
    embedStandAloneMuon = cms.bool(True),
    embedTcMETMuonCorrs = cms.bool(False),
    embedTpfmsMuon = cms.bool(True),
    embedTrack = cms.bool(False),
    embedTunePMuonBestTrack = cms.bool(True),
    forceBestTrackEmbedding = cms.bool(False),
    genParticleMatch = cms.InputTag("muonMatch"),
    hltCollectionFilters = cms.vstring('*'),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("muPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04PAT"),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PAT"),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PAT"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04PAT"),
        pfPhotons = cms.InputTag("muPFIsoValueGamma04PAT")
    ),
    lowPtmvaTrainingFile = cms.FileInPath('RecoMuon/MuonIdentification/data/mu_lowpt_BDTG.weights.xml'),
    miniIsoParams = cms.vdouble(
        0.05, 0.2, 10.0, 0.5, 0.0001, 
        0.01, 0.01, 0.01, 0.0
    ),
    muonSimInfo = cms.InputTag("muonSimClassifier"),
    muonSource = cms.InputTag("muons"),
    mvaDrMax = cms.double(0.4),
    mvaJetTag = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    mvaL1Corrector = cms.InputTag("ak4PFCHSL1FastjetCorrector"),
    mvaL1L2L3ResCorrector = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    mvaTrainingFile = cms.FileInPath('RecoMuon/MuonIdentification/data/mu_2017_BDTG.weights.xml'),
    mvaUseJec = cms.bool(True),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    pfMuonSource = cms.InputTag("particleFlow"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    recomputeBasicSelectors = cms.bool(True),
    resolutions = cms.PSet(

    ),
    rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    softMvaTrainingFile = cms.FileInPath('RecoMuon/MuonIdentification/data/TMVA-muonid-bmm4-B-25.weights.xml'),
    tcMETMuonCorrs = cms.InputTag("muonTCMETValueMapProducer","muCorrData"),
    triggerObjects = cms.InputTag("slimmedPatTrigger"),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    useParticleFlow = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patOOTPhotons = cms.EDProducer("PATPhotonProducer",
    addEfficiencies = cms.bool(False),
    addGenMatch = cms.bool(True),
    addPFClusterIso = cms.bool(False),
    addPhotonID = cms.bool(False),
    addPuppiIsolation = cms.bool(False),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    efficiencies = cms.PSet(

    ),
    electronSource = cms.InputTag("gedGsfElectrons"),
    embedBasicClusters = cms.bool(False),
    embedGenMatch = cms.bool(False),
    embedPreshowerClusters = cms.bool(False),
    embedRecHits = cms.bool(False),
    embedSeedCluster = cms.bool(False),
    embedSuperCluster = cms.bool(False),
    genParticleMatch = cms.InputTag("ootPhotonMatch"),
    isoDeposits = cms.PSet(

    ),
    isolationValues = cms.PSet(

    ),
    photonIDSources = cms.PSet(

    ),
    photonSource = cms.InputTag("ootPhotons"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    saveRegressionData = cms.bool(True),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patPFMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetExtractor"),
    metSource = cms.InputTag("pfMet"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag(cms.InputTag("slimmedElectrons"), cms.InputTag("slimmedMuons"), cms.InputTag("slimmedPhotons")),
    srcPFCands = cms.InputTag("packedPFCandidates"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patPFMetPuppi = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetExtractorPuppi"),
    metSource = cms.InputTag("pfMetPuppi"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFPuppi_phi'),
    srcJetResPt = cms.string('AK4PFPuppi_pt'),
    srcJetSF = cms.string('AK4PFPuppi'),
    srcJets = cms.InputTag("cleanedPatJetsPuppi"),
    srcLeptons = cms.VInputTag(cms.InputTag("slimmedElectrons"), cms.InputTag("slimmedMuons"), cms.InputTag("slimmedPhotons")),
    srcPFCands = cms.InputTag("puppiForMET"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patPFMetT0Corr = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.patPFMetT0CorrPuppi = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0Puppi"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.patPFMetT0pcT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1Puppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2CorrPuppi","type1"), cms.InputTag("patPFMetT0CorrPuppi"))
)


process.patPFMetT0pcT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Puppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2CorrPuppi","type1"), cms.InputTag("patPFMetT2CorrPuppi","type2"), cms.InputTag("patPFMetT0CorrPuppi"))
)


process.patPFMetT0pcT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"))
)


process.patPFMetT1ElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1ElectronEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDownPuppi"))
)


process.patPFMetT1ElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1ElectronEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUpPuppi"))
)


process.patPFMetT1JetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1JetEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDownPuppi"))
)


process.patPFMetT1JetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1JetEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUpPuppi"))
)


process.patPFMetT1JetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResDown"))
)


process.patPFMetT1JetResDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResDownPuppi"))
)


process.patPFMetT1JetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResUp"))
)


process.patPFMetT1JetResUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResUpPuppi"))
)


process.patPFMetT1MuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1MuonEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDownPuppi"))
)


process.patPFMetT1MuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1MuonEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUpPuppi"))
)


process.patPFMetT1PhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1PhotonEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDownPuppi"))
)


process.patPFMetT1PhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1PhotonEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUpPuppi"))
)


process.patPFMetT1Puppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2CorrPuppi","type1"))
)


process.patPFMetT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"))
)


process.patPFMetT1SmearElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1SmearElectronEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDownPuppi"))
)


process.patPFMetT1SmearElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1SmearElectronEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUpPuppi"))
)


process.patPFMetT1SmearJetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1SmearJetEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDownPuppi"))
)


process.patPFMetT1SmearJetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1SmearJetEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUpPuppi"))
)


process.patPFMetT1SmearJetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResDown"))
)


process.patPFMetT1SmearJetResDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResDownPuppi"))
)


process.patPFMetT1SmearJetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResUp"))
)


process.patPFMetT1SmearJetResUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResUpPuppi"))
)


process.patPFMetT1SmearMuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1SmearMuonEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDownPuppi"))
)


process.patPFMetT1SmearMuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1SmearMuonEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUpPuppi"))
)


process.patPFMetT1SmearPhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1SmearPhotonEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDownPuppi"))
)


process.patPFMetT1SmearPhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1SmearPhotonEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUpPuppi"))
)


process.patPFMetT1SmearPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorrPuppi","type1"))
)


process.patPFMetT1SmearTauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1SmearTauEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDownPuppi"))
)


process.patPFMetT1SmearTauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1SmearTauEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUpPuppi"))
)


process.patPFMetT1SmearUnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1SmearUnclusteredEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDownPuppi"))
)


process.patPFMetT1SmearUnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT1SmearUnclusteredEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1SmearPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUpPuppi"))
)


process.patPFMetT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"))
)


process.patPFMetT1T2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2CorrPuppi = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag(""),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Puppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2CorrPuppi","type1"), cms.InputTag("patPFMetT2CorrPuppi","type2"))
)


process.patPFMetT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"))
)


process.patPFMetT1T2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2SmearCorrPuppi = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag(""),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT1T2SmearCorrPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1TauEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDownPuppi"))
)


process.patPFMetT1TauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1TauEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUpPuppi"))
)


process.patPFMetT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TxyPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2CorrPuppi","type1"), cms.InputTag("patPFMetTxyCorrPuppi"))
)


process.patPFMetT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1UnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1UnclusteredEnDownPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDownPuppi"))
)


process.patPFMetT1UnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT1UnclusteredEnUpPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Puppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUpPuppi"))
)


process.patPFMetT2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT2CorrPuppi = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag(""),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT2SmearCorrPuppi = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT2SmearCorrPuppi"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetTxy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetTxyCorr = cms.EDProducer("MultShiftMETcorrInputProducer",
    parameters = cms.VPSet(
        cms.PSet(
            etaMax = cms.double(2.7),
            etaMin = cms.double(0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaPlus'),
            px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
            py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(0),
            etaMin = cms.double(-2.7),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaMinus'),
            px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
            py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.392),
            etaMin = cms.double(-1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0Barrel'),
            px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
            py = cms.vdouble(0.00798098092474, -0.000103998219585),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3),
            etaMin = cms.double(1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapPlus'),
            px = cms.vdouble(-0.00305719113962, -0.00032676418359),
            py = cms.vdouble(-0.00345131507897, 0.000164816815994),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.392),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapMinus'),
            px = cms.vdouble(-0.000159031461755, 0.00012231873804),
            py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.479),
            etaMin = cms.double(-1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaBarrel'),
            px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
            py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3.0),
            etaMin = cms.double(1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapPlus'),
            px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
            py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.479),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapMinus'),
            px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
            py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFPlus'),
            px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
            py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFMinus'),
            px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
            py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFPlus'),
            px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
            py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFMinus'),
            px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
            py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        )
    ),
    srcPFlow = cms.InputTag("packedPFCandidates"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.patPFMetTxyCorrPuppi = cms.EDProducer("MultShiftMETcorrInputProducer",
    parameters = cms.VPSet(
        cms.PSet(
            etaMax = cms.double(2.7),
            etaMin = cms.double(0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaPlus'),
            px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
            py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(0),
            etaMin = cms.double(-2.7),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaMinus'),
            px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
            py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.392),
            etaMin = cms.double(-1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0Barrel'),
            px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
            py = cms.vdouble(0.00798098092474, -0.000103998219585),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3),
            etaMin = cms.double(1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapPlus'),
            px = cms.vdouble(-0.00305719113962, -0.00032676418359),
            py = cms.vdouble(-0.00345131507897, 0.000164816815994),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.392),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapMinus'),
            px = cms.vdouble(-0.000159031461755, 0.00012231873804),
            py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.479),
            etaMin = cms.double(-1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaBarrel'),
            px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
            py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3.0),
            etaMin = cms.double(1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapPlus'),
            px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
            py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.479),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapMinus'),
            px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
            py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFPlus'),
            px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
            py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFMinus'),
            px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
            py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFPlus'),
            px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
            py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFMinus'),
            px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
            py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        )
    ),
    srcPFlow = cms.InputTag("packedPFCandidates"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.patPFMetTxyPuppi = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetPuppi"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetTxyCorrPuppi"))
)


process.patPhotons = cms.EDProducer("PATPhotonProducer",
    addEfficiencies = cms.bool(False),
    addGenMatch = cms.bool(True),
    addPFClusterIso = cms.bool(False),
    addPhotonID = cms.bool(False),
    addPuppiIsolation = cms.bool(False),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    efficiencies = cms.PSet(

    ),
    electronSource = cms.InputTag("gedGsfElectrons"),
    embedBasicClusters = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedPreshowerClusters = cms.bool(True),
    embedRecHits = cms.bool(True),
    embedSeedCluster = cms.bool(True),
    embedSuperCluster = cms.bool(True),
    genParticleMatch = cms.InputTag("photonMatch"),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("phPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("phPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("phPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("phPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("phPFIsoValueChargedAll04PFIdPAT"),
        pfChargedHadrons = cms.InputTag("phPFIsoValueCharged04PFIdPAT"),
        pfNeutralHadrons = cms.InputTag("phPFIsoValueNeutral04PFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("phPFIsoValuePU04PFIdPAT"),
        pfPhotons = cms.InputTag("phPFIsoValueGamma04PFIdPAT")
    ),
    photonIDSources = cms.PSet(

    ),
    photonSource = cms.InputTag("gedPhotons"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    saveRegressionData = cms.bool(True),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patSmearedJets = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(0)
)


process.patSmearedJetsPuppi = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(0)
)


process.patTaus = cms.EDProducer("PATTauProducer",
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenMatch = cms.bool(True),
    addResolutions = cms.bool(False),
    addTauID = cms.bool(True),
    addTauJetCorrFactors = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    embedGenJetMatch = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedIsolationPFCands = cms.bool(False),
    embedIsolationPFChargedHadrCands = cms.bool(False),
    embedIsolationPFGammaCands = cms.bool(False),
    embedIsolationPFNeutralHadrCands = cms.bool(False),
    embedIsolationTracks = cms.bool(False),
    embedLeadPFCand = cms.bool(False),
    embedLeadPFChargedHadrCand = cms.bool(False),
    embedLeadPFNeutralCand = cms.bool(False),
    embedLeadTrack = cms.bool(False),
    embedSignalPFCands = cms.bool(False),
    embedSignalPFChargedHadrCands = cms.bool(False),
    embedSignalPFGammaCands = cms.bool(False),
    embedSignalPFNeutralHadrCands = cms.bool(False),
    embedSignalTracks = cms.bool(False),
    genJetMatch = cms.InputTag("tauGenJetMatch"),
    genParticleMatch = cms.InputTag("tauMatch"),
    isoDeposits = cms.PSet(

    ),
    resolutions = cms.PSet(

    ),
    skipMissingTauID = cms.bool(False),
    tauIDSources = cms.PSet(
        againstElectronLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6LooseElectronRejection"),
        againstElectronMVA6Raw = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection"),
        againstElectronMVA6category = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection","category"),
        againstElectronMediumMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6MediumElectronRejection"),
        againstElectronTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6TightElectronRejection"),
        againstElectronVLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VLooseElectronRejection"),
        againstElectronVTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VTightElectronRejection"),
        againstMuonLoose3 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),
        againstMuonTight3 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
        byIsolationMVArun2v1DBdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw"),
        byIsolationMVArun2v1DBnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw"),
        byIsolationMVArun2v1DBoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw"),
        byIsolationMVArun2v1PWdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw"),
        byIsolationMVArun2v1PWnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw"),
        byIsolationMVArun2v1PWoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw"),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
        byLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT"),
        byLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT"),
        byLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT"),
        byLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT"),
        byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
        byMediumIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT"),
        byMediumIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT"),
        byMediumIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT"),
        byMediumIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT"),
        byMediumIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT"),
        byMediumIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT"),
        byPhotonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"),
        byTightCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
        byTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT"),
        byTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT"),
        byTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT"),
        byTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT"),
        byVLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT"),
        byVLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT"),
        byVLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT"),
        byVLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT"),
        byVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT"),
        byVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT"),
        byVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT"),
        byVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT"),
        byVVLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byVVLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVLooseIsolationMVArun2v1DBnewDMwLT"),
        byVVLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVLooseIsolationMVArun2v1DBoldDMwLT"),
        byVVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT"),
        byVVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT"),
        byVVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT"),
        byVVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
        chargedIsoPtSumdR03 = cms.InputTag("hpsPFTauChargedIsoPtSumdR03"),
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
        decayModeFindingNewDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
        footprintCorrection = cms.InputTag("hpsPFTauFootprintCorrection"),
        footprintCorrectiondR03 = cms.InputTag("hpsPFTauFootprintCorrectiondR03"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauNeutralIsoPtSum"),
        neutralIsoPtSumWeight = cms.InputTag("hpsPFTauNeutralIsoPtSumWeight"),
        neutralIsoPtSumWeightdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumWeightdR03"),
        neutralIsoPtSumdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumdR03"),
        photonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalCone"),
        photonPtSumOutsideSignalConedR03 = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalConedR03"),
        puCorrPtSum = cms.InputTag("hpsPFTauPUcorrPtSum")
    ),
    tauJetCorrFactorsSource = cms.VInputTag(cms.InputTag("patTauJetCorrFactors")),
    tauSource = cms.InputTag("hpsPFTauProducer"),
    tauTransverseImpactParameterSource = cms.InputTag("hpsPFTauTransverseImpactParameters"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patTrkMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    computeMETSignificant = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetTrk"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.39, 1.26, 1.21, 1.23, 1.28),
        pjpar = cms.vdouble(-0.2586, 0.6173)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.pfCandMETcorr = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag("pfCandsNotInJetsForMetCorr")
)


process.pfCandMETcorrPuppi = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag("pfCandsNotInJetsForMetCorrPuppi")
)


process.pfCandidateToVertexAssociation = cms.EDProducer("PFCand_AssoMap",
    AssociationType = cms.InputTag("Both"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ConversionsCollection = cms.InputTag("allConversions"),
    FinalAssociation = cms.untracked.int32(1),
    GetCleanedCollections = cms.bool(False),
    MaxNumberOfAssociations = cms.int32(1),
    NIVertexCollection = cms.InputTag("particleFlowDisplacedVertex"),
    PFCandidateCollection = cms.InputTag("particleFlow"),
    UseBeamSpotCompatibility = cms.untracked.bool(True),
    V0KshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    V0LambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    doReassociation = cms.bool(True),
    ignoreMissingCollection = cms.bool(True),
    nTrackWeight = cms.double(0.001)
)


process.pfCandsForUnclusteredUnc = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMuNoTau"),
    veto = cms.InputTag("slimmedPhotons")
)


process.pfCandsForUnclusteredUncPuppi = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMuNoTauPuppi"),
    veto = cms.InputTag("slimmedPhotons")
)


process.pfCandsNoJets = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("packedPFCandidates"),
    veto = cms.InputTag("cleanedPatJets")
)


process.pfCandsNoJetsNoEle = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJets"),
    veto = cms.InputTag("slimmedElectrons")
)


process.pfCandsNoJetsNoEleNoMu = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEle"),
    veto = cms.InputTag("slimmedMuons")
)


process.pfCandsNoJetsNoEleNoMuNoTau = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMu"),
    veto = cms.InputTag("slimmedTaus")
)


process.pfCandsNoJetsNoEleNoMuNoTauPuppi = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMuPuppi"),
    veto = cms.InputTag("slimmedTaus")
)


process.pfCandsNoJetsNoEleNoMuPuppi = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoElePuppi"),
    veto = cms.InputTag("slimmedMuons")
)


process.pfCandsNoJetsNoElePuppi = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsPuppi"),
    veto = cms.InputTag("slimmedElectrons")
)


process.pfCandsNoJetsPuppi = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("packedPFCandidates"),
    veto = cms.InputTag("cleanedPatJetsPuppi")
)


process.pfCandsNotInJetsForMetCorr = cms.EDProducer("PFCandidateFromFwdPtrProducer",
    src = cms.InputTag("pfCandsNotInJetsPtrForMetCorr")
)


process.pfCandsNotInJetsForMetCorrPuppi = cms.EDProducer("PFCandidateFromFwdPtrProducer",
    src = cms.InputTag("pfCandsNotInJetsPtrForMetCorr")
)


process.pfCandsNotInJetsPtrForMetCorr = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('noJet'),
    topCollection = cms.InputTag("pfJetsPtrForMetCorr"),
    verbose = cms.untracked.bool(False)
)


process.pfDeepCMVADiscriminatorsJetTags = cms.EDProducer("BTagProbabilityToDiscriminator",
    discriminators = cms.VPSet(
        cms.PSet(
            denominator = cms.VInputTag(),
            name = cms.string('BvsAll'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCMVAJetTags","probb"), cms.InputTag("pfDeepCMVAJetTags","probbb"))
        ), 
        cms.PSet(
            denominator = cms.VInputTag(cms.InputTag("pfDeepCMVAJetTags","probc"), cms.InputTag("pfDeepCMVAJetTags","probb"), cms.InputTag("pfDeepCMVAJetTags","probbb")),
            name = cms.string('CvsB'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCMVAJetTags","probc"))
        ), 
        cms.PSet(
            denominator = cms.VInputTag(cms.InputTag("pfDeepCMVAJetTags","probc"), cms.InputTag("pfDeepCMVAJetTags","probudsg")),
            name = cms.string('CvsL'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCMVAJetTags","probc"))
        )
    )
)


process.pfDeepCMVAJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/Model_DeepCMVA.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCMVATagInfos"),
    toAdd = cms.PSet(

    )
)


process.pfDeepCMVANegativeTagInfos = cms.EDProducer("DeepCMVATagInfoProducer",
    cMVAPtThreshold = cms.double(200),
    deepNNTagInfos = cms.InputTag("pfDeepCSVNegativeTagInfos"),
    elInfoSrc = cms.InputTag("softPFElectronsTagInfos"),
    ipInfoSrc = cms.InputTag("pfImpactParameterTagInfos"),
    jpComputerSrc = cms.string('candidateJetProbabilityComputer'),
    jpbComputerSrc = cms.string('candidateJetBProbabilityComputer'),
    muInfoSrc = cms.InputTag("softPFMuonsTagInfos"),
    softelComputerSrc = cms.string('softPFElectronComputer'),
    softmuComputerSrc = cms.string('softPFMuonComputer')
)


process.pfDeepCMVAPositiveTagInfos = cms.EDProducer("DeepCMVATagInfoProducer",
    cMVAPtThreshold = cms.double(200),
    deepNNTagInfos = cms.InputTag("pfDeepCSVPositiveTagInfos"),
    elInfoSrc = cms.InputTag("softPFElectronsTagInfos"),
    ipInfoSrc = cms.InputTag("pfImpactParameterTagInfos"),
    jpComputerSrc = cms.string('candidateJetProbabilityComputer'),
    jpbComputerSrc = cms.string('candidateJetBProbabilityComputer'),
    muInfoSrc = cms.InputTag("softPFMuonsTagInfos"),
    softelComputerSrc = cms.string('softPFElectronComputer'),
    softmuComputerSrc = cms.string('softPFMuonComputer')
)


process.pfDeepCMVATagInfos = cms.EDProducer("DeepCMVATagInfoProducer",
    cMVAPtThreshold = cms.double(200),
    deepNNTagInfos = cms.InputTag("pfDeepCSVTagInfos"),
    elInfoSrc = cms.InputTag("softPFElectronsTagInfos"),
    ipInfoSrc = cms.InputTag("pfImpactParameterTagInfos"),
    jpComputerSrc = cms.string('candidateJetProbabilityComputer'),
    jpbComputerSrc = cms.string('candidateJetBProbabilityComputer'),
    muInfoSrc = cms.InputTag("softPFMuonsTagInfos"),
    softelComputerSrc = cms.string('softPFElectronComputer'),
    softmuComputerSrc = cms.string('softPFMuonComputer')
)


process.pfDeepCSVDiscriminatorsJetTags = cms.EDProducer("BTagProbabilityToDiscriminator",
    discriminators = cms.VPSet(
        cms.PSet(
            denominator = cms.VInputTag(),
            name = cms.string('BvsAll'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags","probb"), cms.InputTag("pfDeepCSVJetTags","probbb"))
        ), 
        cms.PSet(
            denominator = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags","probc"), cms.InputTag("pfDeepCSVJetTags","probb"), cms.InputTag("pfDeepCSVJetTags","probbb")),
            name = cms.string('CvsB'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags","probc"))
        ), 
        cms.PSet(
            denominator = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags","probudsg"), cms.InputTag("pfDeepCSVJetTags","probc")),
            name = cms.string('CvsL'),
            numerator = cms.VInputTag(cms.InputTag("pfDeepCSVJetTags","probc"))
        )
    )
)


process.pfDeepCSVJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourNoSL.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCSVTagInfos"),
    toAdd = cms.PSet(
        probcc = cms.string('probc')
    )
)


process.pfDeepCSVNegativeTagInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(True),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(-2.0),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(10.0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(10.0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(True)
    ),
    svTagInfos = cms.InputTag("pfInclusiveSecondaryVertexFinderNegativeTagInfos")
)


process.pfDeepCSVPositiveTagInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfos")
)


process.pfDeepCSVTagInfos = cms.EDProducer("DeepNNTagInfoProducer",
    computer = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(True),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    svTagInfos = cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfos")
)


process.pfImpactParameterTagInfos = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("particleFlow"),
    computeGhostTrack = cms.bool(True),
    computeProbabilities = cms.bool(True),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("ak4PFJetsCHS"),
    maxDeltaR = cms.double(0.4),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(0),
    minimumNumberOfPixelHits = cms.int32(1),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    useTrackQuality = cms.bool(False)
)


process.pfInclusiveSecondaryVertexFinderTagInfos = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("inclusiveCandidateSecondaryVertices"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(2.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.79),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.pfJetsPtrForMetCorr = cms.EDProducer("PFJetFwdPtrProducer",
    src = cms.InputTag("ak4PFJets")
)


process.pfMETcorrType0 = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.pfMet = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('raw'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.pfMetCHS = cms.EDProducer("PFMETProducer",
    alias = cms.string('pfMet'),
    calculateSignificance = cms.bool(False),
    globalThreshold = cms.double(0.0),
    src = cms.InputTag("pfCHS")
)


process.pfMetPuppi = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('raw'),
    metSource = cms.InputTag("slimmedMETsPuppi","","@skipCurrentProcess")
)


process.pfMetT1 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"))
)


process.pfMetT1T2 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"))
)


process.pfMetTrk = cms.EDProducer("PFMETProducer",
    alias = cms.string('pfMet'),
    calculateSignificance = cms.bool(False),
    globalThreshold = cms.double(0.0),
    src = cms.InputTag("pfTrk")
)


process.pfNegativeDeepCMVAJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/Model_DeepCMVA.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCMVANegativeTagInfos"),
    toAdd = cms.PSet(

    )
)


process.pfNegativeDeepCSVJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourNoSL.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCSVNegativeTagInfos"),
    toAdd = cms.PSet(
        probcc = cms.string('probc')
    )
)


process.pfNoPileUpIsoPFBRECO = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    topCollection = cms.InputTag("pfPileUpIsoPFBRECO"),
    verbose = cms.untracked.bool(False)
)


process.pfNoPileUpPFBRECO = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    topCollection = cms.InputTag("pfPileUpPFBRECO"),
    verbose = cms.untracked.bool(False)
)


process.pfPileUpIsoPFBRECO = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlowPtrs"),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.pfPileUpPFBRECO = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlowPtrs"),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.pfPositiveDeepCMVAJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/Model_DeepCMVA.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCMVAPositiveTagInfos"),
    toAdd = cms.PSet(

    )
)


process.pfPositiveDeepCSVJetTags = cms.EDProducer("DeepFlavourJetTagsProducer",
    NNConfig = cms.FileInPath('RecoBTag/Combined/data/DeepFlavourNoSL.json'),
    checkSVForDefaults = cms.bool(False),
    meanPadding = cms.bool(False),
    src = cms.InputTag("pfDeepCSVPositiveTagInfos"),
    toAdd = cms.PSet(
        probcc = cms.string('probc')
    )
)


process.pfSecondaryVertexTagInfos = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("secondaryVertices"),
    extSVDeltaRToJet = cms.double(0.3),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("pfImpactParameterTagInfos"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(False),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(3.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.65),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.4),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.phPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(False),
        SCMatch_Veto = cms.bool(True),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(0.4),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoValueCharged03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.photonCore = cms.EDProducer("PhotonCoreProducer",
    conversionProducer = cms.InputTag(""),
    minSCEt = cms.double(10.0),
    photonCoreCollection = cms.string(''),
    pixelSeedProducer = cms.InputTag("electronMergedSeeds"),
    risolveConversionAmbiguity = cms.bool(True),
    scHybridBarrelProducer = cms.InputTag("correctedHybridSuperClusters"),
    scIslandEndcapProducer = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower")
)


process.photonCoreFromMultiCl = cms.EDProducer("PhotonCoreProducer",
    conversionProducer = cms.InputTag(""),
    minSCEt = cms.double(10.0),
    photonCoreCollection = cms.string(''),
    pixelSeedProducer = cms.InputTag("electronMergedSeedsFromMultiCl"),
    risolveConversionAmbiguity = cms.bool(True),
    scHybridBarrelProducer = cms.InputTag("correctedHybridSuperClusters"),
    scIslandEndcapProducer = cms.InputTag("particleFlowSuperClusterHGCalFromMultiCl")
)


process.photonIDValueMapProducer = cms.EDProducer("PhotonIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
    pfCandidates = cms.InputTag("particleFlow"),
    pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag(""),
    srcMiniAOD = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.photonMVAValueMapProducer = cms.EDProducer("PhotonMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(
        cms.PSet(
            categoryCuts = cms.vstring(
                'abs(superCluster.eta) <  1.479', 
                'abs(superCluster.eta) >= 1.479'
            ),
            effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased_3bins.txt'),
            mvaName = cms.string('PhotonMVAEstimator'),
            mvaTag = cms.string('Run2Spring16NonTrigV1'),
            nCategories = cms.int32(2),
            phoIsoCutoff = cms.double(2.5),
            phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
            variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesSpring16.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EB_V1.weights.xml.gz', 
                'RecoEgamma/PhotonIdentification/data/MVA/Spring16/EE_V1.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'abs(superCluster.eta) <  1.479', 
                'abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('PhotonMVAEstimator'),
            mvaTag = cms.string('RunIIFall17v1p1'),
            nCategories = cms.int32(2),
            variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V1.weights.xml.gz', 
                'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V1.weights.xml.gz'
            )
        ), 
        cms.PSet(
            categoryCuts = cms.vstring(
                'abs(superCluster.eta) <  1.479', 
                'abs(superCluster.eta) >= 1.479'
            ),
            mvaName = cms.string('PhotonMVAEstimator'),
            mvaTag = cms.string('RunIIFall17v2'),
            nCategories = cms.int32(2),
            variableDefinition = cms.string('RecoEgamma/PhotonIdentification/data/PhotonMVAEstimatorRun2VariablesFall17V1p1.txt'),
            weightFileNames = cms.vstring(
                'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EB_V2.weights.xml.gz', 
                'RecoEgamma/PhotonIdentification/data/MVA/Fall17/EE_V2.weights.xml.gz'
            )
        )
    ),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
)


process.photonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(1.0),
    maxDeltaR = cms.double(0.2),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("gedPhotons")
)


process.puppi = cms.EDProducer("PuppiProducer",
    DeltaZCut = cms.double(0.3),
    MinPuppiWeight = cms.double(0.01),
    PtMaxNeutrals = cms.double(200.0),
    UseDeltaZCut = cms.bool(True),
    algos = cms.VPSet(
        cms.PSet(
            EtaMaxExtrap = cms.double(2.0),
            MedEtaSF = cms.vdouble(1.0),
            MinNeutralPt = cms.vdouble(0.2),
            MinNeutralPtSlope = cms.vdouble(0.015),
            RMSEtaSF = cms.vdouble(1.0),
            etaMax = cms.vdouble(2.5),
            etaMin = cms.vdouble(0.0),
            ptMin = cms.vdouble(0.0),
            puppiAlgos = cms.VPSet(cms.PSet(
                algoId = cms.int32(5),
                applyLowPUCorr = cms.bool(True),
                combOpt = cms.int32(0),
                cone = cms.double(0.4),
                rmsPtMin = cms.double(0.1),
                rmsScaleFactor = cms.double(1.0),
                useCharged = cms.bool(True)
            ))
        ), 
        cms.PSet(
            EtaMaxExtrap = cms.double(2.0),
            MedEtaSF = cms.vdouble(0.9, 0.75),
            MinNeutralPt = cms.vdouble(1.7, 2.0),
            MinNeutralPtSlope = cms.vdouble(0.08, 0.08),
            RMSEtaSF = cms.vdouble(1.2, 0.95),
            etaMax = cms.vdouble(3.0, 10.0),
            etaMin = cms.vdouble(2.5, 3.0),
            ptMin = cms.vdouble(0.0, 0.0),
            puppiAlgos = cms.VPSet(cms.PSet(
                algoId = cms.int32(5),
                applyLowPUCorr = cms.bool(True),
                combOpt = cms.int32(0),
                cone = cms.double(0.4),
                rmsPtMin = cms.double(0.5),
                rmsScaleFactor = cms.double(1.0),
                useCharged = cms.bool(False)
            ))
        )
    ),
    applyCHS = cms.bool(True),
    candName = cms.InputTag("packedPFCandidates"),
    clonePackedCands = cms.bool(True),
    invertPuppi = cms.bool(False),
    puppiDiagnostics = cms.bool(False),
    puppiForLeptons = cms.bool(False),
    useExistingWeights = cms.bool(True),
    useExp = cms.bool(False),
    useWeightsNoLep = cms.bool(False),
    vertexName = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vtxNdofCut = cms.int32(4),
    vtxZCut = cms.double(24)
)


process.puppiForMET = cms.EDProducer("PuppiPhoton",
    candName = cms.InputTag("packedPFCandidates"),
    dRMatch = cms.vdouble(0.005, 0.005, 0.005, 0.005),
    eta = cms.double(2.5),
    pdgids = cms.vint32(22, 11, 211, 130),
    photonId = cms.InputTag(""),
    photonName = cms.InputTag("slimmedPhotons"),
    pt = cms.double(20),
    puppiCandName = cms.InputTag("puppiMerged"),
    recoToPFMap = cms.InputTag("reducedEgamma","reducedPhotonPfCandMap"),
    runOnMiniAOD = cms.bool(True),
    usePFphotons = cms.bool(True),
    useRefs = cms.bool(False),
    useValueMap = cms.bool(False),
    weight = cms.double(1.0),
    weightsName = cms.InputTag("puppi")
)


process.puppiMerged = cms.EDProducer("CandViewMerger",
    src = cms.VInputTag("puppiNoLep", "pfLeptonsPUPPET")
)


process.puppiNoLep = cms.EDProducer("PuppiProducer",
    DeltaZCut = cms.double(0.3),
    MinPuppiWeight = cms.double(0.01),
    PtMaxNeutrals = cms.double(200.0),
    UseDeltaZCut = cms.bool(True),
    algos = cms.VPSet(
        cms.PSet(
            EtaMaxExtrap = cms.double(2.0),
            MedEtaSF = cms.vdouble(1.0),
            MinNeutralPt = cms.vdouble(0.2),
            MinNeutralPtSlope = cms.vdouble(0.015),
            RMSEtaSF = cms.vdouble(1.0),
            etaMax = cms.vdouble(2.5),
            etaMin = cms.vdouble(0.0),
            ptMin = cms.vdouble(0.0),
            puppiAlgos = cms.VPSet(cms.PSet(
                algoId = cms.int32(5),
                applyLowPUCorr = cms.bool(True),
                combOpt = cms.int32(0),
                cone = cms.double(0.4),
                rmsPtMin = cms.double(0.1),
                rmsScaleFactor = cms.double(1.0),
                useCharged = cms.bool(True)
            ))
        ), 
        cms.PSet(
            EtaMaxExtrap = cms.double(2.0),
            MedEtaSF = cms.vdouble(0.9, 0.75),
            MinNeutralPt = cms.vdouble(1.7, 2.0),
            MinNeutralPtSlope = cms.vdouble(0.08, 0.08),
            RMSEtaSF = cms.vdouble(1.2, 0.95),
            etaMax = cms.vdouble(3.0, 10.0),
            etaMin = cms.vdouble(2.5, 3.0),
            ptMin = cms.vdouble(0.0, 0.0),
            puppiAlgos = cms.VPSet(cms.PSet(
                algoId = cms.int32(5),
                applyLowPUCorr = cms.bool(True),
                combOpt = cms.int32(0),
                cone = cms.double(0.4),
                rmsPtMin = cms.double(0.5),
                rmsScaleFactor = cms.double(1.0),
                useCharged = cms.bool(False)
            ))
        )
    ),
    applyCHS = cms.bool(True),
    candName = cms.InputTag("pfNoLepPUPPI"),
    clonePackedCands = cms.bool(True),
    invertPuppi = cms.bool(False),
    puppiDiagnostics = cms.bool(False),
    puppiForLeptons = cms.bool(False),
    useExistingWeights = cms.bool(True),
    useExp = cms.bool(False),
    useWeightsNoLep = cms.bool(True),
    vertexName = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vtxNdofCut = cms.int32(4),
    vtxZCut = cms.double(24)
)


process.puppiPhoton = cms.EDProducer("PuppiPhoton",
    candName = cms.InputTag("particleFlow"),
    dRMatch = cms.vdouble(0.005, 0.005, 0.005, 0.005),
    eta = cms.double(2.5),
    pdgids = cms.vint32(22, 11, 211, 130),
    photonId = cms.InputTag(""),
    photonName = cms.InputTag("reducedEgamma","reducedGedPhotons"),
    pt = cms.double(20),
    puppiCandName = cms.InputTag("puppi"),
    recoToPFMap = cms.InputTag("reducedEgamma","reducedPhotonPfCandMap"),
    runOnMiniAOD = cms.bool(False),
    usePFphotons = cms.bool(True),
    useRefs = cms.bool(True),
    useValueMap = cms.bool(False),
    weight = cms.double(1.0),
    weightsName = cms.InputTag("puppi")
)


process.shiftedPatElectronEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatElectronEnDownPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfElectronsPuppi"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatElectronEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatElectronEnUpPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfElectronsPuppi"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatJetEnDown = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetEnDownPuppi = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFPuppiL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFPuppiL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFPuppi'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("cleanedPatJetsPuppi")
)


process.shiftedPatJetEnUp = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetEnUpPuppi = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFPuppiL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFPuppiL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFPuppi'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("cleanedPatJetsPuppi")
)


process.shiftedPatJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-101)
)


process.shiftedPatJetResDownPuppi = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFPuppi'),
    algopt = cms.string('AK4PFPuppi_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-101)
)


process.shiftedPatJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(101)
)


process.shiftedPatJetResUpPuppi = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFPuppi'),
    algopt = cms.string('AK4PFPuppi_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(101)
)


process.shiftedPatMETCorrElectronEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnDown")
)


process.shiftedPatMETCorrElectronEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfElectronsPuppi"),
    srcShifted = cms.InputTag("shiftedPatElectronEnDownPuppi")
)


process.shiftedPatMETCorrElectronEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnUp")
)


process.shiftedPatMETCorrElectronEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfElectronsPuppi"),
    srcShifted = cms.InputTag("shiftedPatElectronEnUpPuppi")
)


process.shiftedPatMETCorrJetEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnDown")
)


process.shiftedPatMETCorrJetEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJetsPuppi"),
    srcShifted = cms.InputTag("shiftedPatJetEnDownPuppi")
)


process.shiftedPatMETCorrJetEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnUp")
)


process.shiftedPatMETCorrJetEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJetsPuppi"),
    srcShifted = cms.InputTag("shiftedPatJetEnUpPuppi")
)


process.shiftedPatMETCorrJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResDown")
)


process.shiftedPatMETCorrJetResDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJetsPuppi"),
    srcShifted = cms.InputTag("shiftedPatJetResDownPuppi")
)


process.shiftedPatMETCorrJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResUp")
)


process.shiftedPatMETCorrJetResUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJetsPuppi"),
    srcShifted = cms.InputTag("shiftedPatJetResUpPuppi")
)


process.shiftedPatMETCorrMuonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnDown")
)


process.shiftedPatMETCorrMuonEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfMuonsPuppi"),
    srcShifted = cms.InputTag("shiftedPatMuonEnDownPuppi")
)


process.shiftedPatMETCorrMuonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnUp")
)


process.shiftedPatMETCorrMuonEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfMuonsPuppi"),
    srcShifted = cms.InputTag("shiftedPatMuonEnUpPuppi")
)


process.shiftedPatMETCorrPhotonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnDown")
)


process.shiftedPatMETCorrPhotonEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfPhotonsPuppi"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnDownPuppi")
)


process.shiftedPatMETCorrPhotonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnUp")
)


process.shiftedPatMETCorrPhotonEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfPhotonsPuppi"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnUpPuppi")
)


process.shiftedPatMETCorrSmearedJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResDown")
)


process.shiftedPatMETCorrSmearedJetResDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorrPuppi"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResDownPuppi")
)


process.shiftedPatMETCorrSmearedJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResUp")
)


process.shiftedPatMETCorrSmearedJetResUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorrPuppi"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResUpPuppi")
)


process.shiftedPatMETCorrTauEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnDown")
)


process.shiftedPatMETCorrTauEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfTausPuppi"),
    srcShifted = cms.InputTag("shiftedPatTauEnDownPuppi")
)


process.shiftedPatMETCorrTauEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnUp")
)


process.shiftedPatMETCorrTauEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfTausPuppi"),
    srcShifted = cms.InputTag("shiftedPatTauEnUpPuppi")
)


process.shiftedPatMETCorrUnclusteredEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnDown")
)


process.shiftedPatMETCorrUnclusteredEnDownPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUncPuppi"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnDownPuppi")
)


process.shiftedPatMETCorrUnclusteredEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnUp")
)


process.shiftedPatMETCorrUnclusteredEnUpPuppi = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUncPuppi"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnUpPuppi")
)


process.shiftedPatMuonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatMuonEnDownPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfMuonsPuppi"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatMuonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatMuonEnUpPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfMuonsPuppi"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatPhotonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatPhotonEnDownPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfPhotonsPuppi"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatPhotonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatPhotonEnUpPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfPhotonsPuppi"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatSmearedJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-1)
)


process.shiftedPatSmearedJetResDownPuppi = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFPuppi'),
    algopt = cms.string('AK4PFPuppi_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-1)
)


process.shiftedPatSmearedJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(1)
)


process.shiftedPatSmearedJetResUpPuppi = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFPuppi'),
    algopt = cms.string('AK4PFPuppi_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJetsPuppi"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(1)
)


process.shiftedPatTauEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatTauEnDownPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfTausPuppi"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatTauEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatTauEnUpPuppi = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfTausPuppi"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatUnclusteredEnDown = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(
        cms.PSet(
            binSelection = cms.string('charge!=0'),
            binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )
    ),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.shiftedPatUnclusteredEnDownPuppi = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(
        cms.PSet(
            binSelection = cms.string('charge!=0'),
            binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )
    ),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfCandsForUnclusteredUncPuppi")
)


process.shiftedPatUnclusteredEnUp = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(
        cms.PSet(
            binSelection = cms.string('charge!=0'),
            binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )
    ),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.shiftedPatUnclusteredEnUpPuppi = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(
        cms.PSet(
            binSelection = cms.string('charge!=0'),
            binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )
    ),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfCandsForUnclusteredUncPuppi")
)


process.slimmedElectrons = cms.EDProducer("ModifiedElectronProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(
            cms.PSet(
                electron_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPreCorr"),
                    ecalTrkEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPostCorr"),
                    ecalTrkEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPreCorr"),
                    ecalTrkEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPostCorr"),
                    ecalTrkEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPreCorr"),
                    electronSrc = cms.InputTag("updatedElectrons"),
                    energyScaleDown = cms.InputTag("calibratedPatElectrons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatElectrons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatElectrons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatElectrons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatElectrons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatElectrons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatElectrons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatElectrons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatElectrons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatElectrons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatElectrons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatElectrons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatElectrons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatElectrons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatElectrons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatElectrons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatElectrons","energySmearNrSigma")
                ),
                modifierName = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
                overrideExistingValues = cms.bool(True),
                photon_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPreCorr"),
                    energyScaleDown = cms.InputTag("calibratedPatPhotons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatPhotons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatPhotons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatPhotons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatPhotons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatPhotons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatPhotons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatPhotons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatPhotons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatPhotons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatPhotons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatPhotons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatPhotons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatPhotons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatPhotons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatPhotons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatPhotons","energySmearNrSigma"),
                    photonSrc = cms.InputTag("updatedPhotons")
                )
            ), 
            cms.PSet(
                epCombConfig = cms.PSet(
                    ecalTrkRegressionConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMaxHighEt = cms.double(3.0),
                        rangeMaxLowEt = cms.double(3.0),
                        rangeMinHighEt = cms.double(-1.0),
                        rangeMinLowEt = cms.double(-1.0)
                    ),
                    ecalTrkRegressionUncertConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMaxHighEt = cms.double(0.5),
                        rangeMaxLowEt = cms.double(0.5),
                        rangeMinHighEt = cms.double(0.0002),
                        rangeMinLowEt = cms.double(0.0002)
                    ),
                    maxEPDiffInSigmaForComb = cms.double(15.0),
                    maxEcalEnergyForComb = cms.double(200.0),
                    maxRelTrkMomErrForComb = cms.double(10.0),
                    minEOverPForComb = cms.double(0.025)
                ),
                modifierName = cms.string('EGEtScaleSysModifier'),
                overrideExistingValues = cms.bool(True),
                uncertFunc = cms.PSet(
                    highEt = cms.double(46.5),
                    highEtUncert = cms.double(-0.002),
                    lowEt = cms.double(43.5),
                    lowEtUncert = cms.double(0.002),
                    name = cms.string('UncertFuncV1')
                )
            )
        )
    ),
    src = cms.InputTag("updatedElectrons")
)


process.slimmedMETs = cms.EDProducer("PATMETSlimmer",
    caloMET = cms.InputTag("patCaloMet"),
    chsMET = cms.InputTag("patCHSMet"),
    rawVariation = cms.InputTag("patPFMet"),
    runningOnMiniAOD = cms.bool(True),
    src = cms.InputTag("patPFMetT1"),
    t01Variation = cms.InputTag("slimmedMETs","","@skipCurrentProcess"),
    t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%s"),
    t1Uncertainties = cms.InputTag("patPFMetT1%s"),
    tXYUncForT1 = cms.InputTag("patPFMetT1Txy"),
    trkMET = cms.InputTag("patTrkMet")
)


process.slimmedMETsPuppi = cms.EDProducer("PATMETSlimmer",
    caloMET = cms.InputTag("patCaloMet"),
    chsMET = cms.InputTag("patCHSMet"),
    rawVariation = cms.InputTag("patPFMetPuppi"),
    runningOnMiniAOD = cms.bool(True),
    src = cms.InputTag("patPFMetT1Puppi"),
    t01Variation = cms.InputTag("slimmedMETsPuppi","","@skipCurrentProcess"),
    t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%sPuppi"),
    t1Uncertainties = cms.InputTag("patPFMetT1%sPuppi"),
    tXYUncForT1 = cms.InputTag("patPFMetT1TxyPuppi"),
    trkMET = cms.InputTag("patTrkMet")
)


process.slimmedPhotons = cms.EDProducer("ModifiedPhotonProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(
            cms.PSet(
                electron_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPreCorr"),
                    ecalTrkEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPostCorr"),
                    ecalTrkEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPreCorr"),
                    ecalTrkEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPostCorr"),
                    ecalTrkEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPreCorr"),
                    electronSrc = cms.InputTag("updatedElectrons"),
                    energyScaleDown = cms.InputTag("calibratedPatElectrons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatElectrons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatElectrons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatElectrons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatElectrons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatElectrons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatElectrons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatElectrons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatElectrons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatElectrons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatElectrons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatElectrons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatElectrons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatElectrons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatElectrons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatElectrons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatElectrons","energySmearNrSigma")
                ),
                modifierName = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
                overrideExistingValues = cms.bool(True),
                photon_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPreCorr"),
                    energyScaleDown = cms.InputTag("calibratedPatPhotons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatPhotons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatPhotons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatPhotons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatPhotons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatPhotons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatPhotons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatPhotons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatPhotons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatPhotons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatPhotons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatPhotons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatPhotons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatPhotons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatPhotons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatPhotons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatPhotons","energySmearNrSigma"),
                    photonSrc = cms.InputTag("updatedPhotons")
                )
            ), 
            cms.PSet(
                epCombConfig = cms.PSet(
                    ecalTrkRegressionConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMaxHighEt = cms.double(3.0),
                        rangeMaxLowEt = cms.double(3.0),
                        rangeMinHighEt = cms.double(-1.0),
                        rangeMinLowEt = cms.double(-1.0)
                    ),
                    ecalTrkRegressionUncertConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMaxHighEt = cms.double(0.5),
                        rangeMaxLowEt = cms.double(0.5),
                        rangeMinHighEt = cms.double(0.0002),
                        rangeMinLowEt = cms.double(0.0002)
                    ),
                    maxEPDiffInSigmaForComb = cms.double(15.0),
                    maxEcalEnergyForComb = cms.double(200.0),
                    maxRelTrkMomErrForComb = cms.double(10.0),
                    minEOverPForComb = cms.double(0.025)
                ),
                modifierName = cms.string('EGEtScaleSysModifier'),
                overrideExistingValues = cms.bool(True),
                uncertFunc = cms.PSet(
                    highEt = cms.double(46.5),
                    highEtUncert = cms.double(-0.002),
                    lowEt = cms.double(43.5),
                    lowEtUncert = cms.double(0.002),
                    name = cms.string('UncertFuncV1')
                )
            )
        )
    ),
    src = cms.InputTag("updatedPhotons")
)


process.slimmedTausNewID = cms.EDProducer("PATTauIDEmbedder",
    src = cms.InputTag("slimmedTaus"),
    tauIDSources = cms.PSet(
        byDeepTau2017v2p1VSeraw = cms.InputTag("deepTau2017v2p1","VSe"),
        byDeepTau2017v2p1VSjetraw = cms.InputTag("deepTau2017v2p1","VSjet"),
        byDeepTau2017v2p1VSmuraw = cms.InputTag("deepTau2017v2p1","VSmu"),
        byLooseDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeLoose"),
        byLooseDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetLoose"),
        byLooseDeepTau2017v2p1VSmu = cms.InputTag("deepTau2017v2p1","VSmuLoose"),
        byMediumDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeMedium"),
        byMediumDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetMedium"),
        byMediumDeepTau2017v2p1VSmu = cms.InputTag("deepTau2017v2p1","VSmuMedium"),
        byTightDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeTight"),
        byTightDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetTight"),
        byTightDeepTau2017v2p1VSmu = cms.InputTag("deepTau2017v2p1","VSmuTight"),
        byVLooseDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeVLoose"),
        byVLooseDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetVLoose"),
        byVLooseDeepTau2017v2p1VSmu = cms.InputTag("deepTau2017v2p1","VSmuVLoose"),
        byVTightDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeVTight"),
        byVTightDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetVTight"),
        byVVLooseDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeVVLoose"),
        byVVLooseDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetVVLoose"),
        byVVTightDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeVVTight"),
        byVVTightDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetVVTight"),
        byVVVLooseDeepTau2017v2p1VSe = cms.InputTag("deepTau2017v2p1","VSeVVVLoose"),
        byVVVLooseDeepTau2017v2p1VSjet = cms.InputTag("deepTau2017v2p1","VSjetVVVLoose")
    )
)


process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("tauGenJetsSelectorAllHadrons"),
    maxDPtRel = cms.double(3.0),
    maxDeltaR = cms.double(0.1),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer")
)


process.tauGenJets = cms.EDProducer("TauGenJetProducer",
    GenParticles = cms.InputTag("genParticles"),
    includeNeutrinos = cms.bool(False),
    verbose = cms.untracked.bool(False)
)


process.tauIsoDepositPFCandidates = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("particleFlow"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFChargedHadrons = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(0.1),
        Diff_z = cms.double(0.2),
        candidateSource = cms.InputTag("pfAllChargedHadronsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFGammas = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("pfAllPhotonsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFNeutralHadrons = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("pfAllNeutralHadronsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(999.9),
    maxDeltaR = cms.double(999.9),
    mcPdgId = cms.vint32(15),
    mcStatus = cms.vint32(2),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer")
)


process.updatedElectrons = cms.EDProducer("ModifiedElectronProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(cms.PSet(
            allowGsfTrackForConvs = cms.bool(False),
            beamspot = cms.InputTag("offlineBeamSpot"),
            conversions = cms.InputTag("reducedEgamma","reducedConversions"),
            ecalRecHitsEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
            ecalRecHitsEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
            eleCollVMsAreKeyedTo = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
            eleTrkIso = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
            eleTrkIso04 = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso04"),
            modifierName = cms.string('EG9X105XObjectUpdateModifier'),
            phoChargedHadIso = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoChargedHadPFPVIso = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
            phoChargedHadWorstVtxConeVetoIso = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationConeVeto"),
            phoChargedHadWorstVtxIso = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
            phoCollVMsAreKeyedTo = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
            phoNeutralHadIso = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
            phoPhotonIso = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            updateChargedHadPFPVIso = cms.bool(True)
        ))
    ),
    src = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.updatedPatJetsUpdatedJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(False),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsUpdatedJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    printWarning = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag()
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.updatedPhotons = cms.EDProducer("ModifiedPhotonProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(cms.PSet(
            allowGsfTrackForConvs = cms.bool(False),
            beamspot = cms.InputTag("offlineBeamSpot"),
            conversions = cms.InputTag("reducedEgamma","reducedConversions"),
            ecalRecHitsEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
            ecalRecHitsEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
            eleCollVMsAreKeyedTo = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
            eleTrkIso = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
            eleTrkIso04 = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso04"),
            modifierName = cms.string('EG9X105XObjectUpdateModifier'),
            phoChargedHadIso = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoChargedHadPFPVIso = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
            phoChargedHadWorstVtxConeVetoIso = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationConeVeto"),
            phoChargedHadWorstVtxIso = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
            phoCollVMsAreKeyedTo = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
            phoNeutralHadIso = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
            phoPhotonIso = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            updateChargedHadPFPVIso = cms.bool(True)
        ))
    ),
    src = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
)


process.countPatElectrons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    src = cms.InputTag("cleanPatElectrons")
)


process.countPatJets = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    src = cms.InputTag("cleanPatJets")
)


process.countPatLeptons = cms.EDFilter("PATLeptonCountFilter",
    countElectrons = cms.bool(True),
    countMuons = cms.bool(True),
    countTaus = cms.bool(False),
    electronSource = cms.InputTag("cleanPatElectrons"),
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    muonSource = cms.InputTag("cleanPatMuons"),
    tauSource = cms.InputTag("cleanPatTaus")
)


process.countPatMuons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    src = cms.InputTag("cleanPatMuons")
)


process.countPatPhotons = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    src = cms.InputTag("cleanPatPhotons")
)


process.countPatTaus = cms.EDFilter("PATCandViewCountFilter",
    maxNumber = cms.uint32(999999),
    minNumber = cms.uint32(0),
    src = cms.InputTag("cleanPatTaus")
)


process.ecalBadCalibFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
    baddetEcal = cms.vuint32(),
    debug = cms.bool(False),
    ecalMinEt = cms.double(50.0),
    taggingMode = cms.bool(False)
)


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits"),
    baddetEcal = cms.vuint32(
        872439604, 872422825, 872420274, 872423218, 872423215, 
        872416066, 872435036, 872439336, 872420273, 872436907, 
        872420147, 872439731, 872436657, 872420397, 872439732, 
        872439339, 872439603, 872422436, 872439861, 872437051, 
        872437052, 872420649, 872421950, 872437185, 872422564, 
        872421566, 872421695, 872421955, 872421567, 872437184, 
        872421951, 872421694, 872437056, 872437057, 872437313, 
        872438182, 872438951, 872439990, 872439864, 872439609, 
        872437181, 872437182, 872437053, 872436794, 872436667, 
        872436536, 872421541, 872421413, 872421414, 872421031, 
        872423083, 872421439
    ),
    debug = cms.bool(False),
    ecalMinEt = cms.double(50.0),
    taggingMode = cms.bool(True)
)


process.jetSelectorForMet = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt>15 && abs(eta)<9.9'),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("basicJetsForMet")
)


process.jetSelectorForMetPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt>15 && abs(eta)<9.9'),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("basicJetsForMetPuppi")
)


process.pfAllChargedHadronsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(
        211, -211, 321, -321, 999211, 
        2212, -2212
    ),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllChargedParticlesPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(
        211, -211, 321, -321, 999211, 
        2212, -2212, 11, -11, 13, 
        -13
    ),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllNeutralHadronsAndPhotonsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(22, 111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllNeutralHadronsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllPhotonsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(22),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfCHS = cms.EDFilter("CandPtrSelector",
    cut = cms.string('fromPV(0)>0'),
    src = cms.InputTag("packedPFCandidates")
)


process.pfElectrons = cms.EDFilter("CandPtrSelector",
    cut = cms.string("pt > 5 && isPF && gsfTrack.isAvailable() && gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2"),
    src = cms.InputTag("slimmedElectrons")
)


process.pfElectronsPuppi = cms.EDFilter("CandPtrSelector",
    cut = cms.string("pt > 5 && isPF && gsfTrack.isAvailable() && gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2"),
    src = cms.InputTag("slimmedElectrons")
)


process.pfLeptonsPUPPET = cms.EDFilter("CandPtrSelector",
    cut = cms.string('abs(pdgId) == 13 || abs(pdgId) == 11 || abs(pdgId) == 15'),
    src = cms.InputTag("packedPFCandidates")
)


process.pfMuons = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 5.0 && isPFMuon && abs(eta) < 2.4'),
    src = cms.InputTag("slimmedMuons")
)


process.pfMuonsPuppi = cms.EDFilter("CandPtrSelector",
    cut = cms.string('pt > 5.0 && isPFMuon && abs(eta) < 2.4'),
    src = cms.InputTag("slimmedMuons")
)


process.pfNoLepPUPPI = cms.EDFilter("CandPtrSelector",
    cut = cms.string('abs(pdgId) != 13 && abs(pdgId) != 11 && abs(pdgId) != 15'),
    src = cms.InputTag("packedPFCandidates")
)


process.pfNoPileUp = cms.EDFilter("CandPtrSelector",
    cut = cms.string('fromPV > 1'),
    src = cms.InputTag("packedPFCandidates")
)


process.pfNoPileUpPuppi = cms.EDFilter("CandPtrSelector",
    cut = cms.string(''),
    src = cms.InputTag("packedPFCandidates")
)


process.pfPhotons = cms.EDFilter("CandPtrSelector",
    cut = cms.string('abs(pdgId) = 22'),
    src = cms.InputTag("pfNoPileUp")
)


process.pfPhotonsPuppi = cms.EDFilter("CandPtrSelector",
    cut = cms.string('abs(pdgId) = 22'),
    src = cms.InputTag("pfNoPileUpPuppi")
)


process.pfPileUpAllChargedParticlesPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(
        211, -211, 321, -321, 999211, 
        2212, -2212, 11, -11, 13, 
        -13
    ),
    src = cms.InputTag("pfPileUpIsoPFBRECO")
)


process.pfTaus = cms.EDFilter("PATTauRefSelector",
    cut = cms.string('pt > 18.0 & abs(eta) < 2.6 & tauID("decayModeFinding") > 0.5 & isPFTau'),
    src = cms.InputTag("slimmedTaus")
)


process.pfTausPuppi = cms.EDFilter("PATTauRefSelector",
    cut = cms.string('pt > 18.0 & abs(eta) < 2.6 & tauID("decayModeFinding") > 0.5 & isPFTau'),
    src = cms.InputTag("slimmedTaus")
)


process.pfTrk = cms.EDFilter("CandPtrSelector",
    cut = cms.string('fromPV(0) > 0 && charge()!=0'),
    src = cms.InputTag("packedPFCandidates")
)


process.selectedPatElectrons = cms.EDFilter("PATElectronSelector",
    cut = cms.string(''),
    src = cms.InputTag("patElectrons")
)


process.selectedPatJets = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2CorrPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPatJetsForMetT1T2SmearCorrPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJetsPuppi")
)


process.selectedPatJetsForMetT2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT2CorrPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPatJetsForMetT2SmearCorrPuppi = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJetsPuppi")
)


process.selectedPatLowPtElectrons = cms.EDFilter("PATElectronSelector",
    cut = cms.string(''),
    src = cms.InputTag("patLowPtElectrons")
)


process.selectedPatMuons = cms.EDFilter("PATMuonSelector",
    cut = cms.string(''),
    src = cms.InputTag("patMuons")
)


process.selectedPatOOTPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string(''),
    src = cms.InputTag("patOOTPhotons")
)


process.selectedPatPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string(''),
    src = cms.InputTag("patPhotons")
)


process.selectedPatTaus = cms.EDFilter("PATTauSelector",
    cut = cms.string(''),
    src = cms.InputTag("patTaus")
)


process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0 = cms.EDFilter("PATSingleVertexSelector",
    filter = cms.bool(False),
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag("selectedVerticesForPFMEtCorrType0")
)


process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0Puppi = cms.EDFilter("PATSingleVertexSelector",
    filter = cms.bool(False),
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag("selectedVerticesForPFMEtCorrType0Puppi")
)


process.selectedUpdatedPatJetsUpdatedJEC = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("updatedPatJetsUpdatedJEC")
)


process.selectedVerticesForPFMEtCorrType0 = cms.EDFilter("VertexSelector",
    cut = cms.string('isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.'),
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices")
)


process.selectedVerticesForPFMEtCorrType0Puppi = cms.EDFilter("VertexSelector",
    cut = cms.string('isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.'),
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices")
)


process.tauGenJetsSelectorAllHadrons = cms.EDFilter("TauGenJetDecayModeSelector",
    filter = cms.bool(False),
    select = cms.vstring(
        'oneProng0Pi0', 
        'oneProng1Pi0', 
        'oneProng2Pi0', 
        'oneProngOther', 
        'threeProng0Pi0', 
        'threeProng1Pi0', 
        'threeProngOther', 
        'rare'
    ),
    src = cms.InputTag("tauGenJets")
)


process.cleanPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("cleanPatElectrons"), cms.InputTag("cleanPatMuons"), cms.InputTag("cleanPatTaus"), cms.InputTag("cleanPatPhotons"), cms.InputTag("cleanPatJets")),
    logName = cms.untracked.string('cleanPatCandidates|PATSummaryTables')
)


process.content = cms.EDAnalyzer("EventContentAnalyzer")


process.demo = cms.EDAnalyzer("TM",
    HLTriggerResults_ = cms.untracked.InputTag("TriggerResults","","HLT"),
    debug_ = cms.untracked.bool(True),
    ecalBadCalibLabel_ = cms.untracked.InputTag("ecalBadCalibReducedMINIAODFilter"),
    electronLabel_ = cms.untracked.InputTag("slimmedElectrons"),
    fillPhotInfo_ = cms.untracked.bool(False),
    fillelectronInfo_ = cms.untracked.bool(True),
    filleventInfo_ = cms.untracked.bool(True),
    fillgenparticleInfo_ = cms.untracked.bool(True),
    fillmuonInfo_ = cms.untracked.bool(True),
    fillpfjetInfo_ = cms.untracked.bool(True),
    fillpfmetInfo_ = cms.untracked.bool(True),
    fillpileUpInfo_ = cms.untracked.bool(True),
    fillrhoInfo_ = cms.untracked.bool(True),
    filltauInfo_ = cms.untracked.bool(True),
    filltriggerInfo_ = cms.untracked.bool(True),
    fillvertexInfo_ = cms.untracked.bool(True),
    genEventLabel_ = cms.untracked.InputTag("generator"),
    genParticleLabel_ = cms.untracked.InputTag("prunedGenParticles"),
    generatorlheLabel_ = cms.untracked.InputTag("source"),
    hltlabel_ = cms.untracked.string('HLT'),
    isMC_ = cms.untracked.bool(True),
    metfilterspatLabel_ = cms.untracked.InputTag("TriggerResults","","PAT"),
    metfiltersrecoLabel_ = cms.untracked.InputTag("TriggerResults","","RECO"),
    muonLabel_ = cms.untracked.InputTag("slimmedMuons"),
    pfjetLabel_ = cms.untracked.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
    pfmetLabel_ = cms.untracked.InputTag("slimmedMETs"),
    photonLabel_ = cms.untracked.InputTag("slimmedPhotons"),
    pileUpLabel_ = cms.untracked.InputTag("slimmedAddPileupInfo"),
    puppimetLabel_ = cms.untracked.InputTag("slimmedMETsPuppi"),
    rhoLabel_ = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
    tauLabel_ = cms.untracked.InputTag("slimmedTausNewID"),
    triggerEventTag_ = cms.untracked.InputTag("hltTriggerSummaryAOD","","HLT"),
    vertexLabel_ = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
)


process.patCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(
        cms.InputTag("patElectrons"), cms.InputTag("patLowPtElectrons"), cms.InputTag("patMuons"), cms.InputTag("patTaus"), cms.InputTag("patPhotons"), 
        cms.InputTag("patOOTPhotons"), cms.InputTag("patJets"), cms.InputTag("patMETs")
    ),
    logName = cms.untracked.string('patCandidates|PATSummaryTables')
)


process.selectedPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(
        cms.InputTag("selectedPatElectrons"), cms.InputTag("selectedPatLowPtElectrons"), cms.InputTag("selectedPatMuons"), cms.InputTag("selectedPatTaus"), cms.InputTag("selectedPatPhotons"), 
        cms.InputTag("selectedPatOOTPhotons"), cms.InputTag("selectedPatJets")
    ),
    logName = cms.untracked.string('selectedPatCanddiates|PATSummaryTables')
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring(
        'FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'
    ),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring(
        'warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'
    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    closeFileFast = cms.untracked.bool(True),
    fileName = cms.string('histo.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
    MapAuto = cms.untracked.bool(False),
    MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz'),
    SkipHE = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.ak10PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL1Fastjet', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL1Offset', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute'
    )
)


process.ak10PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual'
    )
)


process.ak10PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2Relative')
)


process.ak10PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L3Absolute')
)


process.ak10PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak10PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL1Fastjet', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual'
    )
)


process.ak10PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL1Offset', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual'
    )
)


process.ak10PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL2Relative', 
        'ak10PFL3Absolute'
    )
)


process.ak10PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual'
    )
)


process.ak10PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2Relative')
)


process.ak10PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L3Absolute')
)


process.ak10PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2L3Residual')
)


process.ak1PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL1Fastjet', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL1Offset', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute'
    )
)


process.ak1PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual'
    )
)


process.ak1PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2Relative')
)


process.ak1PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L3Absolute')
)


process.ak1PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak1PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL1Fastjet', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual'
    )
)


process.ak1PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL1Offset', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual'
    )
)


process.ak1PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL2Relative', 
        'ak1PFL3Absolute'
    )
)


process.ak1PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual'
    )
)


process.ak1PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2Relative')
)


process.ak1PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L3Absolute')
)


process.ak1PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2L3Residual')
)


process.ak2PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL1Fastjet', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL1Offset', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute'
    )
)


process.ak2PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual'
    )
)


process.ak2PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2Relative')
)


process.ak2PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L3Absolute')
)


process.ak2PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak2PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL1Fastjet', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual'
    )
)


process.ak2PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL1Offset', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual'
    )
)


process.ak2PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL2Relative', 
        'ak2PFL3Absolute'
    )
)


process.ak2PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual'
    )
)


process.ak2PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2Relative')
)


process.ak2PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L3Absolute')
)


process.ak2PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2L3Residual')
)


process.ak3PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL1Fastjet', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL1Offset', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute'
    )
)


process.ak3PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual'
    )
)


process.ak3PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2Relative')
)


process.ak3PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L3Absolute')
)


process.ak3PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak3PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL1Fastjet', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual'
    )
)


process.ak3PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL1Offset', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual'
    )
)


process.ak3PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL2Relative', 
        'ak3PFL3Absolute'
    )
)


process.ak3PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual'
    )
)


process.ak3PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2Relative')
)


process.ak3PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L3Absolute')
)


process.ak3PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2L3Residual')
)


process.ak4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB'
    )
)


process.ak4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual'
    )
)


process.ak4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual'
    )
)


process.ak4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute'
    )
)


process.ak4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB'
    )
)


process.ak4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual'
    )
)


process.ak4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTFastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTFastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual'
    )
)


process.ak4JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual'
    )
)


process.ak4JPTL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute'
    )
)


process.ak4JPTL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual'
    )
)


process.ak4JPTL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Fastjet')
)


process.ak4L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Offset')
)


process.ak4PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute'
    )
)


process.ak4PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual'
    )
)


process.ak4PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute'
    )
)


process.ak4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB'
    )
)


process.ak4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual'
    )
)


process.ak4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute'
    )
)


process.ak4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual'
    )
)


process.ak4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative', 
        'ak4PFL3Absolute'
    )
)


process.ak4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB'
    )
)


process.ak4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual'
    )
)


process.ak4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ak4TrackL2Relative', 
        'ak4TrackL3Absolute'
    )
)


process.ak4TrackL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4TrackL2Relative', 
        'ak4TrackL3Absolute'
    )
)


process.ak4TrackL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL1Fastjet', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL1Offset', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute'
    )
)


process.ak5PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual'
    )
)


process.ak5PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2Relative')
)


process.ak5PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual'
    )
)


process.ak5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL1Offset', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual'
    )
)


process.ak5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL2Relative', 
        'ak5PFL3Absolute'
    )
)


process.ak5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual'
    )
)


process.ak5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)


process.ak5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)


process.ak5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.ak6PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL1Fastjet', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL1Offset', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute'
    )
)


process.ak6PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual'
    )
)


process.ak6PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2Relative')
)


process.ak6PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L3Absolute')
)


process.ak6PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL1Fastjet', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual'
    )
)


process.ak6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL1Offset', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual'
    )
)


process.ak6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL2Relative', 
        'ak6PFL3Absolute'
    )
)


process.ak6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual'
    )
)


process.ak6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2Relative')
)


process.ak6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L3Absolute')
)


process.ak6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2L3Residual')
)


process.ak7CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB'
    )
)


process.ak7CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual'
    )
)


process.ak7CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak7CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual'
    )
)


process.ak7CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute'
    )
)


process.ak7CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB'
    )
)


process.ak7CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual'
    )
)


process.ak7CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2Relative')
)


process.ak7CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7CaloJetsSoftMuonTagInfos")
)


process.ak7CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2L3Residual')
)


process.ak7JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTFastjet', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute'
    )
)


process.ak7JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTFastjet', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual'
    )
)


process.ak7JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute'
    )
)


process.ak7JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual'
    )
)


process.ak7JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute'
    )
)


process.ak7L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Fastjet')
)


process.ak7L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Offset')
)


process.ak7PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute'
    )
)


process.ak7PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL1Offset', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute'
    )
)


process.ak7PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual'
    )
)


process.ak7PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2Relative')
)


process.ak7PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L3Absolute')
)


process.ak7PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak7PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute'
    )
)


process.ak7PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB'
    )
)


process.ak7PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual'
    )
)


process.ak7PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute'
    )
)


process.ak7PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual'
    )
)


process.ak7PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative', 
        'ak7PFL3Absolute'
    )
)


process.ak7PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB'
    )
)


process.ak7PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual'
    )
)


process.ak7PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2Relative')
)


process.ak7PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L3Absolute')
)


process.ak7PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7PFJetsSoftMuonTagInfos")
)


process.ak7PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2L3Residual')
)


process.ak8PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL1Fastjet', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL1Offset', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute'
    )
)


process.ak8PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual'
    )
)


process.ak8PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2Relative')
)


process.ak8PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L3Absolute')
)


process.ak8PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak8PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL1Fastjet', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual'
    )
)


process.ak8PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL1Offset', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual'
    )
)


process.ak8PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL2Relative', 
        'ak8PFL3Absolute'
    )
)


process.ak8PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual'
    )
)


process.ak8PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2Relative')
)


process.ak8PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L3Absolute')
)


process.ak8PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2L3Residual')
)


process.ak9PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL1Fastjet', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL1Offset', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute'
    )
)


process.ak9PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual'
    )
)


process.ak9PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2Relative')
)


process.ak9PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L3Absolute')
)


process.ak9PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak9PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL1Fastjet', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual'
    )
)


process.ak9PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL1Offset', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual'
    )
)


process.ak9PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL2Relative', 
        'ak9PFL3Absolute'
    )
)


process.ak9PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual'
    )
)


process.ak9PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2Relative')
)


process.ak9PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L3Absolute')
)


process.ak9PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2L3Residual')
)


process.ecalSeverityLevel = cms.ESProducer("EcalSeverityLevelESProducer",
    dbstatusMask = cms.PSet(
        kBad = cms.vstring(
            'kNonRespondingIsolated', 
            'kDeadVFE', 
            'kDeadFE', 
            'kNoDataNoTP'
        ),
        kGood = cms.vstring('kOk'),
        kProblematic = cms.vstring(
            'kDAC', 
            'kNoLaser', 
            'kNoisy', 
            'kNNoisy', 
            'kNNNoisy', 
            'kNNNNoisy', 
            'kNNNNNoisy', 
            'kFixedG6', 
            'kFixedG1', 
            'kFixedG0'
        ),
        kRecovered = cms.vstring(),
        kTime = cms.vstring(),
        kWeird = cms.vstring()
    ),
    flagMask = cms.PSet(
        kBad = cms.vstring(
            'kFaultyHardware', 
            'kDead', 
            'kKilled'
        ),
        kGood = cms.vstring('kGood'),
        kProblematic = cms.vstring(
            'kPoorReco', 
            'kPoorCalib', 
            'kNoisy', 
            'kSaturated'
        ),
        kRecovered = cms.vstring(
            'kLeadingEdgeRecovered', 
            'kTowerRecovered'
        ),
        kTime = cms.vstring('kOutOfTime'),
        kWeird = cms.vstring(
            'kWeird', 
            'kDiWeird'
        )
    ),
    timeThresh = cms.double(2.0)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.ic5CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB'
    )
)


process.ic5CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual'
    )
)


process.ic5CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ic5CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual'
    )
)


process.ic5CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute'
    )
)


process.ic5CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB'
    )
)


process.ic5CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual'
    )
)


process.ic5CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2Relative')
)


process.ic5CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L3Absolute')
)


process.ic5CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5CaloJetsSoftMuonTagInfos")
)


process.ic5CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2L3Residual')
)


process.ic5PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute'
    )
)


process.ic5PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB'
    )
)


process.ic5PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual'
    )
)


process.ic5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ic5PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute'
    )
)


process.ic5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual'
    )
)


process.ic5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative', 
        'ic5PFL3Absolute'
    )
)


process.ic5PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB'
    )
)


process.ic5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual'
    )
)


process.ic5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2Relative')
)


process.ic5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L3Absolute')
)


process.ic5PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5PFJetsSoftMuonTagInfos")
)


process.ic5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2L3Residual')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.kt4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB'
    )
)


process.kt4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual'
    )
)


process.kt4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual'
    )
)


process.kt4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute'
    )
)


process.kt4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB'
    )
)


process.kt4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual'
    )
)


process.kt4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2Relative')
)


process.kt4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L3Absolute')
)


process.kt4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4CaloJetsSoftMuonTagInfos")
)


process.kt4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2L3Residual')
)


process.kt4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute'
    )
)


process.kt4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB'
    )
)


process.kt4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual'
    )
)


process.kt4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute'
    )
)


process.kt4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual'
    )
)


process.kt4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative', 
        'kt4PFL3Absolute'
    )
)


process.kt4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB'
    )
)


process.kt4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual'
    )
)


process.kt4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2Relative')
)


process.kt4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L3Absolute')
)


process.kt4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4PFJetsSoftMuonTagInfos")
)


process.kt4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2L3Residual')
)


process.kt6CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB'
    )
)


process.kt6CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual'
    )
)


process.kt6CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt6CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual'
    )
)


process.kt6CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute'
    )
)


process.kt6CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB'
    )
)


process.kt6CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual'
    )
)


process.kt6CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2Relative')
)


process.kt6CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L3Absolute')
)


process.kt6CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6CaloJetsSoftMuonTagInfos")
)


process.kt6CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2L3Residual')
)


process.kt6PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute'
    )
)


process.kt6PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB'
    )
)


process.kt6PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual'
    )
)


process.kt6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt6PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute'
    )
)


process.kt6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual'
    )
)


process.kt6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative', 
        'kt6PFL3Absolute'
    )
)


process.kt6PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB'
    )
)


process.kt6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual'
    )
)


process.kt6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2Relative')
)


process.kt6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L3Absolute')
)


process.kt6PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6PFJetsSoftMuonTagInfos")
)


process.kt6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2L3Residual')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    ),
    siPixelQualityLabel = cms.string('')
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('102X_upgrade2018_realistic_v18'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ), 
        cms.PSet(
            slope = cms.double(-1.5610227),
            tmax = cms.double(10.0),
            tzero = cms.double(11.977461)
        ), 
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.essourceEcalSev = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalSeverityLevelAlgoRcd')
)


process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('sqlite:Autumn18_V19_MC.db'),
    toGet = cms.VPSet(
        cms.PSet(
            label = cms.untracked.string('AK4PFchs'),
            record = cms.string('JetCorrectionsRecord'),
            tag = cms.string('JetCorrectorParametersCollection_Autumn18_V19_MC_AK4PFchs')
        ), 
        cms.PSet(
            label = cms.untracked.string('AK4PFPuppi'),
            record = cms.string('JetCorrectionsRecord'),
            tag = cms.string('JetCorrectorParametersCollection_Autumn18_V19_MC_AK4PFPuppi')
        )
    )
)


process.loadRecoTauTagMVAsFromPrepDB = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string(''),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet( (
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwoLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwoLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAoldDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAoldDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAnewDMwoLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAnewDMwoLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBoldDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBoldDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBnewDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBnewDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWoldDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWnewDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAPWdR03oldDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVADBdR03oldDMwLTv1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2016v1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2016v1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff95'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff95')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff95'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff95')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff50'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff50')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff95'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff95')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff40'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff40')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_mvaOutput_normalization')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_wGwoGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_woGwoGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_wGwoGSF_BL_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA5v1_gbr_NoEleMatch_woGwoGSF_EC_WPeff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff99'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff99')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff96'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff96')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff91'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff91')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff85'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff85')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff79'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC_WPEff79')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_EC_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_BL_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_EC_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_BL_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_woGwGSF_BL_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_wGwGSF_BL_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_wGwoGSF_EC_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff98'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff98')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff90'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff90')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff80'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff80')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff70'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff70')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff60'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_antiElectronMVA6v3_noeveto_gbr_NoEleMatch_woGwoGSF_EC_WPeff60')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_againstMuonMVAv1'),
            record = cms.string('GBRWrapperRcd'),
            tag = cms.string('RecoTauTag_againstMuonMVAv1')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_againstMuonMVAv1_WPeff99_5'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_againstMuonMVAv1_WPeff99_5')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_againstMuonMVAv1_WPeff99_0'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_againstMuonMVAv1_WPeff99_0')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_againstMuonMVAv1_WPeff98_0'),
            record = cms.string('PhysicsTGraphPayloadRcd'),
            tag = cms.string('RecoTauTag_againstMuonMVAv1_WPeff98_0')
        ), 
        cms.PSet(
            label = cms.untracked.string('RecoTauTag_againstMuonMVAv1_mvaOutput_normalization'),
            record = cms.string('PhysicsTFormulaPayloadRcd'),
            tag = cms.string('RecoTauTag_againstMuonMVAv1_mvaOutput_normalization')
        )
     ) )
)


process.prefer("es_hardcode")

process.prefer("jec")

process.ak4CaloL2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL2L3ResidualCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.pfNoPileUpPFBRECOTask = cms.Task(process.pfNoPileUpPFBRECO, process.pfPileUpPFBRECO)


process.ootPhotonTask = cms.Task(process.ootPhotonCore, process.ootPhotons, process.ootPhotonsTmp, process.particleFlowClusterOOTECAL, process.particleFlowClusterOOTECALUncorrected, process.particleFlowRecHitOOTECAL, process.particleFlowSuperClusterOOTECAL)


process.patPFTauIsolationTask = cms.Task(process.tauIsoDepositPFCandidates, process.tauIsoDepositPFChargedHadrons, process.tauIsoDepositPFGammas, process.tauIsoDepositPFNeutralHadrons)


process.photonPFIsolationValuesPATTask = cms.Task(process.phPFIsoValueCharged03PFIdPAT, process.phPFIsoValueCharged04PFIdPAT, process.phPFIsoValueChargedAll03PFIdPAT, process.phPFIsoValueChargedAll04PFIdPAT, process.phPFIsoValueGamma03PFIdPAT, process.phPFIsoValueGamma04PFIdPAT, process.phPFIsoValueNeutral03PFIdPAT, process.phPFIsoValueNeutral04PFIdPAT, process.phPFIsoValuePU03PFIdPAT, process.phPFIsoValuePU04PFIdPAT)


process.ak4PFL1L2L3CorrectorTask = cms.Task(process.ak4PFL1L2L3Corrector, process.ak4PFL1OffsetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.patJetCorrectionsTask = cms.Task(process.patJetCorrFactors)


process.egammaPostRecoPatUpdatorTask = cms.Task(process.slimmedElectrons, process.slimmedPhotons)


process.updateHPSPFTausTask = cms.Task(process.hpsPFTauChargedIsoPtSum, process.hpsPFTauDiscriminationByLoosePileupWeightedIsolation3Hits, process.hpsPFTauDiscriminationByMediumPileupWeightedIsolation3Hits, process.hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone, process.hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits, process.hpsPFTauDiscriminationByRawPileupWeightedIsolation3Hits, process.hpsPFTauDiscriminationByTightPileupWeightedIsolation3Hits, process.hpsPFTauFootprintCorrection, process.hpsPFTauNeutralIsoPtSum, process.hpsPFTauNeutralIsoPtSumWeight, process.hpsPFTauPUcorrPtSum, process.hpsPFTauPhotonPtSumOutsideSignalCone)


process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL1FastL2L3ResidualCorrector, process.ak4PFPuppiL1FastjetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4PFL1FastL2L3CorrectorTask = cms.Task(process.ak4PFL1FastL2L3Corrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.ak4TrackL2L3CorrectorTask = cms.Task(process.ak4TrackL2L3Corrector, process.ak4TrackL2RelativeCorrector, process.ak4TrackL3AbsoluteCorrector)


process.makePatOOTPhotonsTask = cms.Task(process.ootPhotonMatch, process.patOOTPhotons)


process.electronPFIsolationDepositsPATTask = cms.Task(process.elPFIsoDepositChargedAllPAT, process.elPFIsoDepositChargedPAT, process.elPFIsoDepositGammaPAT, process.elPFIsoDepositNeutralPAT, process.elPFIsoDepositPUPAT)


process.ak4PFPuppiL1L2L3CorrectorTask = cms.Task(process.ak4PFPuppiL1L2L3Corrector, process.ak4PFPuppiL1OffsetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.pfNoPileUpIsoPFBRECOTask = cms.Task(process.pfNoPileUpIsoPFBRECO, process.pfPileUpIsoPFBRECO)


process.patPFMetT2CorrTask = cms.Task(process.patPFMetT2Corr)


process.pfDeepCSVTask = cms.Task(process.pfDeepCMVAJetTags, process.pfDeepCMVATagInfos, process.pfDeepCSVJetTags, process.pfDeepCSVTagInfos)


process.ak4CaloL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL1FastL2L3ResidualCorrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.ak4PFCHSL1FastL2L3CorrectorTask = cms.Task(process.ak4PFCHSL1FastL2L3Corrector, process.ak4PFCHSL1FastjetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4PFL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFL1L2L3ResidualCorrector, process.ak4PFL1OffsetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.patJetFlavourIdLegacyTask = cms.Task(process.patJetFlavourAssociationLegacy, process.patJetPartonAssociationLegacy, process.patJetPartonsLegacy)


process.ak4PFCHSL2L3CorrectorTask = cms.Task(process.ak4PFCHSL2L3Corrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4CaloL2L3CorrectorTask = cms.Task(process.ak4CaloL2L3Corrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.pfSortByTypePFBRECOTask = cms.Task(process.pfAllChargedHadronsPFBRECO, process.pfAllChargedParticlesPFBRECO, process.pfAllNeutralHadronsAndPhotonsPFBRECO, process.pfAllNeutralHadronsPFBRECO, process.pfAllPhotonsPFBRECO, process.pfPileUpAllChargedParticlesPFBRECO)


process.muonPFIsolationDepositsPATTask = cms.Task(process.muPFIsoDepositChargedAllPAT, process.muPFIsoDepositChargedPAT, process.muPFIsoDepositGammaPAT, process.muPFIsoDepositNeutralPAT, process.muPFIsoDepositPUPAT)


process.ak4PFPuppiL1FastL2L3CorrectorTask = cms.Task(process.ak4PFPuppiL1FastL2L3Corrector, process.ak4PFPuppiL1FastjetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.ak4PFL1FastL2L3L6CorrectorTask = cms.Task(process.ak4PFL1FastL2L3L6Corrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFL6SLBCorrector)


process.ak4PFCHSL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL1L2L3ResidualCorrector, process.ak4PFCHSL1OffsetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.ak4PFPuppiL2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL2L3ResidualCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4PFPuppiL2L3CorrectorTask = cms.Task(process.ak4PFPuppiL2L3Corrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.ak4PFCHSL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL1FastL2L3ResidualCorrector, process.ak4PFCHSL1FastjetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.correctionTermsPfMetType1Type2Task = cms.Task(process.ak4PFCHSL1FastL2L3CorrectorTask, process.ak4PFCHSL1FastL2L3ResidualCorrectorTask, process.corrPfMetType1, process.corrPfMetType2, process.particleFlowPtrs, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.pfCandsNotInJetsPtrForMetCorr, process.pfJetsPtrForMetCorr)


process.patJetFlavourIdTask = cms.Task(process.patJetFlavourAssociation, process.patJetPartons)


process.cleanPatCandidatesTask = cms.Task(process.cleanPatElectrons, process.cleanPatJets, process.cleanPatMuons, process.cleanPatPhotons, process.cleanPatTaus)


process.ak4CaloL1FastL2L3L6CorrectorTask = cms.Task(process.ak4CaloL1FastL2L3L6Corrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloL6SLBCorrector)


process.type0PFMEtCorrectionPFCandToVertexAssociationTask = cms.Task(process.particleFlowDisplacedVertex, process.pfCandidateToVertexAssociation, process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0, process.selectedVerticesForPFMEtCorrType0)


process.ak4PFL2L3L6CorrectorTask = cms.Task(process.ak4PFL2L3L6Corrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFL6SLBCorrector)


process.ak4CaloL1L2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL1L2L3ResidualCorrector, process.ak4CaloL1OffsetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.ak4PFPuppiL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL1L2L3ResidualCorrector, process.ak4PFPuppiL1OffsetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4L1JPTOffsetCorrectorTask = cms.Task(process.ak4CaloL1OffsetCorrector, process.ak4L1JPTOffsetCorrector)


process.egammaVIDTask = cms.Task()


process.makePatLowPtElectronsTask = cms.Task(process.lowPtElectronMatch, process.patLowPtElectrons)


process.ak4JPTL1L2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL1L2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.patPFMetTxyCorrTask = cms.Task(process.patPFMetTxyCorr)


process.producePatPFMETCorrectionsTask = cms.Task(process.patPFMet, process.patPFMetT0Corr, process.patPFMetT0pcT1, process.patPFMetT0pcT1T2, process.patPFMetT1, process.patPFMetT1T2, process.patPFMetT1T2Corr, process.patPFMetT2Corr, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.selectedPatJetsForMetT1T2Corr, process.selectedPatJetsForMetT2Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.egmPhotonIDTask = cms.Task(process.egmPhotonIDs, process.photonMVAValueMapProducer)


process.ak4PFCHSL2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL2L3ResidualCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.ak4CaloL1FastL2L3CorrectorTask = cms.Task(process.ak4CaloL1FastL2L3Corrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.ak4CaloL2L3L6CorrectorTask = cms.Task(process.ak4CaloL2L3L6Corrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloL6SLBCorrector)


process.ak4L1JPTFastjetCorrectorTask = cms.Task(process.ak4CaloL1FastjetCorrector, process.ak4L1JPTFastjetCorrector)


process.patPFMetT2SmearCorrTask = cms.Task(process.patPFMetT1T2SmearCorr, process.patPFMetT2SmearCorr, process.patSmearedJets, process.selectedPatJetsForMetT1T2SmearCorr, process.selectedPatJetsForMetT2SmearCorr)


process.correctionTermsCaloMetTask = cms.Task(process.ak4CaloL2L3CorrectorTask, process.ak4CaloL2L3ResidualCorrectorTask, process.corrCaloMetType1, process.corrCaloMetType2, process.muCaloMetCorr)


process.ak4JPTL2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.ak4JPTL2L3CorrectorTask = cms.Task(process.ak4JPTL2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.muonPFIsolationValuesPATTask = cms.Task(process.muPFIsoValueCharged03PAT, process.muPFIsoValueCharged04PAT, process.muPFIsoValueChargedAll03PAT, process.muPFIsoValueChargedAll04PAT, process.muPFIsoValueGamma03PAT, process.muPFIsoValueGamma04PAT, process.muPFIsoValueGammaHighThreshold03PAT, process.muPFIsoValueGammaHighThreshold04PAT, process.muPFIsoValueNeutral03PAT, process.muPFIsoValueNeutral04PAT, process.muPFIsoValueNeutralHighThreshold03PAT, process.muPFIsoValueNeutralHighThreshold04PAT, process.muPFIsoValuePU03PAT, process.muPFIsoValuePU04PAT, process.muPFMeanDRIsoValueCharged03PAT, process.muPFMeanDRIsoValueCharged04PAT, process.muPFMeanDRIsoValueChargedAll03PAT, process.muPFMeanDRIsoValueChargedAll04PAT, process.muPFMeanDRIsoValueGamma03PAT, process.muPFMeanDRIsoValueGamma04PAT, process.muPFMeanDRIsoValueGammaHighThreshold03PAT, process.muPFMeanDRIsoValueGammaHighThreshold04PAT, process.muPFMeanDRIsoValueNeutral03PAT, process.muPFMeanDRIsoValueNeutral04PAT, process.muPFMeanDRIsoValueNeutralHighThreshold03PAT, process.muPFMeanDRIsoValueNeutralHighThreshold04PAT, process.muPFMeanDRIsoValuePU03PAT, process.muPFMeanDRIsoValuePU04PAT, process.muPFSumDRIsoValueCharged03PAT, process.muPFSumDRIsoValueCharged04PAT, process.muPFSumDRIsoValueChargedAll03PAT, process.muPFSumDRIsoValueChargedAll04PAT, process.muPFSumDRIsoValueGamma03PAT, process.muPFSumDRIsoValueGamma04PAT, process.muPFSumDRIsoValueGammaHighThreshold03PAT, process.muPFSumDRIsoValueGammaHighThreshold04PAT, process.muPFSumDRIsoValueNeutral03PAT, process.muPFSumDRIsoValueNeutral04PAT, process.muPFSumDRIsoValueNeutralHighThreshold03PAT, process.muPFSumDRIsoValueNeutralHighThreshold04PAT, process.muPFSumDRIsoValuePU03PAT, process.muPFSumDRIsoValuePU04PAT)


process.egmPhotonIsolationMiniAODTask = cms.Task(process.egmPhotonIsolation)


process.egammaUpdatorTask = cms.Task(process.egmPhotonIsolationMiniAODTask, process.heepIDVarValueMaps, process.photonIDValueMapProducer, process.updatedElectrons, process.updatedPhotons)


process.rerunMvaIsolationTask = cms.Task(process.deepTau2017v2p1)


process.photonPFIsolationDepositsPATTask = cms.Task(process.phPFIsoDepositChargedAllPAT, process.phPFIsoDepositChargedPAT, process.phPFIsoDepositGammaPAT, process.phPFIsoDepositNeutralPAT, process.phPFIsoDepositPUPAT)


process.patPFMetT0CorrTask = cms.Task(process.patPFMetT0Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFL1FastL2L3ResidualCorrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.egammaScaleSmearTask = cms.Task(process.calibratedPatElectrons, process.calibratedPatPhotons)


process.patPFMetSmearCorrTask = cms.Task(process.patPFMetT1T2SmearCorr, process.patSmearedJets, process.selectedPatJetsForMetT1T2SmearCorr)


process.selectedPatCandidatesTask = cms.Task(process.selectedPatElectrons, process.selectedPatJets, process.selectedPatLowPtElectrons, process.selectedPatMuons, process.selectedPatOOTPhotons, process.selectedPatPhotons, process.selectedPatTaus)


process.electronPFIsolationValuesPATTask = cms.Task(process.elPFIsoValueCharged03NoPFIdPAT, process.elPFIsoValueCharged03PFIdPAT, process.elPFIsoValueCharged04NoPFIdPAT, process.elPFIsoValueCharged04PFIdPAT, process.elPFIsoValueChargedAll03NoPFIdPAT, process.elPFIsoValueChargedAll03PFIdPAT, process.elPFIsoValueChargedAll04NoPFIdPAT, process.elPFIsoValueChargedAll04PFIdPAT, process.elPFIsoValueGamma03NoPFIdPAT, process.elPFIsoValueGamma03PFIdPAT, process.elPFIsoValueGamma04NoPFIdPAT, process.elPFIsoValueGamma04PFIdPAT, process.elPFIsoValueNeutral03NoPFIdPAT, process.elPFIsoValueNeutral03PFIdPAT, process.elPFIsoValueNeutral04NoPFIdPAT, process.elPFIsoValueNeutral04PFIdPAT, process.elPFIsoValuePU03NoPFIdPAT, process.elPFIsoValuePU03PFIdPAT, process.elPFIsoValuePU04NoPFIdPAT, process.elPFIsoValuePU04PFIdPAT)


process.patPFMetT1T2CorrTask = cms.Task(process.patPFMetT1T2Corr, process.selectedPatJetsForMetT1T2Corr)


process.ak4PFCHSL1L2L3CorrectorTask = cms.Task(process.ak4PFCHSL1L2L3Corrector, process.ak4PFCHSL1OffsetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4CaloL1L2L3CorrectorTask = cms.Task(process.ak4CaloL1L2L3Corrector, process.ak4CaloL1OffsetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.ak4PFL2L3ResidualCorrectorTask = cms.Task(process.ak4PFL2L3ResidualCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.ak4PFL2L3CorrectorTask = cms.Task(process.ak4PFL2L3Corrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.pfElectronIsolationPATTask = cms.Task(process.electronPFIsolationDepositsPATTask, process.electronPFIsolationValuesPATTask)


process.pfPhotonIsolationPATTask = cms.Task(process.photonPFIsolationDepositsPATTask, process.photonPFIsolationValuesPATTask)


process.type0PFMEtCorrectionTask = cms.Task(process.pfMETcorrType0, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4JPTL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL1FastL2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTFastjetCorrectorTask)


process.pfParticleSelectionPFBRECOTask = cms.Task(process.pfNoPileUpIsoPFBRECOTask, process.pfNoPileUpPFBRECOTask, process.pfSortByTypePFBRECOTask)


process.patHPSPFTauDiscriminationTask = cms.Task(process.updateHPSPFTausTask)


process.muonPFIsolationPATTask = cms.Task(process.muonPFIsolationDepositsPATTask, process.muonPFIsolationValuesPATTask)


process.patPFCandidateIsoDepositSelectionTask = cms.Task(process.pfNoPileUpIsoPFBRECOTask, process.pfSortByTypePFBRECOTask)


process.ak4JPTL1L2L3CorrectorTask = cms.Task(process.ak4JPTL1L2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.ak4JPTL1FastL2L3CorrectorTask = cms.Task(process.ak4JPTL1FastL2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTFastjetCorrectorTask)


process.patMETCorrectionsTask = cms.Task(process.caloMetT1, process.caloMetT1T2, process.correctionTermsCaloMetTask, process.correctionTermsPfMetType1Type2Task, process.pfMetT1, process.pfMetT1T2)


process.makePatJetsTask = cms.Task(process.patJetCharge, process.patJetCorrectionsTask, process.patJetFlavourIdLegacyTask, process.patJetFlavourIdTask, process.patJetGenJetMatch, process.patJetPartonMatch, process.patJets)


process.producePatPFMETCorrectionsUncTask = cms.Task(process.patPFMet, process.patPFMetT0Corr, process.patPFMetT1T2Corr, process.patPFMetT2Corr, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.selectedPatJetsForMetT1T2Corr, process.selectedPatJetsForMetT2Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.jetCorrectorsTask = cms.Task(process.ak4CaloL1FastL2L3CorrectorTask, process.ak4CaloL1FastL2L3L6CorrectorTask, process.ak4CaloL1FastL2L3ResidualCorrectorTask, process.ak4CaloL1L2L3CorrectorTask, process.ak4CaloL1L2L3ResidualCorrectorTask, process.ak4CaloL2L3CorrectorTask, process.ak4CaloL2L3L6CorrectorTask, process.ak4CaloL2L3ResidualCorrectorTask, process.ak4JPTL1FastL2L3CorrectorTask, process.ak4JPTL1FastL2L3ResidualCorrectorTask, process.ak4JPTL1L2L3CorrectorTask, process.ak4JPTL1L2L3ResidualCorrectorTask, process.ak4JPTL2L3CorrectorTask, process.ak4JPTL2L3ResidualCorrectorTask, process.ak4L1JPTFastjetCorrectorTask, process.ak4L1JPTOffsetCorrectorTask, process.ak4PFCHSL1FastL2L3CorrectorTask, process.ak4PFCHSL1FastL2L3ResidualCorrectorTask, process.ak4PFCHSL1L2L3CorrectorTask, process.ak4PFCHSL1L2L3ResidualCorrectorTask, process.ak4PFCHSL2L3CorrectorTask, process.ak4PFCHSL2L3ResidualCorrectorTask, process.ak4PFL1FastL2L3CorrectorTask, process.ak4PFL1FastL2L3L6CorrectorTask, process.ak4PFL1FastL2L3ResidualCorrectorTask, process.ak4PFL1L2L3CorrectorTask, process.ak4PFL1L2L3ResidualCorrectorTask, process.ak4PFL2L3CorrectorTask, process.ak4PFL2L3L6CorrectorTask, process.ak4PFL2L3ResidualCorrectorTask, process.ak4PFPuppiL1FastL2L3CorrectorTask, process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask, process.ak4PFPuppiL1L2L3CorrectorTask, process.ak4PFPuppiL1L2L3ResidualCorrectorTask, process.ak4PFPuppiL2L3CorrectorTask, process.ak4PFPuppiL2L3ResidualCorrectorTask, process.ak4TrackL2L3CorrectorTask)


process.makePatMETsTask = cms.Task(process.patMETCorrectionsTask, process.patMETs)


process.pfParticleSelectionForIsoTask = cms.Task(process.particleFlowPtrs, process.pfParticleSelectionPFBRECOTask)


process.makePatMuonsTask = cms.Task(process.muonMatch, process.muonPFIsolationPATTask, process.patMuons, process.pfParticleSelectionForIsoTask)


process.makePatTausTask = cms.Task(process.patPFCandidateIsoDepositSelectionTask, process.patPFTauIsolationTask, process.patTaus, process.tauGenJetMatch, process.tauGenJets, process.tauGenJetsSelectorAllHadrons, process.tauMatch)


process.makePatPhotonsTask = cms.Task(process.patPhotons, process.pfParticleSelectionForIsoTask, process.pfPhotonIsolationPATTask, process.photonMatch)


process.patAlgosToolsTask = cms.Task(process.basicJetsForMet, process.basicJetsForMetPuppi, process.cleanedPatJets, process.cleanedPatJetsPuppi, process.corrPfMetType1Puppi, process.egmPhotonIDTask, process.genMetExtractor, process.genMetExtractorPuppi, process.jetCorrectorsTask, process.jetSelectorForMet, process.jetSelectorForMetPuppi, process.metrawCalo, process.metrawCaloPuppi, process.patCHSMet, process.patCaloMet, process.patJetCorrFactorsReapplyJEC, process.patJetCorrFactorsReapplyJECPuppi, process.patJetCorrFactorsUpdatedJEC, process.patJetsReapplyJEC, process.patJetsReapplyJECPuppi, process.patPFMet, process.patPFMetPuppi, process.patPFMetT0Corr, process.patPFMetT0CorrPuppi, process.patPFMetT0pcT1, process.patPFMetT0pcT1Puppi, process.patPFMetT0pcT1T2, process.patPFMetT0pcT1T2Puppi, process.patPFMetT1, process.patPFMetT1ElectronEnDown, process.patPFMetT1ElectronEnDownPuppi, process.patPFMetT1ElectronEnUp, process.patPFMetT1ElectronEnUpPuppi, process.patPFMetT1JetEnDown, process.patPFMetT1JetEnDownPuppi, process.patPFMetT1JetEnUp, process.patPFMetT1JetEnUpPuppi, process.patPFMetT1JetResDown, process.patPFMetT1JetResDownPuppi, process.patPFMetT1JetResUp, process.patPFMetT1JetResUpPuppi, process.patPFMetT1MuonEnDown, process.patPFMetT1MuonEnDownPuppi, process.patPFMetT1MuonEnUp, process.patPFMetT1MuonEnUpPuppi, process.patPFMetT1PhotonEnDown, process.patPFMetT1PhotonEnDownPuppi, process.patPFMetT1PhotonEnUp, process.patPFMetT1PhotonEnUpPuppi, process.patPFMetT1Puppi, process.patPFMetT1Smear, process.patPFMetT1SmearElectronEnDown, process.patPFMetT1SmearElectronEnDownPuppi, process.patPFMetT1SmearElectronEnUp, process.patPFMetT1SmearElectronEnUpPuppi, process.patPFMetT1SmearJetEnDown, process.patPFMetT1SmearJetEnDownPuppi, process.patPFMetT1SmearJetEnUp, process.patPFMetT1SmearJetEnUpPuppi, process.patPFMetT1SmearJetResDown, process.patPFMetT1SmearJetResDownPuppi, process.patPFMetT1SmearJetResUp, process.patPFMetT1SmearJetResUpPuppi, process.patPFMetT1SmearMuonEnDown, process.patPFMetT1SmearMuonEnDownPuppi, process.patPFMetT1SmearMuonEnUp, process.patPFMetT1SmearMuonEnUpPuppi, process.patPFMetT1SmearPhotonEnDown, process.patPFMetT1SmearPhotonEnDownPuppi, process.patPFMetT1SmearPhotonEnUp, process.patPFMetT1SmearPhotonEnUpPuppi, process.patPFMetT1SmearPuppi, process.patPFMetT1SmearTauEnDown, process.patPFMetT1SmearTauEnDownPuppi, process.patPFMetT1SmearTauEnUp, process.patPFMetT1SmearTauEnUpPuppi, process.patPFMetT1SmearUnclusteredEnDown, process.patPFMetT1SmearUnclusteredEnDownPuppi, process.patPFMetT1SmearUnclusteredEnUp, process.patPFMetT1SmearUnclusteredEnUpPuppi, process.patPFMetT1T2, process.patPFMetT1T2Corr, process.patPFMetT1T2CorrPuppi, process.patPFMetT1T2Puppi, process.patPFMetT1T2SmearCorrPuppi, process.patPFMetT1TauEnDown, process.patPFMetT1TauEnDownPuppi, process.patPFMetT1TauEnUp, process.patPFMetT1TauEnUpPuppi, process.patPFMetT1Txy, process.patPFMetT1TxyPuppi, process.patPFMetT1UnclusteredEnDown, process.patPFMetT1UnclusteredEnDownPuppi, process.patPFMetT1UnclusteredEnUp, process.patPFMetT1UnclusteredEnUpPuppi, process.patPFMetT2Corr, process.patPFMetT2CorrPuppi, process.patPFMetT2SmearCorrPuppi, process.patPFMetT2SmearCorrTask, process.patPFMetTxy, process.patPFMetTxyCorrPuppi, process.patPFMetTxyCorrTask, process.patPFMetTxyPuppi, process.patSmearedJetsPuppi, process.patTrkMet, process.pfCHS, process.pfCandMETcorr, process.pfCandMETcorrPuppi, process.pfCandsForUnclusteredUnc, process.pfCandsForUnclusteredUncPuppi, process.pfCandsNoJets, process.pfCandsNoJetsNoEle, process.pfCandsNoJetsNoEleNoMu, process.pfCandsNoJetsNoEleNoMuNoTau, process.pfCandsNoJetsNoEleNoMuNoTauPuppi, process.pfCandsNoJetsNoEleNoMuPuppi, process.pfCandsNoJetsNoElePuppi, process.pfCandsNoJetsPuppi, process.pfCandsNotInJetsForMetCorr, process.pfCandsNotInJetsForMetCorrPuppi, process.pfElectrons, process.pfElectronsPuppi, process.pfLeptonsPUPPET, process.pfMet, process.pfMetCHS, process.pfMetPuppi, process.pfMetTrk, process.pfMuons, process.pfMuonsPuppi, process.pfNoLepPUPPI, process.pfNoPileUp, process.pfNoPileUpPuppi, process.pfPhotons, process.pfPhotonsPuppi, process.pfTaus, process.pfTausPuppi, process.pfTrk, process.puppi, process.puppiForMET, process.puppiMerged, process.puppiNoLep, process.puppiPhoton, process.selectedPatJetsForMetT1T2Corr, process.selectedPatJetsForMetT1T2CorrPuppi, process.selectedPatJetsForMetT1T2SmearCorrPuppi, process.selectedPatJetsForMetT2Corr, process.selectedPatJetsForMetT2CorrPuppi, process.selectedPatJetsForMetT2SmearCorrPuppi, process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0Puppi, process.selectedUpdatedPatJetsUpdatedJEC, process.selectedVerticesForPFMEtCorrType0Puppi, process.shiftedPatElectronEnDown, process.shiftedPatElectronEnDownPuppi, process.shiftedPatElectronEnUp, process.shiftedPatElectronEnUpPuppi, process.shiftedPatJetEnDown, process.shiftedPatJetEnDownPuppi, process.shiftedPatJetEnUp, process.shiftedPatJetEnUpPuppi, process.shiftedPatJetResDown, process.shiftedPatJetResDownPuppi, process.shiftedPatJetResUp, process.shiftedPatJetResUpPuppi, process.shiftedPatMETCorrElectronEnDown, process.shiftedPatMETCorrElectronEnDownPuppi, process.shiftedPatMETCorrElectronEnUp, process.shiftedPatMETCorrElectronEnUpPuppi, process.shiftedPatMETCorrJetEnDown, process.shiftedPatMETCorrJetEnDownPuppi, process.shiftedPatMETCorrJetEnUp, process.shiftedPatMETCorrJetEnUpPuppi, process.shiftedPatMETCorrJetResDown, process.shiftedPatMETCorrJetResDownPuppi, process.shiftedPatMETCorrJetResUp, process.shiftedPatMETCorrJetResUpPuppi, process.shiftedPatMETCorrMuonEnDown, process.shiftedPatMETCorrMuonEnDownPuppi, process.shiftedPatMETCorrMuonEnUp, process.shiftedPatMETCorrMuonEnUpPuppi, process.shiftedPatMETCorrPhotonEnDown, process.shiftedPatMETCorrPhotonEnDownPuppi, process.shiftedPatMETCorrPhotonEnUp, process.shiftedPatMETCorrPhotonEnUpPuppi, process.shiftedPatMETCorrSmearedJetResDown, process.shiftedPatMETCorrSmearedJetResDownPuppi, process.shiftedPatMETCorrSmearedJetResUp, process.shiftedPatMETCorrSmearedJetResUpPuppi, process.shiftedPatMETCorrTauEnDown, process.shiftedPatMETCorrTauEnDownPuppi, process.shiftedPatMETCorrTauEnUp, process.shiftedPatMETCorrTauEnUpPuppi, process.shiftedPatMETCorrUnclusteredEnDown, process.shiftedPatMETCorrUnclusteredEnDownPuppi, process.shiftedPatMETCorrUnclusteredEnUp, process.shiftedPatMETCorrUnclusteredEnUpPuppi, process.shiftedPatMuonEnDown, process.shiftedPatMuonEnDownPuppi, process.shiftedPatMuonEnUp, process.shiftedPatMuonEnUpPuppi, process.shiftedPatPhotonEnDown, process.shiftedPatPhotonEnDownPuppi, process.shiftedPatPhotonEnUp, process.shiftedPatPhotonEnUpPuppi, process.shiftedPatSmearedJetResDown, process.shiftedPatSmearedJetResDownPuppi, process.shiftedPatSmearedJetResUp, process.shiftedPatSmearedJetResUpPuppi, process.shiftedPatTauEnDown, process.shiftedPatTauEnDownPuppi, process.shiftedPatTauEnUp, process.shiftedPatTauEnUpPuppi, process.shiftedPatUnclusteredEnDown, process.shiftedPatUnclusteredEnDownPuppi, process.shiftedPatUnclusteredEnUp, process.shiftedPatUnclusteredEnUpPuppi, process.slimmedMETs, process.slimmedMETsPuppi, process.type0PFMEtCorrectionPFCandToVertexAssociationTask, process.updatedPatJetsUpdatedJEC)


process.makePatElectronsTask = cms.Task(process.electronMatch, process.patElectrons, process.pfElectronIsolationPATTask, process.pfParticleSelectionForIsoTask)


process.patCandidatesTask = cms.Task(process.makePatElectronsTask, process.makePatJetsTask, process.makePatLowPtElectronsTask, process.makePatMETsTask, process.makePatMuonsTask, process.makePatOOTPhotonsTask, process.makePatPhotonsTask, process.makePatTausTask)


process.pfDeepCSV = cms.Sequence(process.pfDeepCSVTask)


process.patPFMetT2CorrSequencePuppi = cms.Sequence(process.patPFMetT2CorrPuppi)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidation = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4L1JPTOffsetCorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorTask)


process.makePatPhotons = cms.Sequence(process.makePatPhotonsTask)


process.patPFMetT1T2CorrSequence = cms.Sequence(cms.Task(process.patPFMetT1T2Corr))


process.correctionTermsPfMetType1Type2 = cms.Sequence(process.correctionTermsPfMetType1Type2Task)


process.ak4L1JPTFastjetCorrectorChain = cms.Sequence(process.ak4L1JPTFastjetCorrectorTask)


process.patPFMetSmearCorrSequence = cms.Sequence(process.patPFMetSmearCorrTask)


process.ak4PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3ResidualCorrectorTask)


process.patJetFlavourId = cms.Sequence(process.patJetFlavourIdTask)


process.ak4PFCHSL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorTask)


process.ak4CaloL1L2L3CorrectorChain = cms.Sequence(process.ak4CaloL1L2L3CorrectorTask)


process.patPFMetTxyCorrSequencePuppi = cms.Sequence(process.patPFMetTxyCorrPuppi)


process.patJetCorrections = cms.Sequence(process.patJetCorrectionsTask)


process.ak4PFL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3CorrectorTask)


process.makePatElectrons = cms.Sequence(process.makePatElectronsTask)


process.egammaScaleSmearSeq = cms.Sequence(process.egammaScaleSmearTask)


process.patMetModuleSequence = cms.Sequence(process.pfMet+process.genMetExtractor+process.patJetCorrFactorsReapplyJEC+process.patJetsReapplyJEC+process.basicJetsForMet+process.jetSelectorForMet+process.cleanedPatJets+process.metrawCalo+process.pfCHS+process.pfMetCHS+process.patCHSMet+process.pfTrk+process.pfMetTrk+process.patTrkMet+process.patPFMet)


process.egammaPostRecoPatUpdatorSeq = cms.Sequence(process.egammaPostRecoPatUpdatorTask)


process.pfSortByTypePFBRECOSequence = cms.Sequence(process.pfSortByTypePFBRECOTask)


process.electronPFIsolationValuesPATSequence = cms.Sequence(process.electronPFIsolationValuesPATTask)


process.cleanPatCandidates = cms.Sequence(process.cleanPatCandidateSummary, process.cleanPatCandidatesTask)


process.rerunMvaIsolationSequence = cms.Sequence(process.deepTau2017v2p1)


process.patPFMetT1SmearpatShiftedModuleSequence = cms.Sequence(process.patPFMetT1SmearJetResDown+process.patPFMetT1SmearJetResUp+process.patPFMetT1SmearMuonEnUp+process.patPFMetT1SmearMuonEnDown+process.patPFMetT1SmearJetEnUp+process.patPFMetT1SmearJetEnDown+process.patPFMetT1SmearTauEnUp+process.patPFMetT1SmearTauEnDown+process.patPFMetT1SmearPhotonEnUp+process.patPFMetT1SmearPhotonEnDown+process.patPFMetT1SmearElectronEnDown+process.patPFMetT1SmearElectronEnUp+process.patPFMetT1SmearUnclusteredEnUp+process.patPFMetT1SmearUnclusteredEnDown)


process.patHPSPFTauDiscrimination = cms.Sequence(process.patHPSPFTauDiscriminationTask)


process.ak4PFL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3ResidualCorrectorTask)


process.correctionTermsCaloMet = cms.Sequence(process.correctionTermsCaloMetTask)


process.ak4JPTL1L2L3CorrectorChain = cms.Sequence(process.ak4JPTL1L2L3CorrectorTask)


process.patPFMetSmearCorrSequencePuppi = cms.Sequence(process.patSmearedJetsPuppi+process.selectedPatJetsForMetT1T2SmearCorrPuppi+process.patPFMetT1T2SmearCorrPuppi)


process.muonPFIsolationValuesPATSequence = cms.Sequence(process.muonPFIsolationValuesPATTask)


process.patPFMetT2SmearCorrSequencePuppi = cms.Sequence(process.patSmearedJetsPuppi+process.selectedPatJetsForMetT1T2SmearCorrPuppi+process.selectedPatJetsForMetT2SmearCorrPuppi+process.patPFMetT1T2SmearCorrPuppi+process.patPFMetT2SmearCorrPuppi)


process.patJetFlavourIdLegacy = cms.Sequence(process.patJetFlavourIdLegacyTask)


process.pfNoPileUpIsoPFBRECOSequence = cms.Sequence(process.pfNoPileUpIsoPFBRECOTask)


process.ak4PFCHSL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL2L3ResidualCorrectorTask)


process.patPFCandidateIsoDepositSelection = cms.Sequence(process.patPFCandidateIsoDepositSelectionTask)


process.patPFMetT1SmearPuppipatShiftedModuleSequencePuppi = cms.Sequence(process.patPFMetT1SmearJetResDownPuppi+process.patPFMetT1SmearJetResUpPuppi+process.patPFMetT1SmearMuonEnUpPuppi+process.patPFMetT1SmearMuonEnDownPuppi+process.patPFMetT1SmearJetEnUpPuppi+process.patPFMetT1SmearJetEnDownPuppi+process.patPFMetT1SmearTauEnDownPuppi+process.patPFMetT1SmearTauEnUpPuppi+process.patPFMetT1SmearPhotonEnUpPuppi+process.patPFMetT1SmearPhotonEnDownPuppi+process.patPFMetT1SmearElectronEnDownPuppi+process.patPFMetT1SmearElectronEnUpPuppi+process.patPFMetT1SmearUnclusteredEnUpPuppi+process.patPFMetT1SmearUnclusteredEnDownPuppi)


process.ak4CaloL2L3CorrectorChain = cms.Sequence(process.ak4CaloL2L3CorrectorTask)


process.patCandidates = cms.Sequence(process.patCandidateSummary, process.patCandidatesTask)


process.updateHPSPFTaus = cms.Sequence(process.updateHPSPFTausTask)


process.patCaloTauDiscrimination = cms.Sequence()


process.ak4JPTL1FastL2L3CorrectorChain = cms.Sequence(process.ak4JPTL1FastL2L3CorrectorTask)


process.countPatCandidates = cms.Sequence(process.countPatElectrons+process.countPatMuons+process.countPatTaus+process.countPatLeptons+process.countPatPhotons+process.countPatJets)


process.ak4CaloL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3L6CorrectorTask)


process.ak4PFCHSL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL2L3CorrectorTask)


process.patPFMetT0CorrSequence = cms.Sequence(process.patPFMetT0CorrTask)


process.ak4PFL1L2L3CorrectorChain = cms.Sequence(process.ak4PFL1L2L3CorrectorTask)


process.selectedPatCandidates = cms.Sequence(process.selectedPatCandidateSummary, process.selectedPatCandidatesTask)


process.ak4CaloL1FastL2L3CorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3CorrectorTask)


process.patPFMetT2CorrSequence = cms.Sequence(process.patPFMetT2CorrTask)


process.producePatPFMETCorrectionsPuppi = cms.Sequence(process.patPFMetPuppi+process.pfCandsNotInJetsForMetCorrPuppi+process.selectedPatJetsForMetT1T2CorrPuppi+process.selectedPatJetsForMetT2CorrPuppi+process.patPFMetT1T2CorrPuppi+process.patPFMetT2CorrPuppi+process.selectedVerticesForPFMEtCorrType0Puppi+process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0Puppi+process.particleFlowDisplacedVertex+process.pfCandidateToVertexAssociation+process.patPFMetT0CorrPuppi+process.pfCandMETcorrPuppi+process.patPFMetT1Puppi+process.patPFMetT1T2Puppi+process.patPFMetT0pcT1Puppi+process.patPFMetT0pcT1T2Puppi)


process.ak4PFL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1L2L3ResidualCorrectorTask)


process.electronPFIsolationDepositsPATSequence = cms.Sequence(process.electronPFIsolationDepositsPATTask)


process.ak4JPTL2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL2L3ResidualCorrectorTask)


process.ak4PFL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3L6CorrectorTask)


process.photonPFIsolationDepositsPATSequence = cms.Sequence(process.photonPFIsolationDepositsPATTask)


process.patPFMetT1SmearpatMetUncertaintySequence = cms.Sequence(process.shiftedPatSmearedJetResDown+process.shiftedPatMETCorrSmearedJetResDown+process.shiftedPatSmearedJetResUp+process.shiftedPatMETCorrSmearedJetResUp+process.pfElectrons+process.pfTaus+process.pfMuons+process.pfNoPileUp+process.pfPhotons)


process.patShrinkingConePFTauDiscrimination = cms.Sequence()


process.ak4PFPuppiL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL2L3ResidualCorrectorTask)


process.ootPhotonSequence = cms.Sequence(process.ootPhotonTask)


process.ak4JPTL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL1FastL2L3ResidualCorrectorTask)


process.ak4PFL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL2L3L6CorrectorTask)


process.ak4PFPuppiL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1L2L3ResidualCorrectorTask)


process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask)


process.patDefaultSequence = cms.Sequence(process.patCandidates+process.selectedPatCandidates+process.cleanPatCandidates+process.countPatCandidates)


process.pfNoPileUpPFBRECOSequence = cms.Sequence(process.pfNoPileUpPFBRECOTask)


process.patPFTauIsolation = cms.Sequence(process.patPFTauIsolationTask)


process.makePatTaus = cms.Sequence(process.makePatTausTask)


process.ak4CaloL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3ResidualCorrectorTask)


process.photonPFIsolationValuesPATSequence = cms.Sequence(process.photonPFIsolationValuesPATTask)


process.ak4TrackL2L3CorrectorChain = cms.Sequence(process.ak4TrackL2L3CorrectorTask)


process.ak4CaloL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL2L3L6CorrectorTask)


process.pfParticleSelectionForIsoSequence = cms.Sequence(process.pfParticleSelectionForIsoTask)


process.ak4PFPuppiL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastL2L3CorrectorTask)


process.patPFMetT2SmearCorrSequence = cms.Sequence(process.patPFMetT2SmearCorrTask)


process.ak4JPTL2L3CorrectorChain = cms.Sequence(process.ak4JPTL2L3CorrectorTask)


process.puppiMETSequence = cms.Sequence(process.puppi+process.pfLeptonsPUPPET+process.pfNoLepPUPPI+process.puppiNoLep+process.puppiMerged+process.puppiForMET)


process.ak4JPTL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL1L2L3ResidualCorrectorTask)


process.makePatOOTPhotons = cms.Sequence(process.makePatOOTPhotonsTask)


process.patPFMetTxyCorrSequence = cms.Sequence(process.patPFMetTxyCorrTask)


process.producePatPFMETCorrectionsUnc = cms.Sequence(process.producePatPFMETCorrectionsUncTask)


process.egammaUpdatorSeq = cms.Sequence(process.egammaUpdatorTask)


process.pfPhotonIsolationPATSequence = cms.Sequence(process.pfPhotonIsolationPATTask)


process.patPFMetT0CorrSequencePuppi = cms.Sequence(process.selectedVerticesForPFMEtCorrType0Puppi+process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0Puppi+process.particleFlowDisplacedVertex+process.pfCandidateToVertexAssociation+process.patPFMetT0CorrPuppi)


process.makePatJets = cms.Sequence(process.makePatJetsTask)


process.pfParticleSelectionPFBRECOSequence = cms.Sequence(process.pfParticleSelectionPFBRECOTask)


process.patPFMetT1SmearPuppipatMetUncertaintySequencePuppi = cms.Sequence(process.shiftedPatSmearedJetResDownPuppi+process.shiftedPatMETCorrSmearedJetResDownPuppi+process.shiftedPatSmearedJetResUpPuppi+process.shiftedPatMETCorrSmearedJetResUpPuppi+process.pfElectronsPuppi+process.pfTausPuppi+process.pfMuonsPuppi+process.pfNoPileUpPuppi+process.pfPhotonsPuppi)


process.producePatPFMETCorrections = cms.Sequence(process.producePatPFMETCorrectionsTask)


process.ak4PFL2L3CorrectorChain = cms.Sequence(process.ak4PFL2L3CorrectorTask)


process.patPFMetT1T2CorrSequencePuppi = cms.Sequence(process.patPFMetT1T2CorrPuppi)


process.muonPFIsolationDepositsPATSequence = cms.Sequence(process.muonPFIsolationDepositsPATTask)


process.patFixedConePFTauDiscrimination = cms.Sequence()


process.type0PFMEtCorrection = cms.Sequence(process.type0PFMEtCorrectionTask)


process.muonPFIsolationPATSequence = cms.Sequence(process.muonPFIsolationPATTask)


process.makePatLowPtElectrons = cms.Sequence(process.makePatLowPtElectronsTask)


process.egmPhotonIDSequence = cms.Sequence(process.egmPhotonIDTask)


process.pfElectronIsolationPATSequence = cms.Sequence(process.pfElectronIsolationPATTask)


process.egmPhotonIsolationMiniAODSequence = cms.Sequence(process.egmPhotonIsolationMiniAODTask)


process.egammaVIDSeq = cms.Sequence(process.egammaVIDTask)


process.makePatMETs = cms.Sequence(process.makePatMETsTask)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidationMiniAOD = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFCHSL1L2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1L2L3CorrectorTask)


process.ak4PFPuppiL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL2L3CorrectorTask)


process.patMetModuleSequencePuppi = cms.Sequence(process.pfMetPuppi+process.genMetExtractorPuppi+process.patJetCorrFactorsReapplyJECPuppi+process.patJetsReapplyJECPuppi+process.basicJetsForMetPuppi+process.jetSelectorForMetPuppi+process.cleanedPatJetsPuppi+process.metrawCaloPuppi+process.pfCHS+process.pfMetCHS+process.patCHSMet+process.pfTrk+process.pfMetTrk+process.patTrkMet+process.patPFMetPuppi)


process.ak4CaloL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL2L3ResidualCorrectorTask)


process.ak4PFCHSL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1L2L3ResidualCorrectorTask)


process.type0PFMEtCorrectionPFCandToVertexAssociation = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4CaloL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1L2L3ResidualCorrectorTask)


process.makePatMuons = cms.Sequence(process.makePatMuonsTask)


process.ak4PFL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL2L3ResidualCorrectorTask)


process.ak4PFPuppiL1L2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1L2L3CorrectorTask)


process.patMETCorrections = cms.Sequence(process.patMETCorrectionsTask)


process.patMetUncertaintySequencePuppi = cms.Sequence(process.ak4PFPuppiL1FastL2L3CorrectorChain+process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain+process.shiftedPatJetResDownPuppi+process.shiftedPatMETCorrJetResDownPuppi+process.shiftedPatJetResUpPuppi+process.shiftedPatMETCorrJetResUpPuppi+process.pfCandsNoJetsPuppi+process.pfCandsNoJetsNoElePuppi+process.pfCandsNoJetsNoEleNoMuPuppi+process.pfCandsNoJetsNoEleNoMuNoTauPuppi+process.pfCandsForUnclusteredUncPuppi+process.pfElectronsPuppi+process.pfTausPuppi+process.pfMuonsPuppi+process.pfNoPileUpPuppi+process.pfPhotonsPuppi+process.shiftedPatMuonEnDownPuppi+process.shiftedPatMETCorrMuonEnDownPuppi+process.shiftedPatMuonEnUpPuppi+process.shiftedPatMETCorrMuonEnUpPuppi+process.shiftedPatJetEnDownPuppi+process.shiftedPatMETCorrJetEnDownPuppi+process.shiftedPatJetEnUpPuppi+process.shiftedPatMETCorrJetEnUpPuppi+process.shiftedPatTauEnDownPuppi+process.shiftedPatMETCorrTauEnDownPuppi+process.shiftedPatTauEnUpPuppi+process.shiftedPatMETCorrTauEnUpPuppi+process.shiftedPatPhotonEnDownPuppi+process.shiftedPatMETCorrPhotonEnDownPuppi+process.shiftedPatPhotonEnUpPuppi+process.shiftedPatMETCorrPhotonEnUpPuppi+process.shiftedPatElectronEnDownPuppi+process.shiftedPatMETCorrElectronEnDownPuppi+process.shiftedPatElectronEnUpPuppi+process.shiftedPatMETCorrElectronEnUpPuppi+process.shiftedPatUnclusteredEnDownPuppi+process.shiftedPatMETCorrUnclusteredEnDownPuppi+process.shiftedPatUnclusteredEnUpPuppi+process.shiftedPatMETCorrUnclusteredEnUpPuppi+process.patPFMetT1SmearPuppipatMetUncertaintySequencePuppi)


process.patMetUncertaintySequence = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorChain+process.ak4PFCHSL1FastL2L3ResidualCorrectorChain+process.shiftedPatJetResDown+process.shiftedPatMETCorrJetResDown+process.shiftedPatJetResUp+process.shiftedPatMETCorrJetResUp+process.pfCandsNoJets+process.pfCandsNoJetsNoEle+process.pfCandsNoJetsNoEleNoMu+process.pfCandsNoJetsNoEleNoMuNoTau+process.pfCandsForUnclusteredUnc+process.pfElectrons+process.pfTaus+process.pfMuons+process.pfNoPileUp+process.pfPhotons+process.shiftedPatMuonEnDown+process.shiftedPatMETCorrMuonEnDown+process.shiftedPatMuonEnUp+process.shiftedPatMETCorrMuonEnUp+process.shiftedPatJetEnDown+process.shiftedPatMETCorrJetEnDown+process.shiftedPatJetEnUp+process.shiftedPatMETCorrJetEnUp+process.shiftedPatTauEnDown+process.shiftedPatMETCorrTauEnDown+process.shiftedPatTauEnUp+process.shiftedPatMETCorrTauEnUp+process.shiftedPatPhotonEnDown+process.shiftedPatMETCorrPhotonEnDown+process.shiftedPatPhotonEnUp+process.shiftedPatMETCorrPhotonEnUp+process.shiftedPatElectronEnDown+process.shiftedPatMETCorrElectronEnDown+process.shiftedPatElectronEnUp+process.shiftedPatMETCorrElectronEnUp+process.shiftedPatUnclusteredEnDown+process.shiftedPatMETCorrUnclusteredEnDown+process.shiftedPatUnclusteredEnUp+process.shiftedPatMETCorrUnclusteredEnUp+process.patPFMetT1SmearpatMetUncertaintySequence)


process.patShiftedModuleSequencePuppi = cms.Sequence(process.patPFMetT1JetResDownPuppi+process.patPFMetT1JetResUpPuppi+process.patPFMetT1MuonEnUpPuppi+process.patPFMetT1MuonEnDownPuppi+process.patPFMetT1JetEnUpPuppi+process.patPFMetT1JetEnDownPuppi+process.patPFMetT1TauEnDownPuppi+process.patPFMetT1TauEnUpPuppi+process.patPFMetT1PhotonEnDownPuppi+process.patPFMetT1PhotonEnUpPuppi+process.patPFMetT1ElectronEnUpPuppi+process.patPFMetT1ElectronEnDownPuppi+process.patPFMetT1UnclusteredEnUpPuppi+process.patPFMetT1UnclusteredEnDownPuppi+process.patPFMetT1SmearPuppipatShiftedModuleSequencePuppi)


process.egammaPostRecoSeq = cms.Sequence(process.egammaUpdatorSeq+process.egammaScaleSmearSeq+process.egammaVIDSeq+process.egammaPostRecoPatUpdatorSeq)


process.patShiftedModuleSequence = cms.Sequence(process.patPFMetT1JetResUp+process.patPFMetT1JetResDown+process.patPFMetT1MuonEnUp+process.patPFMetT1MuonEnDown+process.patPFMetT1JetEnUp+process.patPFMetT1JetEnDown+process.patPFMetT1TauEnDown+process.patPFMetT1TauEnUp+process.patPFMetT1PhotonEnUp+process.patPFMetT1PhotonEnDown+process.patPFMetT1ElectronEnUp+process.patPFMetT1ElectronEnDown+process.patPFMetT1UnclusteredEnUp+process.patPFMetT1UnclusteredEnDown+process.patPFMetT1SmearpatShiftedModuleSequence)


process.patMetCorrectionSequencePuppi = cms.Sequence(process.patPFMetT1T2CorrSequencePuppi+process.patPFMetT1Puppi+process.patPFMetTxyCorrSequencePuppi+process.patPFMetT1TxyPuppi+process.patPFMetTxyPuppi+process.patPFMetSmearCorrSequencePuppi+process.patPFMetT1SmearPuppi)


process.patMetCorrectionSequence = cms.Sequence(process.patPFMetT1T2CorrSequence+process.patPFMetT1+process.patPFMetTxyCorrSequence+process.patPFMetT1Txy+process.patPFMetTxy+process.patPFMetSmearCorrSequence+process.patPFMetT1Smear)


process.fullPatMetSequence = cms.Sequence(process.patMetModuleSequence+process.patMetCorrectionSequence+process.patMetUncertaintySequence+process.patShiftedModuleSequence+process.patCaloMet+process.slimmedMETs)


process.fullPatMetSequencePuppi = cms.Sequence(process.patMetModuleSequencePuppi+process.patMetCorrectionSequencePuppi+process.patMetUncertaintySequencePuppi+process.patShiftedModuleSequencePuppi+process.patCaloMet+process.slimmedMETsPuppi)


process.p = cms.Path(process.egammaPostRecoSeq+process.ecalBadCalibReducedMINIAODFilter+process.patJetCorrFactorsUpdatedJEC+process.updatedPatJetsUpdatedJEC+process.selectedUpdatedPatJetsUpdatedJEC+process.fullPatMetSequence+process.puppiMETSequence+process.fullPatMetSequencePuppi+process.rerunMvaIsolationSequence+process.slimmedTausNewID+process.demo)


process.schedule = cms.Schedule(*[ process.p ])
