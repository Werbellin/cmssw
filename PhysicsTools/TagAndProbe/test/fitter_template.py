import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing('analysis')

options.register(
    "isMC",
    XXX_IS_MC_XXX,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Compute MC efficiencies"
    )

options.register(
    "inputFileName",
    "XXX_INPUTFILE_XXX",
#    "TnPTree_mc.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input filename"
    )

options.register(
    "outputFileName",
    "test",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output filename"
    )

options.register(
    "idName",
#    "passingTight",
    "passingALL",
#    "probe_ele_Spring16_HZZ_WP",
    #"passingTrigWP90",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ID variable name as in the fitter_tree"
    )

options.register(
    "dirName",
    "GsfElectronToEleID",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Folder name containing the fitter_tree"
    )

options.register(
    "doCutAndCount",
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Perform cut and count efficiency measurement"
    )

options.parseArguments()


process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

################################################

InputFileName = options.inputFileName
OutputFile = "efficiency-mc-"+options.idName
if (not options.isMC):
    OutputFile = "efficiency-data-"+options.idName

if (options.outputFileName != ""):
    OutputFile = OutputFile+"-"+options.outputFileName+".root"
else:
    OutputFile = OutputFile+".root"

################################################

EfficiencyBins = cms.PSet(
    probe_Ele_abseta = cms.vdouble(XXX_eta_low_XXX ,  XXX_eta_high_XXX),
    probe_Ele_pt = cms.vdouble(XXX_pt_low_XXX, XXX_pt_high_XXX),
    )

EfficiencyBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass", "totWeight"),
    BinnedVariables = cms.PSet(EfficiencyBins,
                               mcTrue = cms.vstring("true")
                               ),
    BinToPDFmap = cms.vstring("pdfSignalPlusBackground")  
    )

if (not options.isMC):
    EfficiencyBinningSpecification.UnbinnedVariables = cms.vstring("mass")
    EfficiencyBinningSpecification.BinnedVariables = cms.PSet(EfficiencyBins)

mcTruthModules = cms.PSet()
if (options.isMC):
    setattr(mcTruthModules, "MCtruth_" + options.idName, cms.PSet(EfficiencyBinningSpecification))
    setattr(getattr(mcTruthModules, "MCtruth_" + options.idName), "EfficiencyCategoryAndState", cms.vstring(options.idName, "pass"))

############################################################################################

process.TnPMeasurement = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                        InputFileNames = cms.vstring(InputFileName),
                                        InputDirectoryName = cms.string(options.dirName),
                                        InputTreeName = cms.string("fitter_tree"), 
                                        OutputFileName = cms.string(OutputFile),
                                        NumCPU = cms.uint32(1),
                                        SaveWorkspace = cms.bool(False), #VERY TIME CONSUMING FOR MC
                                        doCutAndCount = cms.bool(options.doCutAndCount),
                                        floatShapeParameters = cms.bool(True),
                                        binnedFit = cms.bool(True),
                                        binsForFit = cms.uint32(60),
                                        #WeightVariable = cms.string("totWeight"),
                                        #fixVars = cms.vstring("meanP", "meanF", "sigmaP", "sigmaF", "sigmaP_2", "sigmaF_2"),
                                        
                                        # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
                                        Variables = cms.PSet(mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
                                                             #mass = cms.vstring("Tag-Probe Mass 2", "60.0", "120.0", "GeV/c^{2}"),

                                                             #probe_Ele_et = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
                                                             probe_Ele_abseta = cms.vstring("Probe #eta", "-2.5", "2.5", ""), 
                                                             totWeight = cms.vstring("totWeight", "-1000000000", "100000000", ""),
                                                             #probe_Ele_e = cms.vstring("probe_Ele_e", "0", "1000", ""),
                                                             probe_Ele_pt = cms.vstring("probe_Ele_pt", "0", "1000", ""),
                                                             #probe_Ele_trigMVA = cms.vstring("probe_Ele_trigMVA", "-1", "1", ""),
                                                             #passingTrigWP90 = cms.vstring("passingTrigWP90", "-1", "1", ""),
                                                             #probe_Ele_isGap = cms.vstring("probe_ele_isEBEtaGap", "0", "2", ""),
                                                             ),
                                        
                                        # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
                                        Categories = cms.PSet(),

                                        Expressions = cms.PSet(myeop = cms.vstring("myeop", "probe_Ele_e/probe_Ele_pt", "probe_Ele_e", "probe_Ele_pt") 
                                                               ), 
                                        Cuts = cms.PSet(#gap_cut = cms.vstring("probe_Ele_isGap", "0.5", "XXX_GAP_CUT_XXX"),
                                                        #fakeEoPCut = cms.vstring("myeop", "2.", "above")
                                                        ),
                                        # defines all the PDFs that will be available for the efficiency calculations; 
                                        # uses RooFit's "factory" syntax;
                                        # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" 
                                        # and "signalFractionInPassing[0.9]" are used for initial values  
                                        PDFs = cms.PSet(pdfSignalPlusBackground = cms.vstring( XXX_FIT_STRING_XXX,


# "RooCBExGaussShape::signalResPass(mass, X_meanP_X, X_sigmaP_X, X_alphaP_X, X_nP_X, X_sigmaP_2_X)",
# "RooCBExGaussShape::signalResFail(mass, X_meanF_X, X_sigmaF_X, X_alphaF_X, X_nF_X, X_sigmaF_2_X)",
# 
# "RooBreitWigner::signalBWPass(mass,massP[91.1876],gammaP[2.4952])",
# "RooBreitWigner::signalBWFail(mass,massP[91.1876],gammaP[2.4952])",
#
# "FCONV::signalPass(mass, signalBWPass, signalResPass)",
# "FCONV::signalFail(mass, signalBWFail, signalResFail)",
# 
#"RooCMSShape::backgroundPass(mass, alphaPass[45.,40.,90.], betaPass[0.,0,0.1], gammaPass[0.1, 0., 1.], peakPass[91.1876])",
#"RooCMSShape::backgroundFail(mass, alphaFail[45.,40.,90.], betaFail[0.,0,0.1], gammaFail[0.1, 0., 1.], peakFail[91.1876])",


                                                                                              #"RooExponential::signalTailPass(mass, deltaP[0.])",
                                                                                              #"RooExponential::signalTailFail(mass, deltaF[0.])",

                                                                                              #"SUM::signalPass(cPass[0.]*signalTailPass,signal1Pass)",
                                                                                              #"SUM::signalFail(cFail[0.]*signalTailFail,signal1Fail)",

                                                                                              #"RooExponential::backgroundPass(mass, deltaBP[-0.55,-2.,-0.02])",
                                                                                              #"RooExponential::backgroundFail(mass, deltaBF[-0.55,-2.,-0.02])",
                                                                                              "efficiency[0.9,0.,1.]",
                                                                                              "signalFractionInPassing[1.0]",

#            "RooCBExGaussShape::signalResFail(mass,meanF[-0.0,-5.000,5.000],sigmaF[3.331,0.00,15.000],alphaF[1.586, 0.0,50.0],nF[0.464,0.000,20.00], sigmaF_2[1.675,0.500,12.000])",
#            "ZGeneratorLineShape::signalPhy(mass)", 
#            "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], gammaPass[0.1, 0, 1], peakPass[90.0])",
#            "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], gammaFail[0.1, 0, 1], peakFail[90.0])",
#            "FCONV::signalPass(mass, signalPhy, signalResPass)",
#            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
#            "efficiency[0.5,0,1]",
#            "signalFractionInPassing[1.]"     
            ),
                                                        ),
                                        
                                        # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
                                        # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
                                        Efficiencies = cms.PSet(mcTruthModules)
                                        )

if options.isMC :
   process.TnPMeasurement.WeightVariable = cms.string("totWeight") 

setattr(process.TnPMeasurement.Categories, options.idName, cms.vstring(options.idName, "dummy[pass=1,fail=0]"))
setattr(process.TnPMeasurement.Categories, "mcTrue", cms.vstring("MC true", "dummy[true=1,false=0]"))

if (not options.isMC):
    for pdf in process.TnPMeasurement.PDFs.__dict__:
        param =  process.TnPMeasurement.PDFs.getParameter(pdf)
        if (type(param) is not cms.vstring):
            continue
        for i, l in enumerate(getattr(process.TnPMeasurement.PDFs, pdf)):
            if l.find("signalFractionInPassing") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = l.replace("[1.0]","[0.5,0.,1.]")
        
    setattr(process.TnPMeasurement.Efficiencies, options.idName, EfficiencyBinningSpecification)    
    setattr(getattr(process.TnPMeasurement.Efficiencies, options.idName) , "EfficiencyCategoryAndState", cms.vstring(options.idName, "pass"))
else:
    for pdf in process.TnPMeasurement.PDFs.__dict__:
        param =  process.TnPMeasurement.PDFs.getParameter(pdf)
        if (type(param) is not cms.vstring):
            continue
        for i, l in enumerate(getattr(process.TnPMeasurement.PDFs, pdf)):
            if l.find("backgroundPass") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = "RooPolynomial::backgroundPass(mass, a[0.0])"
            if l.find("backgroundFail") != -1:
                getattr(process.TnPMeasurement.PDFs, pdf)[i] = "RooPolynomial::backgroundFail(mass, a[0.0])"

process.fit = cms.Path(
    process.TnPMeasurement  
    )

