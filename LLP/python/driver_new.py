import ROOT
import numpy
import os
import sys
import json
import random
import argparse

ROOT.gStyle.SetPaintTextFormat("3.0f")

NRGBs = 5
NCont = 8

stops = numpy.array( [0.00, 0.34, 0.61, 0.84, 1.00] )
red  = numpy.array( [0.00, 0.00, 0.87, 1.00, 0.51] )
green = numpy.array( [0.00, 0.81, 1.00, 0.20, 0.00] )
blue = numpy.array( [0.51, 1.00, 0.12, 0.00, 0.00] )

colWheelDark = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)

for i in range(NRGBs):
    red[i]=min(1,red[i]*1.1+0.25)
    green[i]=min(1,green[i]*1.1+0.25)
    blue[i]=min(1,blue[i]*1.1+0.25)

colWheel = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
ROOT.gStyle.SetNumberContours(NCont)
ROOT.gRandom.SetSeed(123)

colors=[]
def hex2rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16)/255.0 for i in range(0, lv, lv // 3))

def newColor(red,green,blue):
    newColor.colorindex+=1
    color=ROOT.TColor(newColor.colorindex,red,green,blue)
    colors.append(color)
    return color
    
newColor.colorindex=301

def getDarkerColor(color):
    darkerColor=newColor(color.GetRed()*0.6,color.GetGreen()*0.6,color.GetBlue()*0.6)
    return darkerColor
    
    
    

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--syst', dest='syst', default="nominal", type=str, action='store', help='Systematic')
parser.add_argument('-c','--ctau', dest='ctau', default='1', type=str, action='store', help='Ctau value')
parser.add_argument('-o','--output', dest='output', default='.', type=str, action='store', help='Output dir')
#parser.add_argument('-m','--mass', dest='mass', default=[], type=str, action='append', help='LLP:LSP mass pairs')

args = parser.parse_args()

def htVar(syst="nominal"):
    if syst=="jesUp":
        return "jesTotalUp_ht"
    elif syst=="jesDown":
        return "jesTotalDown_ht"
    elif syst=="jerUp":
        return "jerUp_ht"
    elif syst=="jerDown":
        return "jerDown_ht"
    return "nominal_ht"
    
def mhtVar(syst="nominal"):
    if syst=="jesUp":
        return "jesTotalUp_mht"
    elif syst=="jesDown":
        return "jesTotalDown_mht"
    elif syst=="jerUp":
        return "jerUp_mht"
    elif syst=="jerDown":
        return "jerDown_mht"
    return "nominal_mht"
    
def minPhiVar(syst="nominal"):
    if syst=="jesUp":
        return "jesTotalUp_minPhi"
    elif syst=="jesDown":
        return "jesTotalDown_minPhi"
    elif syst=="jerUp":
        return "jerUp_minPhi"
    elif syst=="jerDown":
        return "jerDown_minPhi"
    return "nominal_minPhi"
    
def nJetVar(syst="nominal"):
    if syst=="jesUp":
        return "nselectedJets_jesTotalUp"
    elif syst=="jesDown":
        return "nselectedJets_jesTotalDown"
    elif syst=="jerUp":
        return "nselectedJets_jerUp"
    elif syst=="jerDown":
        return "nselectedJets_jerDown"
    return "nselectedJets_nominal"
    
def nVetoFwdJetVar(syst="nominal"):
    if syst=="jesUp":
        return "nvetoFwdJets_jesTotalUp"
    elif syst=="jesDown":
        return "nvetoFwdJets_jesTotalDown"
    elif syst=="jerUp":
        return "nvetoFwdJets_jerUp"
    elif syst=="jerDown":
        return "nvetoFwdJets_jerDown"
    return "nvetoFwdJets_nominal"
    
def metVar(syst="nominal"):
    if syst=="jesUp":
        return "jesTotalUp_met"
    elif syst=="jesDown":
        return "jesTotalDown_met"
    elif syst=="jerUp":
        return "jerUp_met"
    elif syst=="jerDown":
        return "jerDown_met"
    elif syst=="unclEnUp":
        return "unclEnUp_met"
    elif syst=="unclEnDown":
        return "unclEnDown_met"
    return "nominal_met"
    
def llpdnnxPrefix(syst,ctau):
    ctauHack = ctau
    if ctauHack=="1":
        ctauHack=="0"
    if syst=="jesUp":
        return "llpdnnx_da_jesUp_%s"%(ctau)
    elif syst=="jesDown":
        return "llpdnnx_da_jesDown_%s"%(ctau)
    elif syst=="jerUp":
        return "llpdnnx_da_jerUp_%s"%(ctau)
    elif syst=="jerDown":
        return "llpdnnx_da_jerDown_%s"%(ctau)
    return "llpdnnx_da_nominal_%s"%(ctau)
    
def llpdnnxVar(syst,ctau,m=0):
    return llpdnnxPrefix(syst,ctau)+"_LLP_min%i"%(m)
    
    
def llpdnnxNLLPtrue(syst,ctau):
    return llpdnnxPrefix(syst,ctau)+"_nLLPTrue"

def llpdnnxNLLPtrueTaggedLLP(syst,ctau):
    return llpdnnxPrefix(syst,ctau)+"_nLLPTrueTaggedLLP"
    

def baseSelection(syst="nominal",region="D"):
    # B | D   
    # A | C   => C/A=D/B => D = B*C/A
    selection = "((MonoCentralPFJet80_PFMETNoMu_PFMHTNoMu_IDTight>0) || (PFHT900>0))"
    selection += "*(nvetoMuons==0)*(nvetoElectrons==0)"
    selection += "*("+nVetoFwdJetVar(syst)+"==0)"
    selection += "*("+nJetVar(syst)+">2)"
    selection += "*(("+mhtVar(syst)+"/"+metVar(syst)+")<1.25)"
    if region=="A":
        selection += "*("+mhtVar(syst)+">270.)*("+mhtVar(syst)+"<300.)"
        selection += "*("+minPhiVar(syst)+"<0.2)"
    elif region=="B":
        selection += "*("+mhtVar(syst)+">270.)*("+mhtVar(syst)+"<300.)"
        selection += "*("+minPhiVar(syst)+">0.2)"
    elif region=="C":
        selection += "*("+mhtVar(syst)+">300.)"
        selection += "*("+minPhiVar(syst)+"<0.2)"
    elif region=="D":
        selection += "*("+mhtVar(syst)+">300.)"
        selection += "*("+minPhiVar(syst)+">0.2)"
    else:
        print "Error - region needs to be A,B,C or D"
        sys.exit(1)
    return selection
    
def mcWeight(sampleName,syst="nominal",lumi="35.822"):
    mcWeight = "(genweight)*"+lumi+"*1000.0"
    if args.syst=="puUp":
        mcWeight+="*(puweightUp)"
    elif args.syst=="puDown":
        mcWeight+="*(puweightDown)"
    else:
        mcWeight+="*(puweight)"
    if sampleName=="WJets":
        if args.syst=="wjetsScaleUp":
            mcWeight+="*TMath::Range(-2,2,scaleweight_0/scaleweight_4)"
        elif args.syst=="wjetsScaleDown":
            mcWeight+="*TMath::Range(-2,2,scaleweight_8/scaleweight_4)"
    elif sampleName=="ttbar":
        if args.syst=="ttbarScaleUp":
            mcWeight+="*TMath::Range(-2,2,scaleweight_0/scaleweight_4)"
        elif args.syst=="ttbarScaleDown":
            mcWeight+="*TMath::Range(-2,2,scaleweight_8/scaleweight_4)"
    elif sampleName=="st":
        if args.syst=="stScaleUp":
            mcWeight+="*TMath::Range(-2,2,scaleweight_0/scaleweight_4)"
        elif args.syst=="stScaleDown":
            mcWeight+="TMath::Range(-2,2,scaleweight_8/scaleweight_4)"
    elif sampleName=="ZNuNu":
        if args.syst=="znunuScaleUp":
            mcWeight+="*TMath::Range(-2,2,scaleweight_0/scaleweight_4)"
        elif args.syst=="znunuScaleDown":
            mcWeight+="*TMath::Range(-2,2,scaleweight_8/scaleweight_4)"
    elif sampleName=="QCDHT":
        pass
    elif sampleName=="SMS-T1qqqq":
        pass
    else:
        raise Exception("Unknown sample name '"+sampleName+"'")
    return mcWeight
    

globalMCWeight = lambda sampleName,region="D": mcWeight(sampleName,args.syst)+"*"+baseSelection(args.syst,region)
globalDataWeight = lambda region="D": baseSelection(args.syst,region)

filePath = "/vols/cms/mkomm/LLP/NANOX_SR"


extensions = {
    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8":"WToLNu_0J_13TeV-amcatnloFXFX-pythia8",
    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8_ext1":"WToLNu_0J_13TeV-amcatnloFXFX-pythia8",
    
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8":"WToLNu_1J_13TeV-amcatnloFXFX-pythia8",
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8_ext1":"WToLNu_1J_13TeV-amcatnloFXFX-pythia8",
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8_ext3":"WToLNu_1J_13TeV-amcatnloFXFX-pythia8",
    
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8_ext1":"WToLNu_2J_13TeV-amcatnloFXFX-pythia8",
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8_ext2":"WToLNu_2J_13TeV-amcatnloFXFX-pythia8",
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8_ext3":"WToLNu_2J_13TeV-amcatnloFXFX-pythia8",
    
    "DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":"DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1":"DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext4":"DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",

    "DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":"DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1":"DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext4":"DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    
    "DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":"DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1":"DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",

    "DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":"DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1":"DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",    

}

xsecs = {

    # Drell-Yan + jets
    # NLO
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    "DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":169.9*1.23,
    "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":147.40*1.23,
    "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":40.99*1.23,
    "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":5.678*1.23,
    "DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":1.367*1.23,
    "DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.6304*1.23,
    "DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.1514*1.23,
    "DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.003565*1.23,
    
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":18610,
    "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":1921.8*3,
    
    # QCD (multijet)
    # LO
    "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":27990000.,
    "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":1712000.,
    "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":347700.,
    "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":32100.,
    "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":6831.,
    "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":1207.,
    "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":119.9,
    "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":25.24,
    
    #https://cms-pdmv.cern.ch/mcm/requests?page=-1&dataset_name=QCD_Pt_*to*_TuneCUETP8M1_13TeV_pythia8&member_of_campaign=RunIIFall14GS
    "QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8":140932000,
    "QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8":19204300,
    "QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8":2762530,
    "QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8":471100,
    "QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8":117276,
    "QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8":7823,
    "QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8":648.2,
    "QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8":186.9,
    "QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8":32.293,
    "QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8":9.4183,
    "QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8":0.84265,
    "QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8":0.114943,
    "QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8":0.00682981,
    "QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8":0.000165445,
    
    
    # Single-top
    # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec#Single_top_t_channel_cross_secti
    # NLO
    "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1": 10.32,
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1": 80.95,
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1": 136.02,
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1": 71.7/2.,
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1": 71.7/2.,
    
    # TTbar
    # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO (mtop=172.5 GeV)
    # missing
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen": 831.76,
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8": 831.76,
    # NLO
    "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8": 831.76,
    # LO
    "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 831.76,
    "TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 2.666535,
    "TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1.098082,
    "TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 0.198748, 
    "TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 0.002368413,
    "TTJets_SingleLeptFromTbar_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8": 831.76*0.5*(1-3*0.1086)*(3*0.1086)*2,
    "TTJets_SingleLeptFromT_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8": 831.76*0.5*(1-3*0.1086)*(3*0.1086)*2,

    # W->l nu + jets
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    # NLO
    "WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1319,#*1.21, 
    "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1345,#*1.21, 
    "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 359.7,#*1.21, 
    "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 48.91,#*1.21, 
    "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 12.05,#*1.21, 
    "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 5.501,#*1.21, 
    "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1.329,#*1.21, 
    "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 0.03216,#*1.21, 
    
    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8": 49670.,
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8": 8264.,
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8": 3226.,
    
    # Z -> nu nu + jets
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
    # LO
    "ZJetsToNuNu_HT-100To200_13TeV-madgraph": 280.35*1.23 ,
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph": 77.67*1.23 ,
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph": 10.73*1.23 ,
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph": 2.559*1.23 ,
    "ZJetsToNuNu_HT-800to1200_13TeV-madgraph": 1.1796*1.23 ,
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph": 0.28833*1.23 ,
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph": 0.006945*1.23 ,
    
    # Z -> nu nu + jets
    # MCM
    # NLO
    "DYJetsToNuNu_Zpt-0To50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 3483.*3,
    "DYJetsToNuNu_PtZ-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 593.9*3,
    "DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 55.03*3,
    "DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 2.082*3,
    "DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 0.2816*3,
    "DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8": 0.02639*3,
    
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    # NNLO
    "WW_TuneCUETP8M1_13TeV-pythia8": 118.7,
    
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    # NLO
    "WZ_TuneCUETP8M1_13TeV-pythia8": 47.13,
    "ZZ_TuneCUETP8M1_13TeV-pythia8": 16.523,
}

with open('eventyields.json',) as f:
    genweights = json.load(f)
    
processDict = {}

for processFolder in os.listdir(filePath):
    print "reading ",processFolder,"...",
    processName = processFolder
    if extensions.has_key(processFolder):
        processName = extensions[processFolder]
        
    if not processDict.has_key(processName):
        processDict[processName] = {
            "files":[],
            "weight":"weightnotfound",
            "xsec":0,
            "integral":0
        }
    for f in os.listdir(os.path.join(filePath,processFolder)):
        fullFilePath = os.path.join(filePath,processFolder,f)
        '''
        if not f.endswith(".root"):
            continue
        
        rootFile = ROOT.TFile(fullFilePath)
        if not rootFile:
            continue
        tree = rootFile.Get("Friends")
        if not tree:
            continue
        '''
        
        
        processDict[processName]["files"].append(fullFilePath)
        '''
        if tree.FindBranch("genweight"):
            h = ROOT.TH1F("nevents"+processFolder+f,"",1,-1,1)
            tree.Project(h.GetName(),"0","genweight")
            processDict[processFolder]["nevents"]+=tree.GetEntries()
            processDict[processFolder]["integral"]+=h.Integral()
        else:
            processDict[processFolder]["nevents"]+=tree.GetEntries()
            processDict[processFolder]["integral"]+=tree.GetEntries()
        rootFile.Close()
        #break
        '''
    print len(processDict[processName]["files"]),
   
    #if xsecs.has_key(processFolder) and genweights.has_key(processFolder):
    if xsecs.has_key(processName):
        print genweights[processFolder]["sum"],genweights[processFolder]["weighted"]
        processDict[processName]["xsec"] = xsecs[processName]
        processDict[processName]["integral"] += genweights[processFolder]["weighted"]
    else: 
        print 'noxsec'
    
for processName in sorted(processDict.keys()):
    if processDict[processName]["integral"]>0:
        processDict[processName]["weight"] = str(processDict[processName]["xsec"]/processDict[processName]["integral"])
    print processName,len(processDict[processName]["files"]),processDict[processName]["weight"]
         
#globalWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(IsoMuTrigger_flag>0)*(nvetoFwdJets_nominal==0)"
#globalWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*((MonoCentralPFJet80_PFMETNoMu_PFMHTNoMu_IDTight>0) || (PFHT900>0))"
#globalWeight += "*()"
#globalMCWeight = globalWeight+"*(genweight)*35.822*1000.0*(puweight)"
#globalMCWeight = globalMCWeight+"*(tightMuons_weight_track_nominal*tightMuons_weight_id_nominal*tightMuons_weight_iso_nominal)"
#globalDataWeight = globalWeight#+"*(IsoMuTrigger_flag==1)"




mcSetDict = {
    "diboson": {
        "processes": [
            "WW_TuneCUETP8M1_13TeV-pythia8",
            "WZ_TuneCUETP8M1_13TeV-pythia8",
            "ZZ_TuneCUETP8M1_13TeV-pythia8"
        ],
        "weight":'1',
        "fill":newColor(0.35,0.95,0.55),
        "title":"Diboson",
        "addtitle":""
    },
    "WJets": {
        "processes": [
            "WToLNu_0J_13TeV-amcatnloFXFX-pythia8",
            "WToLNu_1J_13TeV-amcatnloFXFX-pythia8",
            "WToLNu_2J_13TeV-amcatnloFXFX-pythia8"
        ],
        "weight":'1',
        "fill":newColor(0.33,0.75,0.35),
        "title":"W+jets",
        "addtitle":""
    },
    "WJetsHT":{
        "processes":[
            "WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
        ],
        "weight":'1',
        "fill":newColor(0.33,0.75,0.35),
        "title":"W+jets",
        "addtitle":""
    },
    "ZJets": {
        "processes": [
            "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"
        ],
        "weight":'1',
        "fill":newColor(0.3,0.75,0.95),
        "title":"Z/#gamma*+jets",
        "addtitle":""
    },
    "ZJetsHT": {
        "processes": [
            "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8","(genHt<100.)"],
            ["DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
            ["DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","(znloweight_nominal)"],
        ],
        "weight":'1',
        "fill":newColor(0.3,0.75,0.95),
        "title":"Z/#gamma*+jets",
        "addtitle":"(rew. NLO)"
    },
    "ttbar": {
        "processes": [
            "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",
        ],
        "weight":'1',
        "fill":newColor(1.,0.8,0.0),
        "title":"t#bar{t}",
        "addtitle":""
    },

    "st":{
        "processes": [
            #"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1",
            "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
            "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
            "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
            "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"
        ],
        "weight":'1',
        "fill":newColor(0.5,0.47,0.95),
        "title":"Single t/#bar{t}",
        "addtitle":""
    },
    
    
    "QCD": {
        "processes": [
            "QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8",
            "QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8",
        ],
        "weight":'1',
        "fill":newColor(0.85,0.85,0.85),
        "title":"Multijet",
        "addtitle":""
    },
    

    "QCDHT": {
        "processes": [
            "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        ],
        "weight": '1',
        "fill":newColor(0.85,0.85,0.85),
        "title":"Multijet",
        "addtitle":""
    },
    
    "ZNuNuHT": {
        "processes": [
            "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
            "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
            "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
            "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
            "ZJetsToNuNu_HT-800to1200_13TeV-madgraph",
            "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
            "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph",
        ],
        "weight": '1',
        "fill":newColor(0.3,0.75,0.95),
        "title":"Z#rightarrow#nu#nu",
        "addtitle":""
    },
    
    "ZNuNu": {
        "processes": [
            "DYJetsToNuNu_Zpt-0To50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToNuNu_PtZ-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
        ],
        "weight": '1',
        "fill":newColor(0.3,0.75,0.95),
        "title":"Z#rightarrow#nu#nu",
        "addtitle":""
    },
}



ctauHack = args.ctau
if args.ctau=="1":
    ctauHack = "0"

thresholds = {      
    "0p001": 0.37245617939,
    "0p01": 0.368606702814,
    "0p1": 0.343856976757,
    "0": 0.436759756336,
    "10": 0.557319806347,
    "100": 0.64487092329,
    "1000": 0.565088437328,
    "10000": 0.590826915035
}

signalConfigs = []
massesDict = {}
ctauSampleName = "SMS-T1qqqq_ctau-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"%args.ctau
for signalSample in [ctauSampleName,ctauSampleName+"_extra"]:
    for llpMass in genweights[signalSample].keys():
        for lspMass in genweights[signalSample][llpMass].keys():
            if genweights[signalSample][llpMass][lspMass]["sum"]<10:
                continue
            if int(llpMass)==int(lspMass):
                continue
            if not massesDict.has_key(llpMass):
                massesDict[llpMass] = {}
            if not massesDict[llpMass].has_key(lspMass):
                massesDict[llpMass][lspMass] = 0.
            massesDict[llpMass][lspMass]+=genweights[signalSample][llpMass][lspMass]["weighted"]

for llpMass in sorted(massesDict.keys()):
    print llpMass,":",
    for lspMass in sorted(massesDict[llpMass]):
        genweight = 1./float(massesDict[llpMass][lspMass])
        print lspMass,
        signalConfigs.append({
            "name":"ctau%s_llp%i_lsp%i"%(args.ctau,int(llpMass),int(lspMass)),
            "llp":llpMass,
            "lsp":lspMass,
            "processes":[ctauSampleName,ctauSampleName+"_extra"],
            "weight":"(%7.5e)*(llp==%i)*(lsp==%i)"%(genweight,int(llpMass),int(lspMass))
        })
    print



varBinned=""
binindex = 1
htCut = "900"

if args.ctau=="1000" or args.ctau=="10000":
    htCut = "800."

for cat in range(2):
    if cat==0:
        varBinned += "("+htVar(args.syst)+">"+htCut+")*("
    if cat==1:
        varBinned += "("+htVar(args.syst)+"<"+htCut+")*("
        
    for t in [3,5,6]:
        if t==3:
            for m in [2]:
                taggerCutExp=str(binindex)+"*("+nJetVar(args.syst)+"<4.5)"
                binindex+=1
                if m==2:
                    taggerCutExp += "*("+llpdnnxVar(args.syst,ctauHack,1)+">"+str(thresholds[ctauHack])+")"
                varBinned+=taggerCutExp+"+"
        elif t==5:
            for m in [2]:
                taggerCutExp=str(binindex)+"*("+nJetVar(args.syst)+">4.5)*("+nJetVar(args.syst)+"<5.5)"
                binindex+=1
                if m==2:
                    taggerCutExp += "*("+llpdnnxVar(args.syst,ctauHack,1)+">"+str(thresholds[ctauHack])+")"
                varBinned+=taggerCutExp+"+"
        elif t==6:
            for m in [3]:
                taggerCutExp=str(binindex)+"*("+nJetVar(args.syst)+">5.5)"
                binindex+=1
                if m==3:
                    taggerCutExp += "*("+llpdnnxVar(args.syst,ctauHack,2)+">"+str(thresholds[ctauHack])+")"
                varBinned+=taggerCutExp+"+"
    varBinned+="0)+"
varBinned+="0"

    
nbins = binindex-1


histBinning= {
    "name":"llp",
    "var":varBinned,
    "binning":numpy.linspace(-0.5,nbins+0.5,num=nbins+2),
    "weight":"1"
}


print histBinning


fileName = os.path.join(args.output,"hist_%s_%s_%s.root"%(args.ctau,histBinning["name"],args.syst))
if os.path.exists(fileName):
    print "Output exists -> skip"
    sys.exit(0)
    
histograms = []

def removeNegEntries(hist,avgWeight=-1):
    alpha = 1. - 0.6827
    upErr = ROOT.Math.gamma_quantile_c(alpha/2,1,1)
    avgWeight = hist.Integral()/hist.GetEntries() if hist.GetEntries()>0 else -1
    #print "weight",avgWeight
    for ibin in range(hist.GetNbinsX()):
        c = hist.GetBinContent(ibin+1)
        if c<10**-4:
            hist.SetBinContent(ibin+1,10**-3)
            #note: in case of 0 entries the uncertainty is also small
            #(this is not the case with negative events)
            if hist.GetBinError(ibin+1)<10**-4 and avgWeight>0:
                #set uncertainties for empy bins
                #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
                hist.SetBinError(ibin+1,upErr*avgWeight)
            else:
                hist.SetBinError(ibin+1,10**-4)
        #print "bin%2i, %.1f+-%.1f (+-%.1f%%)"%(ibin,c,hist.GetBinError(ibin+1),100.*hist.GetBinError(ibin+1)/c if c>0 else -1)
                
systPostfix = "_"+args.syst
if args.syst=="nominal":
    systPostfix=""
    
histBinningWeight = histBinning["weight"]
print "bin: ",histBinning["name"]




for region in ["A","B","C","D"]:
    histBackgroundSum = ROOT.TH1F(histBinning["name"]+region+"_sum"+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
    histBackgroundSum.SetDirectory(0)
    histograms.append(histBackgroundSum)
    for backgroundMCSet in ["WJets","st","ttbar","ZNuNu","QCDHT"]:
        print "  background: ",backgroundMCSet,region
        setWeight = globalMCWeight(backgroundMCSet,region=region)
        histBackground = ROOT.TH1F(histBinning["name"]+region+"_"+backgroundMCSet+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
        histBackground.SetDirectory(0)
        histograms.append(histBackground)
        histBackground.Sumw2()
        
        sumPos = 0
        sumNeg = 0
        integralPos = 0
        integralNeg = 0
        
        for process in mcSetDict[backgroundMCSet]["processes"]:
            processWeight = "1"
            if type(process)==type(list()):
                processWeight+="*"+process[1]
                process = process[0]
            processWeight+="*"+processDict[process]["weight"]
            #print histBinningWeight+"*"+setWeight+"*"+processWeight
            for ifile,f in enumerate(processDict[process]["files"]):
                rootFile = ROOT.TFile.Open(f)
                tree = rootFile.Get("Friends")
                if tree:
                    processHistPos = ROOT.TH1F(histBinning["name"]+"_pos_"+backgroundMCSet+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                    tree.Project(processHistPos.GetName(),histBinning["var"],histBinningWeight+"*"+setWeight+"*"+processWeight+"*(genweight>0)")
                    
                    processHistNeg = ROOT.TH1F(histBinning["name"]+"_neg_"+backgroundMCSet+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                    tree.Project(processHistNeg.GetName(),histBinning["var"],histBinningWeight+"*"+setWeight+"*"+processWeight+"*(genweight<0)")
                    
                    processHistPos.SetDirectory(0)
                    processHistNeg.SetDirectory(0)
                    histBackground.Add(processHistPos)
                    histBackground.Add(processHistNeg)
                    
                    sumPos+=processHistPos.GetEntries()
                    sumNeg+=processHistNeg.GetEntries()
                    integralPos+=processHistPos.Integral()
                    integralNeg+=processHistNeg.Integral()
                    
                rootFile.Close()
                break
        
        removeNegEntries(
            histBackground,
            avgWeight=(integralPos-integralNeg)/(sumPos+sumNeg) if (sumPos+sumNeg)>0 else -1
        )
        histBackgroundSum.Add(histBackground)
        print "    -> ",histBackground.GetEntries(),"/",histBackground.Integral()," entries/integral"
        
        
if args.syst=="nominal":
    for region in ["A","B","C","D"]:
        print "  data: ",region
        setWeight = globalDataWeight(region=region)
        histData = ROOT.TH1F(histBinning["name"]+region+"_data"+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
        histData.SetDirectory(0)
        histograms.append(histData)
        histData.Sumw2()
        
        sumPos = 0
        sumNeg = 0
        integralPos = 0
        integralNeg = 0
        
        for process in mcSetDict[backgroundMCSet]["processes"]:
            processWeight = "1"
            if type(process)==type(list()):
                processWeight+="*"+process[1]
                process = process[0]
            processWeight+="*"+processDict[process]["weight"]
            #print histBinningWeight+"*"+setWeight+"*"+processWeight
            for ifile,f in enumerate(processDict[process]["files"]):
                rootFile = ROOT.TFile.Open(f)
                tree = rootFile.Get("Friends")
                if tree:
                    processHistPos = ROOT.TH1F(histBinning["name"]+"_pos_"+backgroundMCSet+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                    tree.Project(processHistPos.GetName(),histBinning["var"],histBinningWeight+"*"+setWeight+"*"+processWeight+"*(genweight>0)")
                    
                    processHistNeg = ROOT.TH1F(histBinning["name"]+"_neg_"+backgroundMCSet+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                    tree.Project(processHistNeg.GetName(),histBinning["var"],histBinningWeight+"*"+setWeight+"*"+processWeight+"*(genweight<0)")
                    
                    processHistPos.SetDirectory(0)
                    processHistNeg.SetDirectory(0)
                    histBackground.Add(processHistPos)
                    histBackground.Add(processHistNeg)
                    
                    sumPos+=processHistPos.GetEntries()
                    sumNeg+=processHistNeg.GetEntries()
                    integralPos+=processHistPos.Integral()
                    integralNeg+=processHistNeg.Integral()
                    
                rootFile.Close()
                break
        
        removeNegEntries(
            histBackground,
            avgWeight=(integralPos-integralNeg)/(sumPos+sumNeg) if (sumPos+sumNeg)>0 else -1
        )
        histBackgroundSum.Add(histBackground)
        print "    -> ",histBackground.GetEntries(),"/",histBackground.Integral()," entries/integral"
        

for signalCfg in signalConfigs:
    for region in ["D"]:
        signalHistSum = ROOT.TH1F(histBinning["name"]+region+"_"+signalCfg["name"]+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
        signalHistSum.SetDirectory(0)
        histograms.append(signalHistSum)
        
        for nLLP in range(7):
            for nLLPTagged in range(nLLP+1):
                print "  signal: ",signalCfg["name"],region,nLLP,nLLPTagged
                signalHist = ROOT.TH1F(histBinning["name"]+region+"_"+signalCfg["name"]+("_%illp%it"%(nLLP,nLLPTagged))+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
                setWeight = globalMCWeight("SMS-T1qqqq",region)
                signalHist.SetDirectory(0)
                histograms.append(signalHist)
                
                signalWeight = signalCfg["weight"]
                if nLLP<6:
                    signalWeight+=("*("+llpdnnxNLLPtrue(args.syst,ctauHack)+"==%i)*("+llpdnnxNLLPtrueTaggedLLP(args.syst,ctauHack)+"==%i)")%(nLLP,nLLPTagged)
                else:
                    if nLLPTagged<nLLP:
                        signalWeight+=("*("+llpdnnxNLLPtrue(args.syst,ctauHack)+">=%i)*("+llpdnnxNLLPtrueTaggedLLP(args.syst,ctauHack)+"==%i)")%(nLLP,nLLPTagged)
                    else:
                        signalWeight+=("*("+llpdnnxNLLPtrue(args.syst,ctauHack)+">=%i)*("+llpdnnxNLLPtrueTaggedLLP(args.syst,ctauHack)+">=%i)")%(nLLP,nLLPTagged)
                #print histBinningWeight+"*"+signalWeight+"*"+setWeight
                for process in signalCfg["processes"]:
                    for ifile,f in enumerate(processDict[process]["files"]):
                        rootFile = ROOT.TFile.Open(f)
                        tree = rootFile.Get("Friends")
                        if tree:
                            processHist = ROOT.TH1F(histBinning["name"]+"_"+signalCfg["name"]+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                            tree.Project(processHist.GetName(),histBinning["var"],histBinningWeight+"*"+signalWeight+"*"+setWeight)
                            processHist.SetDirectory(0)
                            signalHist.Add(processHist)
                            signalHistSum.Add(processHist)
                        rootFile.Close()
                        #break
                print "    -> ",signalHist.GetEntries(),"/",signalHist.Integral()," entries/integral"
                removeNegEntries(
                    signalHist,
                    avgWeight=-1
                )
    removeNegEntries(
        signalHistSum,
        avgWeight=signalHistSum.Integral()/(10**-12+signalHistSum.GetEntries())
    )
outputFile = ROOT.TFile(fileName,"RECREATE")
for hist in histograms:
    print "writing hist: ",hist.GetName()
    hist.SetDirectory(outputFile)
    hist.Write()
outputFile.Write()
outputFile.Close()


