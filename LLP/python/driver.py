import ROOT
import numpy
import os
import sys
import json
import random
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--jetBin', dest='jetBin', default=0, type=int, action='store', help='Jet bin')
parser.add_argument('--htBin', dest='htBin', default=0, type=int, action='store', help='HT bin')
parser.add_argument('--mhtBin', dest='mhtBin', default=0, type=int, action='store', help='MHT bin')
parser.add_argument('--syst', dest='syst', default="nominal", type=str, action='store', help='Systematic')
parser.add_argument('-c','--ctau', dest='ctau', default='1', type=str, action='store', help='Ctau value')
parser.add_argument('-o','--output', dest='output', default='.', type=str, action='store', help='Output dir')
#parser.add_argument('-m','--mass', dest='mass', default=[], type=str, action='append', help='LLP:LSP mass pairs')

args = parser.parse_args()


filePath = "/vols/cms/mkomm/LLP/NANOX_SR/"

#lumi = "35.822" #16 only
lumi = "77.351" #16+17

with open('eventyields.json',) as f:
    genweights = json.load(f)

globalDataWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/nominal_met<1.25)*(nominal_minPhi>0.2)"
if args.syst=="nominal":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/nominal_met<1.25)*(nominal_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="jesUp":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_jesTotalUp==0)*(nselectedJets_jesTotalUp>2)*(jesTotalUp_mht>300.)*(jesTotalUp_mht/jesTotalUp_met<1.25)*(jesTotalUp_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="jesDown":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_jesTotalDown==0)*(nselectedJets_jesTotalDown>2)*(jesTotalDown_mht>300.)*(jesTotalDown_mht/jesTotalDown_met<1.25)*(jesTotalDown_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="jerUp":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_jerUp==0)*(nselectedJets_jerUp>2)*(jerUp_mht>300.)*(jerUp_mht/jerUp_met<1.25)*(jerUp_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="jerDown":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_jerDown==0)*(nselectedJets_jerDown>2)*(jerDown_mht>300.)*(jerDown_mht/jerDown_met<1.25)*(jerDown_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="unclEnUp":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/unclEnUp_met<1.25)*(nominal_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="unclEnDown":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/unclEnDown_met<1.25)*(nominal_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweight)"
elif args.syst=="puUp":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/nominal_met<1.25)*(nominal_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweightUp)"
elif args.syst=="puDown":
    globalMCWeight = "(nvetoMuons==0)*(nvetoElectrons==0)*(nvetoFwdJets_nominal==0)*(nselectedJets_nominal>2)*(nominal_mht>300.)*(nominal_mht/nominal_met<1.25)*(nominal_minPhi>0.2)*(genweight)*"+lumi+"*1000.0*(puweightDown)"




signalConfigs = []
massesDict = {}
ctauSampleName = "SMS-T1qqqq_ctau-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"%args.ctau
for signalSample in [ctauSampleName,ctauSampleName+"_extra"]:
    for llpMass in genweights[signalSample]:
        for lspMass in genweights[signalSample][llpMass]:
            if genweights[signalSample][llpMass][lspMass]<10:
                continue
            if int(llpMass)==int(lspMass):
                continue
            if not massesDict.has_key(llpMass):
                massesDict[llpMass] = {}
            if not massesDict[llpMass].has_key(lspMass):
                massesDict[llpMass][lspMass] = 0.
            massesDict[llpMass][lspMass]+=genweights[signalSample][llpMass][lspMass]

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
            "weight":"(%7.5e)*(llp==%i)*(lsp==%i)*%s"%(genweight,int(llpMass),int(lspMass),globalMCWeight)
        })
    print
        

xsecs = {

    # Drell-Yan + jets
    # NLO
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":147.40,
    "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":40.99,
    "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":5.678,
    "DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":1.367,
    "DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.6304,
    "DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.1514,
    "DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8":0.003565,
    # "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":18610,
    # "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":1921.8*3,
    
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
    "WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1319*1.21, 
    "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1345*1.21,
    "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 359.7*1.21,
    "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 48.91*1.21,
    "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 12.05*1.21,
    "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 5.501*1.21,
    "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 1.329*1.21,
    "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8": 0.03216*1.21 ,
    
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
    
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    # NNLO
    "WW_TuneCUETP8M1_13TeV-pythia8": 118.7,
    
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    # NLO
    "WZ_TuneCUETP8M1_13TeV-pythia8": 47.13,
    "ZZ_TuneCUETP8M1_13TeV-pythia8": 16.523,
}


mcSetDict = {
    "diboson": {
        "processes": [
            "WW_TuneCUETP8M1_13TeV-pythia8",
            "WZ_TuneCUETP8M1_13TeV-pythia8",
            "ZZ_TuneCUETP8M1_13TeV-pythia8"
        ],
        "weight":globalMCWeight,
    },
    "WJets": {
        "processes": [
            "WToLNu_0J_13TeV-amcatnloFXFX-pythia8",
            "WToLNu_1J_13TeV-amcatnloFXFX-pythia8",
            "WToLNu_2J_13TeV-amcatnloFXFX-pythia8"
        ],
        "weight":globalMCWeight,
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
        "weight":globalMCWeight,
    },
    "ZJets": {
        "processes": [
            "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"
        ],
        "weight":globalMCWeight,
    },
    "ttbar": {
        "processes": [
            "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen",
        ],
        "weight":globalMCWeight,
    },
    "ttbarHT": {
        "processes": [
            "TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
            "TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
        ],
        "weight":globalMCWeight,
    },
    "st":{
        "processes": [
            #"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1",
            "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
            "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
            "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
            "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"
        ],
        "weight":globalMCWeight,
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
            #"QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8",
            #"QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8",
            #"QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8",
            #"QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8",
            #"QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8",
        ],
        "weight":globalMCWeight,
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
        "weight":globalMCWeight,
    },
    
    "ZNuNu": {
        "processes": [
            "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
            "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
            "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
            "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
            "ZJetsToNuNu_HT-800to1200_13TeV-madgraph",
            "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
            "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph",
        ],
        "weight":globalMCWeight,
        
    },
}


    

processDict = {}

for processFolder in os.listdir(filePath):
    
    print "reading ",processFolder,"...",
    for f in os.listdir(os.path.join(filePath,processFolder)):
        if not f.endswith(".root"):
            continue
        fullFilePath = os.path.join(filePath,processFolder,f)
        if not processDict.has_key(processFolder):
            processDict[processFolder] = {
                "files":[],
                "weight":"noweightfound",
            }
        processDict[processFolder]["files"].append(fullFilePath)
        
    if processDict.has_key(processFolder):
        if xsecs.has_key(processFolder):
            processDict[processFolder]["weight"] = "("+str(1.*xsecs[processFolder]/genweights[processFolder])+")"
            print processDict[processFolder]["weight"]
        else: 
            print "no xsec available"
    else:
        print "no files found"



binningScheme = []


'''
htBins = [200,700,1300,10000]
mhtBins = [300,600,10000]
jetBins = [3,4,5,50]

if args.htBin>=(len(htBins)-1):
    print "Error - htbin out of range"
    sys.exit(1)
if args.jetBin>=(len(jetBins)-1):
    print "Error - htbin out of range"
    sys.exit(1)
'''
if args.syst=="nominal":
    htVar = "nominal_ht"
    mhtVar = "nominal_mht"
    njetVar = "nselectedJets_nominal"
    taggerSyst = "nominal"
elif args.syst=="jesUp":
    htVar = "jesTotalUp_ht"
    mhtVar = "jesTotalUp_mht"
    njetVar = "nselectedJets_jesTotalUp"
    taggerSyst = "jesUp"
elif args.syst=="jesDown":
    htVar = "jesTotalDown_ht"
    mhtVar = "jesTotalDown_mht"
    njetVar = "nselectedJets_jesTotalDown"
    taggerSyst = "jesDown"
elif args.syst=="jerUp":
    htVar = "jerUp_ht"
    mhtVar = "jerUp_mht"
    njetVar = "nselectedJets_jerUp"
    taggerSyst = "jerUp"
elif args.syst=="jerDown":
    htVar = "jerDown_ht"
    mhtVar = "jerDown_mht"
    njetVar = "nselectedJets_jerDown"
    taggerSyst = "jerDown"
elif args.syst=="unclEnUp":
    htVar = "nominal_ht"
    mhtVar = "nominal_mht"
    njetVar = "nselectedJets_nominal"
    taggerSyst = "nominal"
elif args.syst=="unclEnDown":
    htVar = "nominal_ht"
    mhtVar = "nominal_mht"
    njetVar = "nselectedJets_nominal"
    taggerSyst = "nominal"
elif args.syst=="puUp":
    htVar = "nominal_ht"
    mhtVar = "nominal_mht"
    njetVar = "nselectedJets_nominal"
    taggerSyst = "nominal"
elif args.syst=="puDown":
    htVar = "nominal_ht"
    mhtVar = "nominal_mht"
    njetVar = "nselectedJets_nominal"
    taggerSyst = "nominal"

ctauHack = args.ctau
if args.ctau=="1":
    ctauHack = "0"



histBinning= {
    "name":"llp",#%(jetBins[args.jetBin],htBins[args.htBin],mhtBins[args.mhtBin]),
    "var":"((llpdnnx_da_%s_%s_LLP_min0)>0.98)*1"%(taggerSyst,ctauHack)+\
          "+((llpdnnx_da_%s_%s_LLP_min0)>0.7)*2"%(taggerSyst,ctauHack)+\
          "+((llpdnnx_da_%s_%s_LLP_min1)>0.8)*4"%(taggerSyst,ctauHack)+\
          "+((llpdnnx_da_%s_%s_LLP_min1)>0.5)*8"%(taggerSyst,ctauHack),
    "weight":"(%s>=%i)*(%s<%i)*(%s>=%.1f)*(%s<%.1f)*(%s>=%.1f)*(%s<%.1f)"%(
        njetVar,jetBins[args.jetBin],njetVar,jetBins[args.jetBin+1],
        htVar,htBins[args.htBin],htVar,htBins[args.htBin+1],
        mhtVar,mhtBins[args.mhtBin],mhtVar,mhtBins[args.mhtBin+1]
    ),
    "binning":numpy.linspace(-0.5,3.5,num=5)
}


print histBinning
    
histograms = []

def removeNegEntries(hist):
    for ibin in range(hist.GetNbinsX()+2):
        c = hist.GetBinContent(ibin)
        if c<10**-6:
            hist.SetBinContent(ibin,10**-6)
            if hist.GetBinError(ibin)<10**-5:
                hist.SetBinError(ibin,10**-5)
                
systPostfix = "_"+args.syst
if args.syst=="nominal":
    systPostfix=""
    
histBinningWeight = histBinning["weight"]
print "bin: ",histBinning["name"]
for backgroundMCSet in ["WJetsHT","WJets","st","ttbar","ttbarHT","QCDHT","QCD","ZNuNu"]:
    print "  background: ",backgroundMCSet
    setWeight = mcSetDict[backgroundMCSet]["weight"]
    histBackground = ROOT.TH1F(histBinning["name"]+"_"+backgroundMCSet+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
    histBackground.SetDirectory(0)
    histograms.append(histBackground)
    histBackground.Sumw2()
    for process in mcSetDict[backgroundMCSet]["processes"]:
        processWeight = processDict[process]["weight"]
        #print histBinningWeight+"*"+setWeight+"*"+processWeight
        for ifile,f in enumerate(processDict[process]["files"]):
            rootFile = ROOT.TFile.Open(f)
            tree = rootFile.Get("Friends")
            if tree:
                processHist = ROOT.TH1F(histBinning["name"]+"_"+backgroundMCSet+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                tree.Project(processHist.GetName(),histBinning["var"],histBinningWeight+"*"+setWeight+"*"+processWeight)
                
                processHist.SetDirectory(0)
                histBackground.Add(processHist)
            rootFile.Close()
            #break
    removeNegEntries(histBackground)
    print "    -> ",histBackground.GetEntries(),"/",histBackground.Integral()," entries/integral"
            
for signalCfg in signalConfigs:
    print "  signal: ",signalCfg["name"]
    signalHist = ROOT.TH1F(histBinning["name"]+"_"+signalCfg["name"]+systPostfix,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
    signalHist.SetDirectory(0)
    histograms.append(signalHist)
    signalWeight = signalCfg["weight"]
    #print histBinningWeight+"*"+signalWeight
    for process in signalCfg["processes"]:
        for ifile,f in enumerate(processDict[process]["files"]):
            rootFile = ROOT.TFile.Open(f)
            tree = rootFile.Get("Friends")
            if tree:
                processHist = ROOT.TH1F(histBinning["name"]+"_"+signalCfg["name"]+"_"+process+str(random.random())+str(ifile),"",len(histBinning["binning"])-1,histBinning["binning"])
                tree.Project(processHist.GetName(),histBinning["var"],histBinningWeight+"*"+signalWeight)
                processHist.SetDirectory(0)
                signalHist.Add(processHist)
            rootFile.Close()
            #break
    removeNegEntries(signalHist)
    print "    -> ",signalHist.GetEntries(),"/",signalHist.Integral()," entries/integral"
    
outputFile = ROOT.TFile(os.path.join(args.output,"hist_%s_%s_%s.root"%(args.ctau,histBinning["name"],args.syst)),"RECREATE")
for hist in histograms:
    hist.SetDirectory(outputFile)
    hist.Write()
outputFile.Write()
outputFile.Close()


