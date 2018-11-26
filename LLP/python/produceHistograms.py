import ROOT
import numpy
import os
import sys
import json
import random

filePath = "/vols/cms/mkomm/LLP/NANOX_180618_SR/"


with open('eventyields.json',) as f:
    genweights = json.load(f)

massesDict = {
    600:[0,200,400,500],
    800:[0,200,400,600,700],
    1000:[0,200,400,600,800,900],
    1200:[0,200,400,600,800,900,1000,1100],
    1400:[0,200,400,600,800,900,1000,1200,1300],
    1600:[0,200,400,600,800,900,1000,1200,1400],
    1800:[0,200,400,600,800,1000,1200,1400],
    2000:[0,200,400,600,800,1000,1200,1400],
    2200:[0,200,400,600,800,1000,1200,1400],
    2400:[0,200,400,600,800,1000,1200,1400]
}

globalWeight = "(nLooseMuons==0)*(nLooseElectrons==0)*(nSelectedJets_central>2)*(ht_central>200.0)*(met>200.0)"
globalMCWeight = globalWeight+"*(genweight)*35.822*1000.0*(puweight)"
globalDataWeight = globalWeight


signalConfigs = []

for llpMass in sorted(massesDict.keys()):
    for lspMass in massesDict[llpMass]:
        genweight = 1./float(genweights["SMS-T1qqqq_ctau-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"][str(llpMass)][str(lspMass)])
        signalConfigs.append({
            "name":"ctau1_llp%i_lsp%i"%(llpMass,lspMass),
            "llp":llpMass,
            "lsp":lspMass,
            "processes":["SMS-T1qqqq_ctau-1_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"],
            "weight":"(%7.5e)*(llpmass==%i)*(lspmass==%i)*%s"%(genweight,llpMass,lspMass,globalMCWeight)
        })



 

xsecs = {
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":18610,
    "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8":1921.8*3,
    
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

    
    #https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec#Single_top_t_channel_cross_secti
    "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1": 10.32,
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1": 80.95,
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1": 136.02,
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1": 71.7/2.,
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1": 71.7/2.,
    
    #https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO (mtop=172.5 GeV)
    "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen": 831.76,
    
    #https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=WToLNu_*J_13TeV-amcatnloFXFX-pythia8
    "WToLNu_0J_13TeV-amcatnloFXFX-pythia8": 49670.,
    "WToLNu_1J_13TeV-amcatnloFXFX-pythia8": 8264.,
    "WToLNu_2J_13TeV-amcatnloFXFX-pythia8": 3226.,
    
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    "WW_TuneCUETP8M1_13TeV-pythia8": 118.7,
    
    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
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
    "st":{
        "processes": [
            "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1",
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
                "weight":"1",
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

htBins = [200,400,600,800,1000,1400,1600,2000,2400,10000]
for iht in range(len(htBins)-1):
    binningScheme.append(
        {
            "name":"ht"+str(htBins[iht]),
            "var":"daref_max_fromLLP_ctau1",
            "weight":"(ht_central>%.1f)*(ht_central<%.1f)"%(htBins[iht],htBins[iht+1]),
            "binning":numpy.array([0.,0.075,0.15,0.25,0.4,0.6,0.95,1.0])#numpy.linspace(0,1,num=11)
        }
    )
    
histograms = []

def removeNegEntries(hist):
    for ibin in range(hist.GetNbinsX()+2):
        c = hist.GetBinContent(ibin)
        if c<10**-6:
            hist.SetBinContent(ibin,10**-6)
            if hist.GetBinError(ibin)<10**-5:
                hist.SetBinError(ibin,10**-5)
    
for histBinning in binningScheme:
    histBinningWeight = histBinning["weight"]
    print "bin: ",histBinning["name"]
    for backgroundMCSet in ["WJets","ZJets","st","ttbar","QCD"]:
        print "  background: ",backgroundMCSet
        setWeight = mcSetDict[backgroundMCSet]["weight"]
        histBackground = ROOT.TH1F(histBinning["name"]+"_"+backgroundMCSet,"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
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
        signalHist = ROOT.TH1F(histBinning["name"]+"_"+signalCfg["name"],"",len(histBinning["binning"])-1,histBinning["binning"]) #<channel>_<process>_<syst>
        signalHist.SetDirectory(0)
        histograms.append(signalHist)
        signalWeight = signalCfg["weight"]
        #print signalWeight
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
        
outputFile = ROOT.TFile("hist.root","RECREATE")
for hist in histograms:
    hist.SetDirectory(outputFile)
    hist.Write()
outputFile.Write()
outputFile.Close()


