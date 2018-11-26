import subprocess
import re
import time
import os
import ROOT
import json
import random 
import numpy
basepath = "/vols/build/cms/mkomm/LLP/CMSSW_10_2_0_dev/src/testing"

txtFiles = sorted(os.listdir(basepath))

processDict = {}
pileupHists = {}
massHistDict = {}

for txtFile in txtFiles:
    process = txtFile.split(".")[0]
    
    if process.find("SingleMu")>=0:
        print "skip ",process
        continue
        
    if process.find("SingleEle")>=0:
        print "skip ",process
        continue
    
    if process.find("HTMHT")>=0:
        print "skip ",process
        continue
    
    if process.find("MET")>=0:
        print "skip ",process
        continue
    
    if not processDict.has_key(process):
        pileupHists[process] = ROOT.TH1F(process,"",101,0,100)
        pileupHists[process].SetDirectory(0)
    
    if not processDict.has_key(process):
        if process.find("SMS-T1qqqq_ctau")>=0:
            
            processDict[process] = {}
            
            
            #assume binning in 100 GeV steps
            massHistDict[process] = ROOT.TH2F("mass"+process,"",
                31,
                numpy.linspace(-50.,3050.,32),
                31,
                numpy.linspace(-50.,3050.,32)
            )
        else:
            processDict[process] = 0.
    
    f = open (os.path.join(basepath,txtFile))
    print "reading ",txtFile, "...",
    for l in f:
        if l.startswith("#"):
            continue
        if len(l)>0:
            fileName = l.replace("\n","").replace("\r","")
            print "opening",fileName
            rootFile=None
            while (rootFile==None):
                rootFile = ROOT.TFile.Open(fileName)
                if not rootFile:
                    print "Cannot found file: ",rootFile, "-> retry"
                    continue
                break
            tree = rootFile.Get("Events")
            if not rootFile or not tree:
                print "Cannot found tree in file: ",rootFile
                continue
            hPU = pileupHists[process].Clone(pileupHists[process].GetName()+str(random.random()))
            tree.Project(hPU.GetName(),"Pileup_nTrueInt","genWeight")
            pileupHists[process].Add(hPU)
            
            
            if process.find("SMS-T1qqqq_ctau")>=0:
                massHist = massHistDict[process].Clone(massHistDict[process].GetName()+str(random.random()))
                #note: synatx is 'y:x' so lsp is on x-axis
                tree.Project(massHist.GetName(),"llpinfo_llp_mass[0]:llpinfo_lsp_mass[0]","genWeight*(nllpinfo>0)")
                massHistDict[process].Add(massHist)
            else:
                processDict[process] += hPU.Integral()
            
            rootFile.Close()
            #break
    if process.find("SMS-T1qqqq_ctau")>=0:
        for llpMassBin in range(massHistDict[process].GetYaxis().GetNbins()):
            llpMass = int(massHistDict[process].GetYaxis().GetBinCenter(llpMassBin+1))
            for lspMassBin in range(massHistDict[process].GetXaxis().GetNbins()):
                lspMass = int(massHistDict[process].GetXaxis().GetBinCenter(lspMassBin+1))
                content = massHistDict[process].GetBinContent(lspMassBin+1,llpMassBin+1)
                if content>0.:
                    if not processDict[process].has_key(llpMass):
                        processDict[process][llpMass] = {}
                    processDict[process][llpMass][lspMass] = content
                    
    else:
        print processDict[process]

with open('eventyields.json', 'w') as outfile:
    json.dump(processDict, outfile,ensure_ascii=True,indent=2,sort_keys=True)    
    
    
    

rootFile = ROOT.TFile("pileup.root","RECREATE")
for process in pileupHists.keys():
    pileupHists[process].SetDirectory(rootFile)
    pileupHists[process].SetName(process)
    pileupHists[process].Write()
rootFile.Close()


