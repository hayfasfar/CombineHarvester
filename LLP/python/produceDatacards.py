import os
import sys

import CombineHarvester.CombineTools.ch as ch

    
def makeDatacard(cats,signalProc,outputPath):
    print "Producing datacards for signal '"+signalProc+"' under '"+outputPath+"'"
    if os.path.exists(outputPath):
        pass
    else:
        os.makedirs(outputPath)
        
    cb = ch.CombineHarvester()
    cb.AddProcesses(era=["13TeV2016"],procs=["WJets","ZJets","st","ttbar","QCD"],bin=cats,signal=False)
    cb.AddProcesses(era=["13TeV2016"],procs=[signalProc],bin=cats,signal=True)

    cb.cp().AddSyst(cb,"lumi_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.026))
    cb.cp().process(['QCD']).AddSyst(cb,"qcd_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],2.0))
    cb.cp().process(['WJets','ZJets']).AddSyst(cb,"wzjets_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.3))
    cb.cp().process(['st','ttbar']).AddSyst(cb,"topbkg_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.1))
    
    cb.cp().ExtractShapes(
           "hist.root",
          "$BIN_$PROCESS",
          "$BIN_$PROCESS_$SYSTEMATIC")


    bbFactory = ch.BinByBinFactory()
    bbFactory.SetAddThreshold(0.1)
    #bbFactory.SetMergeThreshold(0.5)
    bbFactory.SetFixNorm(True)
    bbFactory.SetPattern("bb_$BIN_$PROCESS_bin_$#")
    #bbFactory.MergeBinErrors(cb.cp().backgrounds())
    bbFactory.AddBinByBin(cb.cp().backgrounds(), cb)


    #cb.PrintAll()

    cb.cp().WriteDatacard(
        os.path.join(outputPath,"out.txt"),
        os.path.join(outputPath,"out.root")
    )
    

htBins = [200,400,600,800,1000,1400,1600,2000,2400]

categories = []
for i,htBin in enumerate(htBins):
    categories.append((i,"ht%i"%htBin))
print categories

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

basePath = "cards"
if os.path.exists(os.path.join(basePath,'log')):
    pass
else:
    os.makedirs(os.path.join(basePath,'log'))

jobArrayCfg = []
for ctau in [1]:
    for llpMass in sorted(massesDict.keys()):
        for lspMass in massesDict[llpMass]:
            signalProcess = "ctau%i_llp%i_lsp%i"%(ctau,llpMass,lspMass)
            datacardPath = os.path.join(basePath,signalProcess)
            makeDatacard(
                categories,
                signalProcess,
                datacardPath
            )
            
            jobArrayCfg.append([datacardPath])
            
submitFile = open("runCombine.sh","w")
submitFile.write('''#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=00:30:00 
#$ -t 1-'''+str(len(jobArrayCfg))+'''
#$ -e '''+os.path.join(basePath,'log')+'''/log.$TASK_ID.err
#$ -o '''+os.path.join(basePath,'log')+'''/log.$TASK_ID.out
hostname
date
source ~/.bashrc
cd /vols/build/cms/mkomm/LLP/CMSSW_8_1_0/src
eval `scramv1 runtime -sh`
''')

submitFile.write("JOBS=(\n")
submitFile.write(" \"pseudo job\"\n")
for jobCfg in jobArrayCfg:
    submitFile.write(" \"")
    for opt in jobCfg:
        submitFile.write(" "+opt)
    submitFile.write("\"\n")
submitFile.write(")\n")

submitFile.write('''
echo ${JOBS[$SGE_TASK_ID]}
cd ${JOBS[$SGE_TASK_ID]}
combine -M FitDiagnostics --plots --saveWithUncertainties -d out.txt
combine -M AsymptoticLimits --rAbsAcc 0.000001 --saveToys -t -1 -d out.txt
date
''')
submitFile.close()


