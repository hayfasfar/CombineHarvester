import os
import sys
import json
import CombineHarvester.CombineTools.ch as ch

    
def makeDatacard(cats,ctau,signalProc,outputPath,systematics=[]):
    print "Producing datacards for signal '"+signalProc+"' under '"+outputPath+"'"
    if os.path.exists(outputPath):
        pass
    else:
        os.makedirs(outputPath)
        
    cb = ch.CombineHarvester()
    cb.AddProcesses(era=["13TeV2016"],procs=["WJetsHT","st","ttbar","QCDHT","ZNuNu"],bin=cats,signal=False)
    cb.AddProcesses(era=["13TeV2016"],procs=[signalProc],bin=cats,signal=True)

    cb.cp().AddSyst(cb,"lumi_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.026))
    cb.cp().process(['QCDHT']).AddSyst(cb,"qcd_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.3))
    cb.cp().process(['WJetsHT']).AddSyst(cb,"wzjets_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.3))
    cb.cp().process(['st','ttbar']).AddSyst(cb,"topbkg_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.1))
    cb.cp().process(['ZNuNu']).AddSyst(cb,"ZNuNu_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.3))
    
    for syst in systematics:
        cb.cp().AddSyst(cb,syst, "shape", ch.SystMap("era")(["13TeV2016"],1.0))
        
    cb.cp().ExtractShapes(
           "hists2/hist_%s.root"%ctau,
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
    

htBins = [200,700,1300,10000]
mhtBins = [300,600,10000]
jetBins = [3,4,5,50]
ctauValues = ["0p001","0p01","0p1","1","10","100","1000","10000"]
systematics = ["jes","jer","unclEn","pu"]

categories = []
for jetBin in range(len(jetBins)-1):
    for htBin in range(len(htBins)-1):
        for mhtBin in range(len(mhtBins)-1):
            categories.append((len(categories),"jet%i_ht%i_mht%i"%(jetBins[jetBin],htBins[htBin],mhtBins[mhtBin])))
print categories

with open('eventyields.json',) as f:
    genweights = json.load(f)

massesDict = {}
for ctau in ctauValues:
    ctauSampleName = "SMS-T1qqqq_ctau-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"%ctau
    if not massesDict.has_key(ctau):
        massesDict[ctau] = {}
    for signalSample in [ctauSampleName,ctauSampleName+"_extra"]:
        for llpMass in genweights[signalSample]:
            for lspMass in genweights[signalSample][llpMass]:
                if genweights[signalSample][llpMass][lspMass]<10:
                    continue
                if int(llpMass)==int(lspMass):
                    continue
                if not massesDict[ctau].has_key(llpMass):
                    massesDict[ctau][llpMass] = []
                if not lspMass in massesDict[ctau][llpMass]:
                    massesDict[ctau][llpMass].append(lspMass)

basePath = "cards3"
if os.path.exists(os.path.join(basePath,'log')):
    pass
else:
    os.makedirs(os.path.join(basePath,'log'))

jobArrayCfg = []
for ctau in ctauValues:
    for llpMass in sorted(massesDict[ctau].keys()):
        for lspMass in massesDict[ctau][llpMass]:
            signalProcess = "ctau%s_llp%s_lsp%s"%(str(ctau),str(llpMass),str(lspMass))
            datacardPath = os.path.join(basePath,signalProcess)
            makeDatacard(
                categories,
                ctau,
                signalProcess,
                datacardPath,
                systematics=systematics
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
#combine -M FitDiagnostics --plots --saveWithUncertainties -t -1 -d out.txt
#combine -M MultiDimFit --saveFitResult -t -1 -d out.txt
combine -M AsymptoticLimits --rAbsAcc 0.0000001 --saveToys -t -1 -d out.txt
date
''')
submitFile.close()


