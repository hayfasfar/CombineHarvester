import os
import sys
import json
import ROOT
import math
import numpy
import scipy.interpolate
import CombineHarvester.CombineTools.ch as ch

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

def getHist(fileName,histName):
    rootFile = ROOT.TFile(fileName)
    hist = rootFile.Get(histName)
    hist = hist.Clone()
    hist.SetDirectory(0)
    rootFile.Close()
    return hist
    
alpha = 1. - 0.6827
upErr = ROOT.Math.gamma_quantile_c(alpha/2,1,1)
    
def makeDatacard(cats,ctau,signalProc,histPath,outputPath,llpMCEff=1.,systematics=[]):
    print "Producing datacards for signal '"+signalProc+"' under '"+outputPath+"'"
    if os.path.exists(outputPath):
        pass
    else:
        os.makedirs(outputPath)
        
        
    shapeHist = getHist(
        os.path.join(histPath,"hists_%s.root"%ctau),
        "llpD_QCDHT"
    )
    shapeHist.Scale(0.)
    nbins = shapeHist.GetNbinsX()
    binMin = -0.5
    binMax = nbins-0.5
        
    cb = ch.CombineHarvester()
    #cb.SetVerbosity(0)
    
    '''
    cb.ParseDatacard(os.path.join(os.path.dirname(os.path.abspath(__file__)),"dummy.txt"),
		"",
		"",
		"",
		0,
		""
	) 	
    '''
    cb.AddExtArgValue('llpEff', 1.0)
    cb.GetParameter('llpEff').set_range(0.75, 1.25)
    
    '''
    llpEffFakeProc = ch.Process()
    llpEffFakeProc.set_process("llpDummyProc")
    llpEffFakeProc.set_bin('llpA')
    llpEffFakeProc.set_era('13TeV2016')
    hist = ROOT.TH1F("llpDummyHist","",nbins,binMin,binMax)
    hist.SetDirectory(0)
    llpEffFakeProc.set_shape(hist,True)
    cb.InsertProcess(llpEffFakeProc)
    cb.cp().process(['llpDummyProc']).AddSyst(cb,"llpEff","lnN",ch.SystMap("era")(["13TeV2016"],1.5))
    cb.PrintAll()
    '''
    
    cb.AddProcesses(era=["13TeV2016"],procs=["WJets","st","ttbar","ZNuNu"],bin=cats.values(),signal=False)

    
    #add signal only to region 'D'
    for nLLP in range(7):
        for nLLPTagged in range(nLLP+1):
            processSubName = signalProc+("_%illp%it"%(nLLP,nLLPTagged))
            cb.AddProcesses(era=["13TeV2016"],procs=[processSubName],bin=[cats['D']],signal=True)
            cb.cp().process([processSubName]).AddSyst(cb,processSubName+"_rate","rateParam",
                ch.SystMap("era")(["13TeV2016"],(
                    #"1.0",
                    "TMath::Power(TMath::Range(0.01,0.99,@0*%.3f),%i)*TMath::Power(1-TMath::Range(0.01,0.99,@0*%.3f),%i)/(TMath::Power(%.3f,%i)*TMath::Power(1-%.3f,%i))"%(
                        llpMCEff,nLLPTagged,llpMCEff,nLLP-nLLPTagged,
                        llpMCEff,nLLPTagged,llpMCEff,nLLP-nLLPTagged 
                    ),
                    "llpEff"
                ))
            )
    

    cb.cp().AddSyst(cb,"lumi","lnN",ch.SystMap("era")(["13TeV2016"],1.026))
    #cb.cp().process(['QCDHT']).AddSyst(cb,"qcd_$ERA","lnN",ch.SystMap("era")(["13TeV2016"],1.2))
    cb.cp().process(['WJets']).AddSyst(cb,"wzjets_yield","lnN",ch.SystMap("era")(["13TeV2016"],1.2))
    cb.cp().process(['st','ttbar']).AddSyst(cb,"topbkg_yield","lnN",ch.SystMap("era")(["13TeV2016"],1.2))
    cb.cp().process(['ZNuNu']).AddSyst(cb,"znunnu_yield","lnN",ch.SystMap("era")(["13TeV2016"],1.2))
    
    #
        
    for syst in systematics:
        cb.cp().AddSyst(cb,syst, "shape", ch.SystMap("era")(["13TeV2016"],1.0))
        
    cb.cp().ExtractShapes(
           os.path.join(histPath,"hists_%s.root"%ctau),
          "$BIN_$PROCESS",
          "$BIN_$PROCESS_$SYSTEMATIC")
          
    bbFactory = ch.BinByBinFactory()
    bbFactory.SetAddThreshold(0.1)
    #bbFactory.SetMergeThreshold(0.5)
    bbFactory.SetFixNorm(True)
    #bbFactory.SetPoissonErrors(True)
    bbFactory.SetPattern("bb_$BIN_$PROCESS_bin_$#")
    #bbFactory.MergeBinErrors(cb.cp().backgrounds())
    bbFactory.AddBinByBin(cb.cp().process(['WJets','st','ttbar','ZNuNu']), cb)
          
    #required for toys to know it's not unbinned
    
    dummyObs = []
    for cat in cats.keys():
        obs = ch.Observation()
        histObs = getHist(
            os.path.join(histPath,"hists_%s.root"%ctau),
            "llp%s_sum"%cat
        )
        histObs.SetName("llp%s_data_obs"%cats[cat][0])
        histObs.SetDirectory(0)
        #obs.set_name("llp%s_data_obs"%cats[cat][0])
        obs.set_shape(histObs,True)
        obs.set_bin(cats[cat][1])
        obs.set_era('13TeV2016')
        cb.InsertObservation(obs)
        
    #cb.AddObservations(['*'], ['*'], ['13TeV2016'], ['*'], cats.values())
    
    for region in ["A","B","C"]:
        qcdHistNominal = getHist(
            os.path.join(histPath,"hists_%s.root"%ctau),
            "llp%s_QCDHT"%region
        )
        sumHistNominal = getHist(
            os.path.join(histPath,"hists_%s.root"%ctau),
            "llp%s_sum"%region
        )
        for ibin in range(nbins):
            proc = ch.Process()
            processName = 'QCD_llp%s_bin%i'%(region,ibin+1)
            systName = 'qcd_llp%s_bin%i'%(region,ibin+1)
            proc.set_process(processName)
            proc.set_bin('llp%s'%region)
            proc.set_era('13TeV2016')
            hist = ROOT.TH1F("qcdHist_llp%s_bin%i"%(region,ibin+1),"",nbins,binMin,binMax)
            hist.SetBinContent(ibin+1,1)
            hist.SetBinError(ibin+1,0.00000000001)
            hist.SetDirectory(0)
            proc.set_shape(hist,True)
            cb.InsertProcess(proc)
            
            cb.cp().process([processName]).AddSyst(cb,systName,"rateParam",
                ch.SystMap("era")(["13TeV2016"],max(10**-3,qcdHistNominal.GetBinContent(ibin+1)))
            )
            param =  cb.GetParameter(systName)
            avgWeight = qcdHistNominal.Integral()/qcdHistNominal.GetEntries() if hist.GetEntries()>0 else 1.
            content = max(0.,qcdHistNominal.GetBinContent(ibin+1))
            err = max(0.,math.sqrt(sumHistNominal.GetBinContent(ibin+1)))
            param.set_val(content)
            param.set_range(0.,content+2.*err)
    
    qcdProcessNamesInSR = []
    
    for ibin in range(nbins):
        proc = ch.Process()
        processName = 'QCD_llp%s_bin%i'%('D',ibin+1)
        qcdProcessNamesInSR.append(processName)
        systNameA = 'qcd_llp%s_bin%i'%('A',ibin+1)
        systNameB = 'qcd_llp%s_bin%i'%('B',ibin+1)
        systNameC = 'qcd_llp%s_bin%i'%('C',ibin+1)
        systNameD = 'qcd_llp%s_bin%i'%('D',ibin+1)
        proc.set_process(processName)
        proc.set_bin('llpD')
        proc.set_era('13TeV2016')
        hist = ROOT.TH1F("qcdHist_llp%s_bin%i"%('D',ibin+1),"",nbins,binMin,binMax)
        hist.SetBinContent(ibin+1,1.)
        hist.SetBinError(ibin+1,0.00000000001)
        hist.SetDirectory(0)
        proc.set_shape(hist,True)
        cb.InsertProcess(proc)
        # B | D   
        # A | C   => C/A=D/B => D = C*B/A
        cb.cp().process([processName]).AddSyst(cb,systNameD,"rateParam",
            ch.SystMap("era")(["13TeV2016"],(
                "TMath::Max(@0,0.000001)*TMath::Range(0,10.,TMath::Max(@1,0.000001)/TMath::Max(@2,0.000001))",
                systNameB+","+systNameC+","+systNameA
            ))
        )
    
    cb.cp().process(qcdProcessNamesInSR).AddSyst(cb,'qcd_yield',"lnN",ch.SystMap("era")(["13TeV2016"],1.5))
    
    
    
    
    
    
    cb.cp().WriteDatacard(
        os.path.join(outputPath,"out.txt"),
        os.path.join(outputPath,"out.root")
    )
    

ctauValues = ["0p001","0p01","0p1","1","10","100","1000","10000"]
#ctauValues = ["1"]

systematics = ["jes","jer","unclEn","pu","wjetsScale","ttbarScale","stScale","znunuScale"]
#systematics = ["jes","jer","unclEn","pu","wjetsScale","ttbarScale","znunuScale"]

categories = {}
for region in ['A','B','C','D']:
    categories[region] = ((len(categories),"llp"+region))
print categories

with open('eventyields.json',) as f:
    genweights = json.load(f)
    
with open('llpEfficiency.json',) as f:
    llpEfficiency = json.load(f)
    
    
def getTheoryXsecFct(filePath):
    f = open(filePath)
    llpvalues = []
    xsecvalues = []
    xsecvaluesUp = []
    xsecvaluesDown = []
    for l in f:
        if len(l)==0:
            continue
        splitted = l.split(",")
        if len(splitted)!=3:
            print "cannot parse xsec line: ",l
            continue
        llpvalues.append(float(splitted[0]))
        xsec = float(splitted[1])
        xsecRelUnc = float(splitted[2])
        xsecvalues.append(math.log(xsec))
        xsecvaluesUp.append(math.log(xsec*(1+xsecRelUnc)))
        xsecvaluesDown.append(math.log(xsec*(1-xsecRelUnc)))
        
    llpvalues = numpy.array(llpvalues,dtype=numpy.float32)
    xsecvalues = numpy.array(xsecvalues,dtype=numpy.float32)
    xsecvaluesUp = numpy.array(xsecvaluesUp,dtype=numpy.float32)
    xsecvaluesDown = numpy.array(xsecvaluesDown,dtype=numpy.float32)
    
    tckXsec = scipy.interpolate.splrep(llpvalues,xsecvalues,s=1e-3)
    tckXsecUp = scipy.interpolate.splrep(llpvalues,xsecvaluesUp,s=1e-3)
    tckXsecDown = scipy.interpolate.splrep(llpvalues,xsecvaluesDown,s=1e-3)

    def getValue(llpMass):
        xsec = math.exp(scipy.interpolate.splev(llpMass,tckXsec))
        xsecUp = math.exp(scipy.interpolate.splev(llpMass,tckXsecUp))
        xsecDown = math.exp(scipy.interpolate.splev(llpMass,tckXsecDown))
        return xsec,xsecUp,xsecDown
    return getValue
    
theoXsec = getTheoryXsecFct("theory_xsec.dat")

massesDict = {}
for ctau in ctauValues:
    ctauSampleName = "SMS-T1qqqq_ctau-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"%ctau
    if not massesDict.has_key(ctau):
        massesDict[ctau] = {}
    for signalSample in [ctauSampleName,ctauSampleName+"_extra"]:
        for llpMass in genweights[signalSample].keys():
            for lspMass in genweights[signalSample][llpMass].keys():
                if genweights[signalSample][llpMass][lspMass]["sum"]<10:
                    continue
                if int(llpMass)==int(lspMass):
                    continue
                if not massesDict[ctau].has_key(llpMass):
                    massesDict[ctau][llpMass] = []
                if not lspMass in massesDict[ctau][llpMass]:
                    massesDict[ctau][llpMass].append(lspMass)
                

basePath = "cards_new"
histPath = "hists_new"
if os.path.exists(os.path.join(basePath,'log')):
    pass
else:
    os.makedirs(os.path.join(basePath,'log'))

jobLimitHybridArrayCfg = []
jobLimitArrayCfg = []
jobToysArrayCfg = []
jobImpactArrayCfg = []

for ctau in ctauValues:
    for llpMass in sorted(massesDict[ctau].keys()):
        xsec,xsecUp,xsecDown = theoXsec(1.*int(llpMass))
        for lspMass in massesDict[ctau][llpMass]:
            
            signalProcess = "ctau%s_llp%s_lsp%s"%(str(ctau),str(llpMass),str(lspMass))
            datacardPath = os.path.join(basePath,signalProcess)
            ctauHack = ctau
            if ctau =="1":
                ctauHack = "0"
            selectedJets = llpEfficiency[ctauHack][llpMass][lspMass]
            eff = 1.*selectedJets["tagged"]/selectedJets["total"] if selectedJets["total"]>0. else 1.
            '''
            makeDatacard(
                categories,
                ctau,
                signalProcess,
                histPath,
                datacardPath,
                systematics=systematics,
                llpMCEff=eff
            )
            '''
            jobLimitArrayCfg.append({
                "path":datacardPath,
                "cmd": [
                    "combine -M AsymptoticLimits --setParameterRanges llpEff=0.75,1.25 --rAbsAcc 0.000001 --run expected --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=9999999 --cminPreScan -t -1 -d out.txt"
                ]
            })
            
            jobLimitHybridArrayCfg.append({
                "path":datacardPath,
                "cmd": [
                    "combine -M HybridNew -t 1 --frequentist --testStat LHC --singlePoint 0.01 --rule CLs --rAbsAcc 0.000001 --saveToys --saveHybridResult --toysH 100 --setParameterRanges llpEff=0.75,1.25 --cminDefaultMinimizerStrategy 0 -d out.txt"
                ]
            })
            
            for i in range(5):
                jobToysArrayCfg.append({
                    "path":datacardPath,
                    "cmd": [
                        "combine -M Significance --verbose -1 -t 1000 -s 123%02i45%02i --rMax 1000 --setParameterRanges llpEff=0.75,1.25 --maxTries 10 --cminDefaultMinimizerStrategy 0 --expectSignal=0 -d out.txt"%(i*21+1,i*7+3)
                    ]
                })
            
            
            jobImpactArrayCfg.append({
                "path":datacardPath,
                "cmd": [
                    "text2workspace.py -m 120 out.txt -o workspace.root",
                    "combineTool.py -M Impacts -d workspace.root -m 120 --setParameterRanges llpEff=0.5,1.5 --X-rtd MINIMIZER_analytic --robustFit 1 --doInitialFit -t -1 --expectSignal=%6.4e"%(xsec),
                    "combineTool.py -M Impacts -d workspace.root -m 120 --setParameterRanges llpEff=0.5,1.5 --X-rtd MINIMIZER_analytic --robustFit 1 --doFits --exclude 'bb_.*,qcd_llp[ABC].*' -t -1 --expectSignal=%6.4e"%(xsec),
                    "combineTool.py -M Impacts -d workspace.root -m 120 --setParameterRanges llpEff=0.5,1.5 --X-rtd MINIMIZER_analytic --exclude 'bb_.*,qcd_llp[ABC].*' -o impacts.json",
                    "plotImpacts.py -i impacts.json -o impacts"
                ]
            })
            
            #make a fit with 0 signal
            #combine -M MultiDimFit -t -1 --expectSignal=0 --robustFit 1 -d cards_new/ctau1_llp1000_lsp0/out.txt --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=99999999999 --saveFitResult
            #PostFitShapesFromWorkspace -w cards_new/ctau1_llp1000_lsp0/workspace.root -o postfit.root --print 1 -f multidimfit.root:fit_mdf --sampling 1 --postfit 1
            
            #break
        #break

def makeSubmitFile(jobArrayCfg,name):
    if len(jobArrayCfg)==0:
        print "No jobs for '"+name+"' -> skip"
        return
    submitFile = open(name,"w")#"runCombine.sh","w")
    submitFile.write('''#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=01:30:00 
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
        submitFile.write("CMDPATH="+jobCfg["path"]) 
        submitFile.write("; CMDS=(\\\""+jobCfg["cmd"][0]+"\\\"") 
        for cmd in jobCfg["cmd"][1:]:
            submitFile.write(" \\\""+cmd+"\\\"")
        submitFile.write(")\"\n")
    submitFile.write(")\n")

    submitFile.write('''

eval ${JOBS[$SGE_TASK_ID]}
echo "CMDPATH="$CMDPATH
cd $CMDPATH
ls -lh
for cmd in "${CMDS[@]}"
    do
    echo "-----------------------------------------------------"
    echo "command: "$cmd
    echo "-----------------------------------------------------"
    $cmd
    done


date
    ''')
    submitFile.close()
    
    
    '''
    #echo ${JOBS[$SGE_TASK_ID]}
    #cd ${JOBS[$SGE_TASK_ID]}
    #combine -M FitDiagnostics --plots --saveWithUncertainties -t -1 -d out.txt
    #combine -M MultiDimFit --saveFitResult -t -1 -d out.txt

    #expected impacts
    #note: excluding regex only works in hacked combineTool.py
    #text2workspace.py -m 120 out.txt -o workspace.root
    #combineTool.py -M Impacts -d workspace.root -m 120 --X-rtd MINIMIZER_analytic --robustFit 1 --doInitialFit -t -1 --expectSignal=1 
    #combineTool.py -M Impacts -d workspace.root -m 120 --X-rtd MINIMIZER_analytic --robustFit 1 --doFits --exclude "bb_.*,qcd_llp[ABC].*" -t -1 --expectSignal=1
    #combineTool.py -M Impacts -d workspace.root -m 120 --X-rtd MINIMIZER_analytic --exclude "bb_.*,qcd_llp[ABC].*" -o impacts.json
    #plotImpacts.py -i impacts.json -o impacts

    #expected limits
    #combine -M AsymptoticLimits --setParameterRanges llpEff=0.75,1.25 --rAbsAcc 0.000001 --run expected --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=99999999999 -t -1 -d out.txt
    '''
    
makeSubmitFile(jobLimitHybridArrayCfg,"runHybridLimits.sh")
makeSubmitFile(jobLimitArrayCfg,"runLimits.sh")
makeSubmitFile(jobToysArrayCfg,"runToys.sh")
makeSubmitFile(jobImpactArrayCfg,"runImpacts.sh")

