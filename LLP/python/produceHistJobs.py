import os
import sys

basePath = '/vols/build/cms/mkomm/LLP/CMSSW_8_1_0/src/hists_new'
'''
htBins = [200,700,1300,10000]
mhtBins = [300,600,10000]
jetBins = [3,4,5,50]
'''
ctauValues = ["0p001","0p01","0p1","1","10","100","1000","10000"]
systematics = [ 
    "nominal","jesUp","jesDown","jerUp","jerDown","unclEnUp","unclEnDown","puUp","puDown",
    "wjetsScaleUp","wjetsScaleDown","ttbarScaleUp","ttbarScaleDown","stScaleUp","stScaleDown",
    "znunuScaleUp","znunuScaleDown"
]

jobArrayCfg = []
for ctau in ctauValues:
    for syst in systematics:
        jobArrayCfg.append([
            "-o %s"%basePath,
            "--syst %s"%syst,
            "--ctau %s"%ctau
        ])
            
submitFile = open("runHist.sh","w")
submitFile.write('''#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=8:00:00 
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
python CombineHarvester/LLP/python/driver_new.py ${JOBS[$SGE_TASK_ID]}
date
''')
submitFile.close()


