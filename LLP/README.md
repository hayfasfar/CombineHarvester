# Setting limits for HNL analysis

## Step 0: Setup cms environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
```

## Step 1: Produce histograms

* Have a look at ```make_hists.py``` script [here](https://github.com/LLPDNNX/histo/blob/master/limits/make_hists.py) to produce inpute histograms.

## Step 2: Produce datacards & run limits

* Running ```python makeCards.py``` produces a datacard per HNL lifetime, mass, and coupling point. The resulting cards will be located in the ```cards``` directory.
* An SGE batch submission file is created as well: ```runCombine.sh```. It is recommended to test a single datacard with ```SGE_TASK_ID=1 ./runCombine.sh``` before submitting ```qsub runCombine.sh```.

## Step 3: Sanity checks (optional)

Perform maximum likelihood fits with an Asimov data set by artificially setting the signal strength to a given value. This value should be matched by the resulting best-fit value for the signal strength:

* ```combine -M FitDiagnostics -t -1 --expectSignal 0 cards/2016/coupling_1/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt```. Expect to see best fit r=0.
* ```combine -M FitDiagnostics -t -1 --expectSignal 1 cards/2016/coupling_1/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt```. Expect to see best fit r=1.

## Step 4: Collect results & plot limits

Parse results into json files to facilitate plotting of limits:
```
for YEAR in 2016 2017 2018 combined; do combineTool.py -M CollectLimits cards/$YEAR/coupling_*/*/*HNL*.root --use-dirs -o jsons/limits_$YEAR.json; done
python plotLimits.py
```

## Step 5: Impact plots

* Optionally, make a rename dictionary using ```python makePUrenameDict.py``` for aesthetic purposes
* Run the following set of commands. The runtime is very long for all three year combination, so consider only running on one year.
* After unblinding, consider replacing ```-t -1``` by ```--rMin -1```. Otherwise impacts will be one sided as signal strength is forced positive
* You can use the --job-mode SGE option for a shorter runtime.

```
text2workspace.py cards/2016/coupling_1/HNL_majorana_pt20_ctau1p0e01_massHNL8p0_Vall6p702e-04/out.txt -o workspace.root
combineTool.py -M Impacts -d workspace.root -m 1 --doInitialFit --robustFit 1 --expectSignal=1 -t -1
combineTool.py -M Impacts -d workspace.root -m 1 --robustFit 1 --doFits --parallel 16 --expectSignal=1 -t -1
combineTool.py -M Impacts -d workspace.root -m 1 -o impacts.json
plotImpacts.py -i impacts.json -o impacts -t rename.json --units pb
```

## Step 6: Expected & observed yields

* Run FitDiagnostics saving the shapes ```combine -M FitDiagnostics --saveWithUncertainties --saveShapes cards/combined/coupling_1/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt -t -1 --expectSignal 0.07```. For unblinded results remove the Asimov option ```-t -1```.
* Plot the results using ```python postFitPlot.py```


## Step 7: Saturated goodness-of-fit (after unblinding)

```
combine -M GoodnessOfFit cards/2016/coupling_1/HNL_majorana_pt20_ctau1p0e01_massHNL8p0_Vall6p702e-04/out.txt --algo saturated 
combine -M GoodnessOfFit cards/2016/coupling_1/HNL_majorana_pt20_ctau1p0e01_massHNL8p0_Vall6p702e-04/out.txt --algo saturated --toysFreq -t 300 --fixedSignalStrength=0
combineTool.py -M CollectGoodnessOfFit --input higgsCombineTest.GoodnessOfFit.mH120.root higgsCombineTest.GoodnessOfFit.mH120.123456.root -m 120.0 -o gof.json --fixedSignalStrength=0
plotGof.py gof.json --statistic saturated --mass 120.0 -o gof_plot --title-right="Run 2 (signal region)"
```
