# Setting limits for HNL analysis


## Step 1: Produce histograms

* Have a look at ```make_hists.py``` script in [here](https://github.com/LLPDNNX/histo/blob/master/make_hists.py)

## Step 2: Produce datacards and run combine

* Running ```python makeCards.py``` produces a datacard per HNL coupling, mass and scenario point. The cards will be located in the ```cards``` directory.
* a SGE batch submission file is created as well: ```runCombine.sh```. You can test a single datacard with ```SGE_TASK_ID=1 ./runCombine.sh``` before submitting.

## Step 3: Do some sanity checks:

* ```combine -M FitDiagnostics -t -1 --expectSignal 0 cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt```. Expect to see best fit r=0.
* ```combine -M FitDiagnostics -t -1 --expectSignal 1 cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt```. Expect to see best fit r=1.


## Step 4: Collect results

* parse results into json files, e.g. ```for YEAR in 2016 2017 2018 combined; do combineTool.py -M CollectLimits cards/$YEAR/coupling_*/*/*HNL*.root --use-dirs -o jsons/limits_$YEAR.json; done```.
* plot limits using ```python CombineHarvester/LLP/plotLimits.py```

## Step 5: Make impact plots for one benchmark model (takes a while)
* Optionally, make a rename dictionary using makePUrenameDict.py for aesthetic purposes
* Run the following set of commands. The runtime is very long (at least 30 mins) for all three year combination, consider doing only one.
```
text2workspace.py cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt -o workspace.root
combineTool.py -M Impacts -d workspace.root -m 7 -t -1 --doInitialFit --robustFit 1 --expectSignal=1
combineTool.py -M Impacts -d workspace.root -m 7 -t -1 --robustFit 1 --doFits --parallel 16 --expectSignal=1
combineTool.py -M Impacts -d workspace.root -m 7 -o impacts.json
plotImpacts.py -i impacts.json -o impacts -t rename.json --units pb
```

## Step 6: Get expected and observed yields

* Run FitDiagnostics saving the shapes ```combine -M FitDiagnostics --saveWithUncertainties --saveShapes cards/combined/coupling_1/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt -t -1 --expectSignal 0.07```. For unblinded results remove the Asimov option ```-t -1```.
* Plot the results using ```python postFitPlot.py```
