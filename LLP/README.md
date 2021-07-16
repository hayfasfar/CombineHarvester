# Setting limits for HNL analysis

## Step 1: Produce histograms

* Have a look at ```make_hists.py``` script in [here](https://github.com/LLPDNNX/histo/blob/master/make_hists.py)

## Step 2: Produce datacards and run combine

* Running ```python makeCards.py``` produces a datacard per HNL coupling, mass and scenario point. The cards will be located in the ```cards``` directory.
* a SGE batch submission file is created as well: ```runCombine.sh```. You can test a single datacard with ```SGE_TASK_ID=1 ./runCombine.sh``` before submitting.

## Step 3: Collect results

* parse results into json files, e.g. ```for YEAR in 2016 2017 2018 combined; do combineTool.py -M CollectLimits cards/$YEAR/coupling_*/*/*HNL*.root --use-dirs -o jsons/limits_$YEAR.json; done```.
* plot limits using ```python CombineHarvester/LLP/plotLimits.py```

## Step 4: Make impact plots for one benchmark model (takes a while)
* Optionally, make a rename dictionary using makePUrenameDict.py
* Run the following set of commands. The runtime is a few minutes

```
text2workspace.py cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt -o workspace.root
combineTool.py -M Impacts -d workspace.root -m 7 -t -1 --doInitialFit --robustFit 1 --expectSignal=1
combineTool.py -M Impacts -d workspace.root -m 7 -t -1 --robustFit 1 --doFits --parallel 16 --expectSignal=1
combineTool.py -M Impacts -d workspace.root -m 7 -o impacts.json
plotImpacts.py -i impacts.json -o impacts -t rename.json --units pb
```

## Step 5: Get expected and observed (once unblinded) yields

* Make workspace: ```text2workspace.py cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt -o workspace.root --channel-masks```
* Do a bkg-only fit ```combine -M MultiDimFit -t -1 --expectSignal=0 -d cards/2016/coupling_7/HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03/out.txt --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=99999999999 --saveFitResult```
* Make pre and post-fit shape histograms ```PostFitShapesFromWorkspace -w workspace.root -o shapes.root -m 12 -f multidimfit.root:fit_mdf --postfit --sampling 1 --print```
* Plot them using ```python postFitPlot.py```
