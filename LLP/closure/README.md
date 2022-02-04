# Closure checks (optional)

## Step 1: Make cards in parent directory

```
cd ...
python makeCardsClosure.py
cd closure
```

## Step 2: Combine datacards and produce workspace

```
rm combined.txt
combineCards.py 201*.txt>> combined.txt 
```

## Step 3: perform fit, save pre- and post-fit shapes for plotting

```combine -M FitDiagnostics -m 125 --saveShapes --saveWithUncertainties combined.txt```

## Step 4: plot closure

```
python plotClosure.py
```

## Step 5: Saturated goodness-of-fit

```
text2workspace.py combined.txt -o workspace.root
combine -M GoodnessOfFit combined.txt --algo saturated --fixedSignalStrength=0
combine -M GoodnessOfFit combined.txt --algo saturated --toysFreq -t 300 --fixedSignalStrength=0
combineTool.py -M CollectGoodnessOfFit --input higgsCombineTest.GoodnessOfFit.mH120.root higgsCombineTest.GoodnessOfFit.mH120.123456.root -m 120.0 -o gof.json
plotGof.py gof.json --statistic saturated --mass 120.0 -o gof_plot --title-right="Run 2 (VR1)"
```
