import pandas as pd
import numpy as np
import json
import ROOT
import math
from scipy.stats import chi2


category_names = []
category_names_raw = []
nbins = 6


def make_pretty(s):

    """ Make pretty the category text"""
    s = s.replace("mumu_OS", "$\mu^{\pm}\mu^{\mp}$")
    s = s.replace("mue_OS", "$\mu^{\pm}e^{\mp}$")
    s = s.replace("emu_OS", "$e^{\pm}\mu^{\mp}$")
    s = s.replace("ee_OS", "$e^{\pm}e^{\mp}$")
    s = s.replace("mumu_SS", "$\mu^{\pm}\mu^{\pm}$")
    s = s.replace("mue_SS", "$\mu^{\pm}e^{\pm}$")
    s = s.replace("emu_SS", "$e^{\pm}\mu^{\pm}$")
    s = s.replace("ee_SS", "$e^{\pm}e^{\pm}$")
    return s

def get_hist(file_name, hist_name):
    #print(f"Reading {hist_name} from {file_name}")
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)
    hist = hist.Clone()
    hist.SetDirectory(0)
    rootFile.Close()
    return hist


categories = [
    "mumu_OS", "ee_OS", "mue_OS", "emu_OS", 
    "mumu_SS", "ee_SS", "mue_SS", "emu_SS"
]

for year in ["2016"]:
    print(year)
    print()
    s = """\\begin{tabular}{llllll}
$\ell_{1}\ell_{2}$ & topology &  & \multicolumn{3}{c}{$d_{xy}^\\mathrm{sig}$ } \\\\
\\cline{4-6}
&          &  & $d_{xy}^\\mathrm{sig}<1$       & $1<d_{xy}^\\mathrm{sig}<10$       & $d_{xy}^\\mathrm{sig}>10$     \\\\
\\hline"""
    print(s)
    for category in categories:
        category_names_raw.append(category)
        category_names.append(make_pretty(category))
        hist_pred = get_hist("shapes.root", category+"_D_postfit/TotalBkg")
        hist_obs = get_hist("shapes.root", category+"_D_postfit/data_obs")

        obs_strings_merged = []
        pred_strings_merged = []
        obs_strings_resolved = []
        pred_strings_resolved = []

        for idx in [1, 2, 3, 4, 5, 6]:
            pred = hist_pred.GetBinContent(idx)
            pred_up = hist_pred.GetBinErrorUp(idx)
            pred_low = hist_pred.GetBinErrorLow(idx)
            pred_string =  '${:.1f}^{{+{:.1f}}}_{{-{:.1f}}}$'.format(pred, pred_up, pred_low)

            obs = hist_obs.GetBinContent(idx)
            obs_up = hist_obs.GetBinErrorUp(idx)
            obs_low = hist_obs.GetBinErrorLow(idx)
            obs_string =  '${:.1f}^{{+{:.1f}}}_{{-{:.1f}}}$'.format(obs, obs_up, obs_low)
            if idx in [1, 2, 3]:
                pred_strings_merged.append(pred_string)
                obs_strings_merged.append(obs_string)
            elif idx in [4, 5, 6]:
                pred_strings_resolved.append(pred_string)
                obs_strings_resolved.append(obs_string)

        print("{} & merged & pred. & {} & {} & {} \\\\").format(make_pretty(category), pred_strings_merged[0], pred_strings_merged[1], pred_strings_merged[2])
        print("& & obs. & {} & {} & {} \\\\").format(pred_strings_merged[0], pred_strings_merged[1], pred_strings_merged[2])
        print("& resolved & pred. & {} & {} & {} \\\\").format(pred_strings_resolved[0], pred_strings_resolved[1], pred_strings_resolved[2])
        print("& & obs. & {} & {} & {} \\\\").format(pred_strings_resolved[0], pred_strings_resolved[1], pred_strings_resolved[2])
    print("\end{tabular}")
            


# mu1mu2 & 0.4           & obs  &         &         &        \\
#        &               & pred &         &         &       
        #hist_signal = get_hist(f"hists_merged/HNL_majorana_all_ctau1p0e00_massHNL8p0_Vall2p119e-03_{year}.root", f"{category}_D/HNL_coupling_7")

            #hist_signal.Scale(0.1932398841208851)

            #signal_yield = hist_signal.GetBinContent(idx)
            #signal_uncertainty = hist_signal.GetBinError(idx)
            #signal_string = f"${signal_yield:0.1f} \\pm {signal_uncertainty:0.1f}$"

            #signal_yields[year].append(signal_yield)
            #signal_yield_errors[year].append(signal_uncertainty)
            #signal_yields_string[year].append(signal_string)

# df_plot = pd.DataFrame(dict(name=category_names_raw, 
# pred2016=yields["2016"], pred2017=yields["2017"], pred2018=yields["2018"],
# errPred2016=yield_errors["2016"], errPred2017=yield_errors["2017"], errPred2018=yield_errors["2018"],
# signal2016=signal_yields["2016"], signal2017=signal_yields["2017"], signal2018=signal_yields["2018"],
# errSignal2016=signal_yield_errors["2016"], errSignal2017=signal_yield_errors["2017"], errSignal2018=signal_yield_errors["2018"],
# )
# )
# print(df_plot)
# #plot_yields(df_plot, topology=topology)
# df = pd.DataFrame(dict(name=category_names, pred2016=yields["2016"], pred2017=yields["2017"], pred2018=yields["2018"], signal2016=signal_yields["2016"], signal2017=signal_yields["2017"], signal2018=signal_yields["2018"]))
# print(df.to_latex(index=False, escape=False, caption=f"Expected background and signal yields for {topology} categories. The signal model shown is Majorana HNL of $8\\GeV, c\\tau_0=1\\unit{{mm}}, V_\\mu=V_e$.", label=f"tab:yields_{topology}", column_format='c|c|c|c|c|c|c'))
