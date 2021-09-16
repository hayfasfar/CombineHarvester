import pandas as pd
import numpy as np
import json
import ROOT
import math
from scipy.stats import chi2
import style

def make_pretty(s):
    """ Make pretty the category text"""
    s = s.replace("mumu_OS", "$\mu^{\pm}\mu$")
    s = s.replace("mue_OS", "$\mu^{\pm}e$")
    s = s.replace("emu_OS", "$e^{\pm}\mu^{\mp}$")
    s = s.replace("ee_OS", "$e^{\pm}e^{\mp}$")
    s = s.replace("mumu_SS", "$\mu^{\pm}\mu^{\pm}$")
    s = s.replace("mue_SS", "$\mu^{\pm}e^{\pm}$")
    s = s.replace("emu_SS", "$e^{\pm}\mu^{\pm}$")
    s = s.replace("ee_SS", "$e^{\pm}e^{\pm}$")
    return s

def shorten(s):
    """ Make pretty the category text"""
    s = s.replace("mumu_OS", "\mu\mu")
    s = s.replace("mue_OS", "\mu e")
    s = s.replace("emu_OS", "e\mu")
    s = s.replace("ee_OS", "ee")
    s = s.replace("mumu_SS", "\mu\mu")
    s = s.replace("mue_SS", "\mu e")
    s = s.replace("emu_SS", "e\mu")
    s = s.replace("ee_SS", "ee")
    return s

def get_hist(file_name, hist_name):
    #print(f"Reading {hist_name} from {file_name}")
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)
    hist = hist.Clone(hist_name)
    if type(hist) == ROOT.TH1F:
        hist.SetDirectory(0)
    rootFile.Close()
    return hist

    

def plot_yields(obs_hist, pred_hist_prefit, pred_hist, signal_hist, year="2016", topology="boosted"):
    lumi = {"2016": 35.9, "2017": 41.5, "2018": 59.7, "combined": 137.1}

    signal_hist.SetLineWidth(2)
    signal_hist.SetLineColor(ROOT.kRed)
    signal_hist.SetMarkerColor(ROOT.kRed)


    pred_hist.GetYaxis().SetTitle("Events / category")
    pred_hist.GetYaxis().SetTitleOffset(1)

    pred_hist.SetLineColor(ROOT.kAzure+1)
    pred_hist.SetMarkerColor(ROOT.kAzure+1)
    pred_hist.SetFillColor(ROOT.kAzure+1)

    pred_hist_prefit.SetLineColor(ROOT.kAzure-2)
    pred_hist_prefit.SetMarkerSize(0)
    pred_hist_prefit.SetMarkerColor(ROOT.kAzure-2)
    pred_hist_prefit.SetLineWidth(2)
    pred_hist_prefit.SetLineStyle(2)

    pred_hist.SetLineWidth(2)

    obs_hist.SetLineWidth(2)
    pred_hist.SetMaximum(2.5*max(pred_hist.GetMaximum(), obs_hist.GetMaximum()))

    pred_hist.SetMinimum(0)
    obs_hist.SetMinimum(0)
    residuals = obs_hist.Clone("residuals")
    residuals.Divide(pred_hist)
    residuals.GetXaxis().LabelsOption("v")

    err_hist = pred_hist.Clone("")

    cv = style.makeCanvas()
    cv.SetBottomMargin(0.1)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.1)

    pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.13)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)

    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.25)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.13)
    pad2.SetRightMargin(0.1)
    pad2.SetBorderMode(0)

    pad1.Draw()
    pad2.Draw()

    style.makeText(0.2, 0.95, 0.7, 0.95, topology)
    legend = style.makeLegend(0.5, 0.74, 0.8, 0.86)
    legend.AddEntry(obs_hist, 'Data', 'lp')
    legend.AddEntry(pred_hist, 'Background', "lfp")
    legend.AddEntry(pred_hist_prefit, 'Background (prefit)', "l")

    legend_signal = style.makeLegend(0.15, 0.67, 0.7, 0.74)
    legend_signal.AddEntry(signal_hist, "m_{N} = 10 GeV, c#tau_{0} = 1 mm, V_{e} = V_{#mu} = V_{#tau}")


    l_OS_SS = ROOT.TLine(obs_hist.GetNbinsX()/2-0.5, obs_hist.GetMinimum(), obs_hist.GetNbinsX()/2-0.5, 0.6*pred_hist.GetMaximum())
    l_prompt_displaced = ROOT.TLine(obs_hist.GetNbinsX()/6-0.5, obs_hist.GetMinimum(), obs_hist.GetNbinsX()/6-0.5, 0.5*pred_hist.GetMaximum())
    l_prompt_displaced_2 = ROOT.TLine(obs_hist.GetNbinsX()*2/6-0.5, obs_hist.GetMinimum(), obs_hist.GetNbinsX()*2/6-0.5, 0.5*pred_hist.GetMaximum())
    l_prompt_displaced_3 = ROOT.TLine(obs_hist.GetNbinsX()*4/6-0.5, obs_hist.GetMinimum(), obs_hist.GetNbinsX()*4/6-0.5, 0.5*pred_hist.GetMaximum())
    l_prompt_displaced_4 = ROOT.TLine(obs_hist.GetNbinsX()*5/6-0.5, obs_hist.GetMinimum(), obs_hist.GetNbinsX()*5/6-0.5, 0.5*pred_hist.GetMaximum())
    l_prompt_displaced.SetLineWidth(2)
    l_prompt_displaced_2.SetLineWidth(2)
    l_prompt_displaced_3.SetLineWidth(2)
    l_prompt_displaced_4.SetLineWidth(2)
    l_prompt_displaced.SetLineStyle(3)
    l_prompt_displaced_2.SetLineStyle(3)
    l_prompt_displaced_3.SetLineStyle(3)
    l_prompt_displaced_4.SetLineStyle(3)
    l_OS_SS.SetLineStyle(2)
    l_OS_SS.SetLineWidth(3)

    text_OS = ROOT.TLatex(obs_hist.GetNbinsX()*1/4-1, 0.6*pred_hist.GetMaximum(), "OS")
    text_SS = ROOT.TLatex(obs_hist.GetNbinsX()*3/4-1, 0.6*pred_hist.GetMaximum(), "SS")

    text_prompt = ROOT.TLatex(obs_hist.GetNbinsX()*1/6-2.5, 0.45*pred_hist.GetMaximum(), "d^{sig}_{xy}<1")
    text_prompt_2 = ROOT.TLatex(obs_hist.GetNbinsX()*4/6-2.5, 0.45*pred_hist.GetMaximum(), "d^{sig}_{xy}<1")

    text_medium = ROOT.TLatex(obs_hist.GetNbinsX()*2/6-2.5, 0.45*pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")
    text_medium_2 = ROOT.TLatex(obs_hist.GetNbinsX()*5/6-2.5, 0.45*pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")

    text_displaced = ROOT.TLatex(obs_hist.GetNbinsX()*3/6-2.5, 0.45*pred_hist.GetMaximum(), "d^{sig}_{xy}>10")
    text_displaced_2 = ROOT.TLatex(obs_hist.GetNbinsX()*6/6-2.5, 0.45*pred_hist.GetMaximum(), "d^{sig}_{xy}>10")

    text_prompt.SetTextAlign(22)
    text_prompt_2.SetTextAlign(22)
    text_medium.SetTextAlign(22)
    text_medium_2.SetTextAlign(22)
    text_displaced.SetTextAlign(22)
    text_displaced_2.SetTextAlign(22)

    text_OS.SetTextAlign(22)
    text_SS.SetTextAlign(22)

    err_hist.SetFillColor(ROOT.kAzure+2)
    err_hist.SetFillStyle(3345)
    err_hist.SetLineColor(ROOT.kAzure+2)
    #legend.AddEntry(err_hist, 'stat. unc', 'f')

    pad2.cd()
    residuals.GetYaxis().SetTitle("Obs/Pred")
    residuals.GetYaxis().SetTitleOffset(1)
    residuals.GetYaxis().SetNdivisions(504)
    #n1 + 100*n2
    residuals.SetMarkerColor(1)
    residuals.SetLineColor(1)
    residuals.GetYaxis().SetRangeUser(-0.3, 3.5)
    residuals.Draw("P")
    lres = ROOT.TLine(-0.5, 1., n_bins_total-0.5, 1.)
    lres.Draw("SAME")

    pad1.cd()
    pred_hist.Draw("HIST")
    err_hist.Draw("E2 SAME")
    pred_hist_prefit.Draw("HIST SAME")
    obs_hist.Draw("HIST P L E SAME")
    signal_hist.Draw("HIST SAME")
    l_prompt_displaced.Draw("SAME")
    l_prompt_displaced_2.Draw("SAME")
    l_prompt_displaced_3.Draw("SAME")
    l_prompt_displaced_4.Draw("SAME")
    l_OS_SS.Draw("SAME")
    text_prompt.Draw("SAME")
    text_prompt_2.Draw("SAME")
    text_medium.Draw("SAME")
    text_medium_2.Draw("SAME")
    text_displaced.Draw("SAME")
    text_displaced_2.Draw("SAME")
    text_OS.Draw("SAME")
    text_SS.Draw("SAME")

    cv.cd()
    legend.Draw("")
    legend_signal.Draw("")
    style.makeCMSText(0.15, 0.91, additionalText="Preliminary")
    style.makeLumiText(0.6, 0.91, year=year, lumi=lumi[year])

    cv.SaveAs("yields/{}_{}.pdf".format(topology, year))
    cv.SaveAs("yields/{}_{}.png".format(topology, year))



categories = [
    "mumu_OS", "ee_OS", "mue_OS", "emu_OS", 
    "mumu_SS", "ee_SS", "mue_SS", "emu_SS"
]

n_bins_total = len(categories)*3

category_names = []
category_names_raw = []
nbins = 6

toys = True


root_file = "fitDiagnostics.root"


for i, year in enumerate(["2016", "2017", "2018"]):
    hist_pred_merged = ROOT.TH1D("pred_merged", "pred_merged", n_bins_total, -0.5, n_bins_total-0.5)
    hist_obs_merged = hist_pred_merged.Clone("obs_merged")
    hist_pred_resolved = hist_pred_merged.Clone("pred_resolved")
    hist_obs_resolved = hist_pred_merged.Clone("obs_resolved")
    hist_pred_prefit_merged = hist_pred_merged.Clone("pred_merged_prefit")
    hist_pred_prefit_resolved = hist_pred_merged.Clone("pred_resolved_prefit")
    hist_signal_merged = hist_pred_merged.Clone("signal_merged")
    hist_signal_resolved = hist_pred_resolved.Clone("signal_resolved")

    print(year)
    print()
    s = """\\begin{tabular}{llllll}
$\ell_{1}\ell_{2}$ & topology &  & \multicolumn{3}{c}{$d_{xy}^\\mathrm{sig}$ } \\\\
\\cline{4-6}
&          &  & $d_{xy}^\\mathrm{sig}<1$       & $1<d_{xy}^\\mathrm{sig}<10$       & $d_{xy}^\\mathrm{sig}>10$     \\\\
\\hline"""
    print(s)
    for j, category in enumerate(categories):
        category_names_raw.append(category)
        category_names.append(make_pretty(category))

        hist_pred = get_hist(root_file, "shapes_fit_b/ch"+str(i+1)+"_"+category+"_D/total_background")
        hist_pred_prefit = get_hist(root_file, "shapes_prefit/ch"+str(i+1)+"_"+category+"_D/total_background")
        hist_obs = get_hist(root_file, "shapes_fit_b/ch"+str(i+1)+"_"+category+"_D/data")
        hist_signal = get_hist(root_file, "shapes_prefit/ch"+str(i+1)+"_"+category+"_D/HNL")

        obs_strings_merged = []
        pred_strings_merged = []
        obs_strings_resolved = []
        pred_strings_resolved = []

        for idx in [1, 2, 3, 4, 5, 6]:
            pred = hist_pred.GetBinContent(idx)
            pred_up = hist_pred.GetBinErrorUp(idx)
            pred_low = hist_pred.GetBinErrorLow(idx)
            pred_string =  '${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$'.format(pred, pred_up, pred_low)

            pred_prefit = hist_pred_prefit.GetBinContent(idx)
            pred_up_prefit = hist_pred_prefit.GetBinErrorUp(idx)
            pred_low_prefit = hist_pred_prefit.GetBinErrorLow(idx)

            obs = hist_obs.GetY()
            obs = obs[idx-1]
            obs_up = hist_obs.GetErrorYhigh(idx)
            obs_low = hist_obs.GetErrorYlow(idx)
            obs_string =  '${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$'.format(obs, obs_up, obs_low)

            signal = hist_signal.GetBinContent(idx)

            if "SS" in category:
                j_hack = j + 8
            else:
                j_hack = j

            if idx in [1, 2, 3]:
                hist_index = (idx-1)*4 + j_hack + 1
                pred_strings_merged.append(pred_string)
                obs_strings_merged.append(obs_string)
                hist_pred_merged.SetBinContent(hist_index, pred)
                hist_pred_merged.SetBinError(hist_index, pred_up)
                hist_obs_merged.SetBinContent(hist_index, obs)

                hist_pred_prefit_merged.SetBinContent(hist_index, pred_prefit)
                hist_obs_merged.GetXaxis().SetBinLabel(hist_index, shorten(category))
                hist_pred_merged.GetXaxis().SetBinLabel(hist_index, shorten(category))

                hist_signal_merged.SetBinContent(hist_index, signal)


            elif idx in [4, 5, 6]:
                hist_index = (idx-4)*4 + j_hack + 1
                pred_strings_resolved.append(pred_string)
                obs_strings_resolved.append(obs_string)
                hist_pred_resolved.SetBinContent(hist_index, pred)
                hist_pred_resolved.SetBinError(hist_index, pred_up)
                hist_obs_resolved.SetBinContent(hist_index, obs)
                hist_pred_prefit_resolved.SetBinContent(hist_index, pred_prefit)
                hist_obs_resolved.GetXaxis().SetBinLabel(hist_index, shorten(category))
                hist_pred_resolved.GetXaxis().SetBinLabel(hist_index, shorten(category))
                hist_signal_resolved.SetBinContent(hist_index, signal)


        print("{} & boosted & pred. & {} & {} & {} \\\\").format(make_pretty(category), pred_strings_merged[0], pred_strings_merged[1], pred_strings_merged[2])
        print("& & obs. & {} & {} & {} \\\\").format(obs_strings_merged[0], obs_strings_merged[1], obs_strings_merged[2])
        print("& resolved & pred. & {} & {} & {} \\\\").format(pred_strings_resolved[0], pred_strings_resolved[1], pred_strings_resolved[2])
        print("& & obs. & {} & {} & {} \\\\").format(obs_strings_resolved[0], obs_strings_resolved[1], obs_strings_resolved[2])
    print("\end{tabular}")
    plot_yields(hist_obs_merged, hist_pred_prefit_merged, hist_pred_merged, hist_signal_merged, year=year, topology="boosted")
    plot_yields(hist_obs_resolved, hist_pred_prefit_resolved, hist_pred_resolved, hist_signal_resolved, year=year, topology="resolved")