import ROOT
ROOT.gROOT.SetBatch(True)
import style
ROOT.gStyle.SetErrorX(0)

import json

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
    s = s.replace("mue_OS", "#mue")
    s = s.replace("emu_OS", "e\mu")
    s = s.replace("ee_OS", "ee")
    s = s.replace("mumu_SS", "\mu\mu")
    s = s.replace("mue_SS", "#mue")
    s = s.replace("emu_SS", "e\mu")
    s = s.replace("ee_SS", "ee")
    return s

def get_hist(file_name, hist_name):
    #print(f"Reading {hist_name} from {file_name}")
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)
    hist = hist.Clone(hist_name)
    if type(hist) == ROOT.TH1F or type(hist) == ROOT.TH1D:
        hist.SetDirectory(0)
    rootFile.Close()
    return hist


def plot_yields(obs_hist, pred_hist_prefit, pred_hist, signal_hist_1, signal_hist_2, year="2016", topology="boosted"):

    colour_signal_1 = "#9e2a2b"
    colour_signal_2 = "#ff7d00"
    colour_prefit = "#15616d"

    lumi = {"2016": 36, "2017": 42, "2018": 60, "combined": 138}

    signal_hist_1.SetLineWidth(3)
    signal_hist_1.SetLineColor(ROOT.TColor.GetColor(colour_signal_1))
    signal_hist_1.SetMarkerColor(ROOT.TColor.GetColor(colour_signal_1))

    signal_hist_2.SetLineWidth(3)
    signal_hist_2.SetLineColor(ROOT.TColor.GetColor(colour_signal_2))
    signal_hist_2.SetMarkerColor(ROOT.TColor.GetColor(colour_signal_2))

    pred_hist.GetYaxis().SetTitle("Events / category")
    pred_hist.GetYaxis().SetTitleOffset(1)

    pred_hist.SetLineColor(ROOT.kAzure+1)
    pred_hist.SetMarkerSize(0)
    pred_hist.SetFillColor(ROOT.kAzure+1)

    pred_hist_prefit.SetLineColor(ROOT.TColor.GetColor(colour_prefit))
    pred_hist_prefit.SetMarkerSize(0)
    pred_hist_prefit.SetMarkerColor(ROOT.TColor.GetColor(colour_prefit))
    pred_hist_prefit.SetLineWidth(3)
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
    err_hist.SetFillColor(ROOT.kAzure+2)
    err_hist.SetFillStyle(3345)
    err_hist.SetLineColor(ROOT.kAzure+2)
    err_hist_ratio = err_hist.Clone("err_hist_ratio")

    obs_hist.SetBinErrorOption(ROOT.TH1D.kPoisson)
    residuals.SetBinErrorOption(ROOT.TH1D.kNormal)

    for j in range(pred_hist.GetNbinsX()):
        obs = obs_hist.GetBinContent(j+1)
        err_obs = obs_hist.GetBinError(j+1)
        pred = pred_hist.GetBinContent(j+1)

        ratio = obs/pred
        err_pred = pred_hist.GetBinError(j+1)

        err_hist_ratio.SetBinContent(j+1, 1)
        err_hist_ratio.SetBinError(j+1, err_pred/pred)

        residuals.SetBinError(j+1, err_obs*ratio/(obs+1e-9))

    cv = style.makeCanvas()
    cv.SetBottomMargin(0.1)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.1)

    pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.13)
    pad1.SetRightMargin(0.1)
    pad1.SetBorderMode(0)

    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.23)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.13)
    pad2.SetRightMargin(0.1)
    pad2.SetBorderMode(0)

    pad1.Draw()
    pad2.Draw()

    text = style.makeText(0.2, 0.88, 0.7, 0.88, topology)
    text.SetTextFont(63)
    text.SetTextSize(31)

    legend = style.makeLegend(0.54, 0.74, 0.79, 0.86)
    legend.AddEntry(obs_hist, 'Data', 'lp')
    legend.AddEntry(pred_hist, 'Background', "lf")
    legend.AddEntry(pred_hist_prefit, 'Background (prefit)', "l")

    legend_signal = style.makeLegend(0.14, 0.76, 0.5, 0.87)
    legend_signal.SetTextSize(19)
    legend_signal.AddEntry(signal_hist_1, "m#lower[0.3]{#scale[0.7]{N}} = 2 GeV, c#tau#lower[0.3]{#scale[0.7]{0}} = 100 mm")
    #legend_signal.AddEntry(signal_hist_2, "m#lower[0.3]{#scale[0.7]{N}} = 8 GeV, c#tau#lower[0.3]{#scale[0.7]{0}} = 0.1 mm")
    legend_signal.AddEntry(signal_hist_2, "m#lower[0.3]{#scale[0.7]{N}} = 10 GeV, c#tau#lower[0.3]{#scale[0.7]{0}} = 1 mm")


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

    text_OS = ROOT.TLatex(obs_hist.GetNbinsX()*1/4-1, 3*pred_hist.GetMaximum(), "OS")
    text_SS = ROOT.TLatex(obs_hist.GetNbinsX()*3/4-1, 3*pred_hist.GetMaximum(), "SS")

    text_prompt = ROOT.TLatex(obs_hist.GetNbinsX()*1/6-2.5, pred_hist.GetMaximum(), "d^{sig}_{xy}<1")
    text_prompt_2 = ROOT.TLatex(obs_hist.GetNbinsX()*4/6-2.5, pred_hist.GetMaximum(), "d^{sig}_{xy}<1")

    text_medium = ROOT.TLatex(obs_hist.GetNbinsX()*2/6-2.5, pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")
    text_medium_2 = ROOT.TLatex(obs_hist.GetNbinsX()*5/6-2.5, pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")

    text_displaced = ROOT.TLatex(obs_hist.GetNbinsX()*3/6-2.5, pred_hist.GetMaximum(), "d^{sig}_{xy}>10")
    text_displaced_2 = ROOT.TLatex(obs_hist.GetNbinsX()*6/6-2.5, pred_hist.GetMaximum(), "d^{sig}_{xy}>10")

    text_prompt.SetTextAlign(22)
    text_prompt_2.SetTextAlign(22)
    text_medium.SetTextAlign(22)
    text_medium_2.SetTextAlign(22)
    text_displaced.SetTextAlign(22)
    text_displaced_2.SetTextAlign(22)

    text_OS.SetTextAlign(22)
    text_SS.SetTextAlign(22)

    #legend.AddEntry(err_hist, 'stat. unc', 'f')

    pad2.cd()
    residuals.GetYaxis().SetTitle("Obs/Pred")
    residuals.GetYaxis().SetTitleOffset(1)
    residuals.GetYaxis().SetNdivisions(504)
    #n1 + 100*n2
    residuals.SetMarkerColor(1)
    residuals.SetLineColor(1)
    residuals.GetYaxis().SetRangeUser(0, 5.3)
    residuals.Draw("P")
    lres = ROOT.TLine(-0.5, 1., n_bins_total-0.5, 1.)
    lres.Draw("SAME")
    err_hist_ratio.Draw("E2 SAME")

    pad1.cd()
    pad1.SetLogy()
    pred_hist.SetMinimum(1e-3)
    pred_hist_prefit.SetMinimum(1e-3)
    err_hist.SetMinimum(1e-3)
    obs_hist.SetMinimum(1e-3)
    signal_hist_1.SetMinimum(1e-3)
    signal_hist_2.SetMinimum(1e-3)
    pred_hist.GetYaxis().SetRangeUser(0.3, 300*pred_hist.GetMaximum())
    pred_hist.GetXaxis().SetLabelSize(0)

    pred_hist.Draw("HIST")
    err_hist.Draw("E2 SAME")
    pred_hist_prefit.Draw("HIST SAME")
    obs_hist.Draw("HIST P L E SAME")
    signal_hist_1.Draw("HIST SAME")
    signal_hist_2.Draw("HIST SAME")

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
    #style.makeCMSText(0.15, 0.91, additionalText="Preliminary")
    style.makeLumiText(0.9, 0.96, year=year, lumi=lumi[year])

    cv.SaveAs("yields/{}_{}.pdf".format(topology, year))
    cv.SaveAs("yields/{}_{}.png".format(topology, year))

categories = [
    "mumu_OS", "ee_OS", "mue_OS", "emu_OS", 
    "mumu_SS", "ee_SS", "mue_SS", "emu_SS"
]

n_bins_displacement = 3
n_bins_total = len(categories)*n_bins_displacement
category_names = []
category_names_raw = []

root_file = "fitDiagnostics.root"
signal = "HNL_majorana_all_ctau1p0e02_massHNL2p0_Vall9p078e-03"
#signal_2 = "HNL_majorana_all_ctau1p0e-01_massHNL8p0_Vall6p702e-03"
signal_2 = "HNL_majorana_all_ctau1p0e00_massHNL10p0_Vall1p177e-03"

signal_hist_path = "/vols/cms/vc1117/AN-19-207/histo/limits/hists_merged/"

root_file_signal_1 = signal_hist_path+signal+"_YEAR.root"
root_file_signal_2 = signal_hist_path+signal_2+"_YEAR.root"

coupling = "1"

k_factor = 1.1

with open("/vols/cms/LLP/gridpackLookupTable.json") as lookup_table_file:
    lookup_table = json.load(lookup_table_file)
    lu_infos = lookup_table[signal]['weights'][str(int(coupling))]
    lu_infos_2 = lookup_table[signal_2]['weights'][str(int(coupling))]

    signal_xsec = lu_infos['xsec']['nominal']
    signal_xsec_2 = lu_infos_2['xsec']['nominal']

for topology in ["boosted", "resolved"]:
    output_hists_by_year = {}

    for i, year in enumerate(["2016", "2017", "2018", "combined"]):
        output_hists_by_year[year] = {}
        output_hists_by_year[year]["pred"] = []
        output_hists_by_year[year]["pred_prefit"] = []
        output_hists_by_year[year]["obs"] = []
        output_hists_by_year[year]["signal"] = []
        output_hists_by_year[year]["signal2"] = []

        if year != "combined":
            hist_pred = ROOT.TH1D("pred"+topology, "pred"+topology, n_bins_total, -0.5, n_bins_total-0.5)
            hist_obs = hist_pred.Clone("obs"+topology)
            hist_pred_prefit = hist_pred.Clone("pred_prefit"+topology)

            hist_signal = hist_pred.Clone("signal"+topology)
            hist_signal2 = hist_pred.Clone("signal2"+topology)

            print(year, topology)
            if topology == "boosted":
                indices = [1, 2, 3]
            else:
                indices = [4, 5, 6]

            for j, category in enumerate(categories):
                category_names_raw.append(category)
                category_names.append(make_pretty(category))

                hist_pred_input = get_hist(root_file, "shapes_fit_b/ch"+str(i+1)+"_"+category+"_D/total_background")
                hist_pred_prefit_input = get_hist(root_file, "shapes_prefit/ch"+str(i+1)+"_"+category+"_D/total_background")
                hist_obs_input = get_hist(root_file, "shapes_fit_b/ch"+str(i+1)+"_"+category+"_D/data")
                hist_signal_input = get_hist(root_file_signal_1.replace("YEAR", year), category+"_D/HNL_coupling_1")
                hist_signal2_input = get_hist(root_file_signal_2.replace("YEAR", year), category+"_D/HNL_coupling_1")

                for idx in indices:
                    pred = hist_pred_input.GetBinContent(idx)
                    pred_up = hist_pred_input.GetBinErrorUp(idx)
                    pred_low = hist_pred_input.GetBinErrorLow(idx)

                    pred_prefit = hist_pred_prefit_input.GetBinContent(idx)
                    pred_up_prefit = hist_pred_prefit_input.GetBinErrorUp(idx)
                    pred_low_prefit = hist_pred_prefit_input.GetBinErrorLow(idx)

                    obs = hist_obs_input.GetY()
                    obs = obs[idx-1]
                    obs_up = hist_obs_input.GetErrorYhigh(idx)
                    obs_low = hist_obs_input.GetErrorYlow(idx)

                    signal = hist_signal_input.GetBinContent(idx)
                    signal2 = hist_signal2_input.GetBinContent(idx)

                    if "SS" in category:
                        j_hack = j + 8
                    else:
                        j_hack = j

                    if topology == "boosted":
                        hist_index = (idx-1)*4 + j_hack + 1
                    else:
                        hist_index = (idx-4)*4 + j_hack + 1

                    hist_pred.SetBinContent(hist_index, pred)
                    hist_pred.SetBinError(hist_index, pred_up)
                    hist_obs.SetBinContent(hist_index, obs)

                    hist_pred_prefit.SetBinContent(hist_index, pred_prefit)
                    hist_obs.GetXaxis().SetBinLabel(hist_index, shorten(category))
                    hist_pred.GetXaxis().SetBinLabel(hist_index, shorten(category))

                    hist_signal.SetBinContent(hist_index, signal)
                    hist_signal2.SetBinContent(hist_index, signal2)

            hist_signal.Scale(k_factor*signal_xsec)
            hist_signal2.Scale(k_factor*signal_xsec_2)

            output_hists_by_year[year]["pred"] = hist_pred
            output_hists_by_year[year]["pred_prefit"] = hist_pred_prefit
            output_hists_by_year[year]["obs"] = hist_obs
            output_hists_by_year[year]["signal"] = hist_signal
            output_hists_by_year[year]["signal2"] = hist_signal2
        else:
            hist_pred = output_hists_by_year["2016"]["pred"]
            hist_pred_2017 = output_hists_by_year["2017"]["pred"]
            hist_pred_2018 = output_hists_by_year["2018"]["pred"]
            hist_pred.Add(hist_pred_2017)
            hist_pred.Add(hist_pred_2018)

            hist_pred_prefit = output_hists_by_year["2016"]["pred_prefit"]
            hist_pred_prefit_2017 = output_hists_by_year["2017"]["pred_prefit"]
            hist_pred_prefit_2018 = output_hists_by_year["2018"]["pred_prefit"]
            hist_pred_prefit.Add(hist_pred_prefit_2017)
            hist_pred_prefit.Add(hist_pred_prefit_2018)

            hist_obs = output_hists_by_year["2016"]["obs"]
            hist_obs_2017 = output_hists_by_year["2017"]["obs"]
            hist_obs_2018 = output_hists_by_year["2018"]["obs"]
            hist_obs.Add(hist_obs_2017)
            hist_obs.Add(hist_obs_2018)

            hist_signal = output_hists_by_year["2016"]["signal"]
            hist_signal_2017 = output_hists_by_year["2017"]["signal"]
            hist_signal_2018 = output_hists_by_year["2018"]["signal"]
            hist_signal.Add(hist_signal_2017)
            hist_signal.Add(hist_signal_2018)

            hist_signal2 = output_hists_by_year["2016"]["signal2"]
            hist_signal2_2017 = output_hists_by_year["2017"]["signal2"]
            hist_signal2_2018 = output_hists_by_year["2018"]["signal2"]
            hist_signal2.Add(hist_signal2_2017)
            hist_signal2.Add(hist_signal2_2018)

            # add hists for all years

        plot_yields(hist_obs, hist_pred_prefit, hist_pred, hist_signal, hist_signal2, year=year, topology=topology)


#     print()
#     s = """\\begin{tabular}{llllll}
# $\ell_{1}\ell_{2}$ & topology &  & \multicolumn{3}{c}{$d_{xy}^\\mathrm{sig}$ } \\\\
# \\cline{4-6}
# &          &  & $d_{xy}^\\mathrm{sig}<1$       & $1<d_{xy}^\\mathrm{sig}<10$       & $d_{xy}^\\mathrm{sig}>10$     \\\\
# \\hline"""
#     print(s)
#             pred_string =  '${:.2f}^{{+{:.2f}}}_{{-{:.2f}}}$'.format(pred, pred_up, pred_low)

# print("{} & boosted & pred. & {} & {} & {} \\\\").format(make_pretty(category), pred_strings_merged[0], pred_strings_merged[1], pred_strings_merged[2])
# print("& & signal & {} & {} & {} \\\\").format(signal_strings_merged[0], signal_strings_merged[1], signal_strings_merged[2])
# print("& resolved & pred. & {} & {} & {} \\\\").format(pred_strings_resolved[0], pred_strings_resolved[1], pred_strings_resolved[2])
# print("& & signal & {} & {} & {} \\\\").format(signal_strings_resolved[0], signal_strings_resolved[1], signal_strings_resolved[2])
# print("\end{tabular}")

