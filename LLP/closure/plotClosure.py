import ROOT
ROOT.gROOT.SetBatch(True)
import math
import style
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def mean_and_sigma(pulls):
    """
    Takes a list of pulls
    Returns mean and standard deviation and their uncertainties (+- 1 sigma)
    """
    pulls = np.array(pulls)
    mean = np.mean(pulls)
    sigma = np.std(pulls, ddof=1)
    n = len(pulls)
    unc_mean = sigma/math.sqrt(n)
    unc_sigma = sigma/math.sqrt(2*n-2)
    return mean, unc_mean, sigma, unc_sigma

def plot_pulls(pulls, pulls_prefit, save_text):
    meanPull, uncMeanPull, sigmaPull, uncSigmaPull = mean_and_sigma(pulls)
    meanPullPre, uncMeanPullPre, sigmaPullPre, uncSigmaPullPre = mean_and_sigma(pulls_prefit)

    bins = np.linspace(-5.5, 5.5, num=12)
    x = np.linspace(-5, 5, num=100)

    fig, ax = plt.subplots()
    #W, p = scipy.stats.kstest(pulls, 'norm')
    plt.hist(pulls, bins,  histtype='step', color="red", lw=2)

    pdf = stats.norm.pdf(x, meanPull, sigmaPull)*len(pulls)
    plt.plot(x, pdf, color="red", lw=3, ls='--', label="postfit: $\mu = %.2f \pm %.2f, \sigma=%.2f \pm %.2f$" % (meanPull, uncMeanPull, sigmaPull, uncSigmaPull))

    #pdfPre = stats.norm.pdf(x, meanPullPre, sigmaPullPre)*len(pulls)
    #plt.plot(x, pdfPre, color="black", lw=3, ls='--', label="prefit: $\mu = %.2f \pm %.2f, \sigma=%.2f \pm %.2f$" % (meanPullPre, uncMeanPullPre, sigmaPullPre, uncSigmaPullPre))

    ax.set_ylabel('Entries')
    ax.set_xlabel('Pull, p')
    ax.legend(loc=9)
    plt.savefig('closure/pulls_%s_highmass_perbin.pdf' % (save_text))
    plt.savefig('closure/pulls_%s_highmass_perbin.png' % (save_text))
    fig.clf()

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
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)
    hist = hist.Clone(hist_name)
    if type(hist) == ROOT.TH1F or type(hist) == ROOT.TH1D:
        hist.SetDirectory(0)
    rootFile.Close()
    return hist

    
def plot_yields(obs_hist, pred_hist_prefit, pred_hist, year="2016", topology="boosted", extra_text="", kappa_string=""):

    colour_prefit = "#15616d"

    lumi = {"2016": 35.9, "2017": 41.5, "2018": 59.7, "combined": 137.1}

    pred_hist.GetYaxis().SetTitle("Events / category")
    pred_hist.GetYaxis().SetTitleOffset(1.2)
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

    obs_hist.SetBinErrorOption(ROOT.TH1D.kPoisson)


    residuals = obs_hist.Clone("residuals")
    residuals.SetBinErrorOption(ROOT.TH1D.kNormal)
    err_hist = pred_hist.Clone("err_hist")
    pull_hist = residuals.Clone("pull_hist")
    pull_hist_pre = residuals.Clone("pull_hist_pre")

    err_hist.SetFillColor(ROOT.kAzure+2)
    err_hist.SetFillStyle(3345)
    err_hist.SetLineColor(ROOT.kAzure+2)

    err_hist_ratio = err_hist.Clone("err_hist_ratio")

    residuals.Divide(pred_hist)
    residuals.GetXaxis().LabelsOption("v")
    pred_hist.GetXaxis().SetLabelSize(0)
 
    pulls_pre = []
    pulls = []

    for j in range(pred_hist.GetNbinsX()):
        obs = obs_hist.GetBinContent(j+1)
        err_obs = obs_hist.GetBinError(j+1)

        pred = pred_hist.GetBinContent(j+1)
        err_pred = pred_hist.GetBinError(j+1)
        pred_pre = pred_hist_prefit.GetBinContent(j+1)
        err_pred_pre = pred_hist_prefit.GetBinError(j+1)

        ratio = obs/pred

        #pull = (obs-pred)/err_pred
        if obs - err_pred**2 <= 0 :
           pull = 0
        else :
          if obs > 0 :
            #print "predict is " , pred , "error " , err_pred , " obs " , obs
            pull = (obs-pred)/math.sqrt((obs - err_pred**2))
          elif obs == 0 :
            #print "predict is " , pred , "error " , err_pred
            obs = 1.8
            pull = (obs-pred)/math.sqrt((obs - err_pred**2))
        pull_hist.SetBinContent(j+1, pull)
        pull_hist.SetBinError(j+1, 0)

        pull_pre = (obs-pred_pre)/err_pred_pre
        pull_hist_pre.SetBinContent(j+1, pull_pre)
        pull_hist_pre.SetBinError(j+1, 0)
        pulls.append(pull)
        pulls_pre.append(pull_pre)

        err_hist_ratio.SetBinContent(j+1, 1)
        err_hist_ratio.SetBinError(j+1, err_pred/pred)

        residuals.SetBinError(j+1, err_obs*ratio/(obs+1e-9))


    plot_pulls(pulls, pulls_pre, topology+year)

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

    #pad3 = ROOT.TPad("pad3","pad3",0,0.03,1,0.22)
    #pad3.SetTopMargin(0)
    #pad3.SetBottomMargin(0.3)
    #pad3.SetLeftMargin(0.13)
    #pad3.SetRightMargin(0.1)
    #pad3.SetBorderMode(0)

    pad1.Draw()
    pad2.Draw()
    #pad3.Draw()

    text = style.makeText(0.2, 0.88, 0.7, 0.88, topology+", "+extra_text)
    text.SetTextFont(63)
    text.SetTextSize(31)


    text_kappa = style.makeText(0.2, 0.8, 0.7, 0.8, kappa_string)
    text_kappa.SetTextFont(63)
    text_kappa.SetTextSize(31)

    legend = style.makeLegend(0.54, 0.68, 0.79, 0.86)
    legend.AddEntry(obs_hist, 'Observed', 'lp')
    legend.AddEntry(pred_hist, 'Predicted', "lf")
    legend.AddEntry(err_hist, 'Total unc.', "lf")
    legend.AddEntry(pred_hist_prefit, 'Predicted (prefit)', "l")

    legend_signal = style.makeLegend(0.12, 0.76, 0.5, 0.87)
    legend_signal.SetTextSize(19)

    leg_pull = style.makeLegend(0.15, 0.75, 0.4, 0.85)

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

    text_prompt = ROOT.TLatex(obs_hist.GetNbinsX()*1/6-2.5, 0.5*pred_hist.GetMaximum(), "d^{sig}_{xy}<1")
    text_prompt_2 = ROOT.TLatex(obs_hist.GetNbinsX()*4/6-2.5, 0.5*pred_hist.GetMaximum(), "d^{sig}_{xy}<1")

    text_medium = ROOT.TLatex(obs_hist.GetNbinsX()*2/6-2.5, 0.5*pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")
    text_medium_2 = ROOT.TLatex(obs_hist.GetNbinsX()*5/6-2.5, 0.5*pred_hist.GetMaximum(), "1<d^{sig}_{xy}<10")

    text_displaced = ROOT.TLatex(obs_hist.GetNbinsX()*3/6-2.5, 0.5*pred_hist.GetMaximum(), "d^{sig}_{xy}>10")
    text_displaced_2 = ROOT.TLatex(obs_hist.GetNbinsX()*6/6-2.5, 0.5*pred_hist.GetMaximum(), "d^{sig}_{xy}>10")

    text_prompt.SetTextAlign(22)
    text_prompt_2.SetTextAlign(22)
    text_medium.SetTextAlign(22)
    text_medium_2.SetTextAlign(22)
    text_displaced.SetTextAlign(22)
    text_displaced_2.SetTextAlign(22)

    text_OS.SetTextAlign(22)
    text_SS.SetTextAlign(22)

    #legend.AddEntry(err_hist, 'stat. unc', 'f')

    #pad3.cd()
    #pull_hist.GetYaxis().SetRangeUser(-3.5, 3.5)
    #pull_hist.GetYaxis().SetTitleOffset(1)
    #pull_hist.SetMarkerColor(ROOT.kAzure+2)
    #pull_hist.GetYaxis().SetTitle("Pull")
    #pull_hist.GetYaxis().CenterTitle()
    #pull_hist.Draw("P")

    # pull_hist_pre.SetMarkerColor(ROOT.TColor.GetColor(colour_prefit))
    # pull_hist_pre.GetYaxis().SetRangeUser(-3.5, 3.5)
    # pull_hist_pre.SetMarkerStyle(2)
    # pull_hist_pre.SetMarkerSize(2)
    # pull_hist_pre.Draw("P SAME")

    # leg_pull.AddEntry(pull_hist, "post-fit", "p")
    # leg_pull.AddEntry(pull_hist_pre, "pre-fit", "p")

    # lrespull = ROOT.TLine(-0.5, 0., n_bins_total-0.5, 0.)
    # lrespull_plus = ROOT.TLine(-0.5, 2, n_bins_total-0.5, 2.)
    # lrespull_minus = ROOT.TLine(-0.5, -2., n_bins_total-0.5, -2.)

    # lrespull.SetLineWidth(2)
    # lrespull_plus.SetLineStyle(2)
    # lrespull_minus.SetLineStyle(2)

    # lrespull.Draw("SAME")
    
    pad2.cd()
    residuals.GetYaxis().SetTitle("Obs/Pred")
    residuals.GetYaxis().SetTitleOffset(1.2)
    residuals.GetYaxis().SetNdivisions(504)

    #n1 + 100*n2
    residuals.SetMarkerColor(1)
    residuals.SetLineColor(1)
    residuals.GetYaxis().SetRangeUser(0.2, 1.8)
    residuals.Draw("P")
    lres = ROOT.TLine(-0.5, 1., n_bins_total-0.5, 1.)
    lres.Draw("SAME")

    err_hist_ratio.Draw("E2 SAME")

    pad1.cd()
    #pad1.SetLogy()
    #pred_hist.SetMinimum(1e-3)
    #pred_hist_prefit.SetMinimum(1e-3)
    #err_hist.SetMinimum(1e-3)
    #obs_hist.SetMinimum(1e-3)
    #pred_hist.GetYaxis().SetRangeUser(0.3, 300*pred_hist.GetMaximum())

    pred_hist.Draw("HIST")
    err_hist.Draw("E2 SAME")
    pred_hist_prefit.Draw("HIST SAME")
    obs_hist.Draw("HIST P L E SAME")

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

    #leg_pull.Draw("")
    style.makeLumiText(0.6, 0.91, year=year, lumi=lumi[year])

    cv.SaveAs("closure/{}_{}_{}.pdf".format(topology, extra_text, year))
    cv.SaveAs("closure/{}_{}_{}.png".format(topology, extra_text, year))

categories = [
    "mumu_OS", "ee_OS", "mue_OS", "emu_OS", 
    "mumu_SS", "ee_SS", "mue_SS", "emu_SS"
]

n_bins_total = len(categories)*3

category_names = []
category_names_raw = []
nbins = 6

root_file = "fitDiagnostics_highMass_relaxed_perbin.root"

for k, topology in enumerate(["boosted", "resolved"]):

    for i, year in enumerate(["2016", "2017", "2018"]):
        hist_pred = ROOT.TH1D("pred", "pred", n_bins_total, -0.5, n_bins_total-0.5)
        hist_obs = hist_pred.Clone("obs")
        hist_obs.Sumw2(ROOT.kFALSE)

        hist_pred_prefit = hist_pred.Clone("pred_prefit")

        for j, category in enumerate(categories):
            category_names_raw.append(category)
            category_names.append(make_pretty(category))
            kappa_string = ""
            '''
            rootFile = ROOT.TFile(root_file)
            roo_fit_result = rootFile.Get("fit_b")
            kappa_info = roo_fit_result.floatParsFinal().find("uncertainty_"+topology+"_"+year+"_3")
            kappa = 1.15 ** kappa_info.getVal()
            kappa_unc = 1.15 ** (kappa_info.getVal()+kappa_info.getError()) - kappa
            kappa_string = "#kappa = %.2f#pm%.2f" % (kappa, kappa_unc)
            print(kappa_string)
            rootFile.Close()
            '''
            hist_pred_raw = get_hist(root_file, "shapes_fit_b/ch"+str(k+1+i*2)+"_"+category+"_D/total_background")
            hist_pred_prefit_raw = get_hist(root_file, "shapes_prefit/ch"+str(k+1+i*2)+"_"+category+"_D/total_background")
            hist_obs_raw = get_hist(root_file, "shapes_fit_b/ch"+str(k+1+i*2)+"_"+category+"_D/data")

            for idx in [1, 2, 3]:
                pred = hist_pred_raw.GetBinContent(idx)
                pred_up = hist_pred_raw.GetBinErrorUp(idx)
                pred_low = hist_pred_raw.GetBinErrorLow(idx)

                pred_prefit = hist_pred_prefit_raw.GetBinContent(idx)
                pred_up_prefit = hist_pred_prefit_raw.GetBinErrorUp(idx)
                pred_low_prefit = hist_pred_prefit_raw.GetBinErrorLow(idx)

                obs = hist_obs_raw.GetY()
                obs = obs[idx-1]

                if "SS" in category:
                    j_hack = j + 8
                else:
                    j_hack = j

                hist_index = (idx-1)*4 + j_hack + 1
                hist_pred.SetBinContent(hist_index, pred)
                hist_pred.SetBinError(hist_index, pred_up)
                hist_obs.SetBinContent(hist_index, obs)


                hist_pred_prefit.SetBinContent(hist_index, pred_prefit)
                hist_obs.GetXaxis().SetBinLabel(hist_index, shorten(category))
                hist_pred.GetXaxis().SetBinLabel(hist_index, shorten(category))

        plot_yields(hist_obs, hist_pred_prefit, hist_pred, year=year, topology=topology, extra_text="VR1_relaxed", kappa_string=kappa_string)
