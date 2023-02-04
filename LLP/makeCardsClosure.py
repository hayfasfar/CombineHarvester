import CombineHarvester.CombineTools.ch as ch
import ROOT
import os
import math
from multiprocessing import Pool
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

YEARS = ["2016", "2017", "2018"]
NBINS = 24
ZERO_BIN_RATE = 0.001



def get_hist(file_name, hist_name):
    rootFile = ROOT.TFile(file_name)
    hist = rootFile.Get(hist_name)

    try:
        hist = hist.Clone()
        hist.SetDirectory(0)
    except:
        print("Could not read hist from file"+hist_name+file_name)
        return -1
    else:
        rootFile.Close()
        return hist


hist_path = "/home/hep/hsfar/private/limits/histo/abcd/highmass_postunblinding"
#hist_path = "/home/hep/hsfar/private/limits/histo/abcd/lowMass"
#hist_path = "/home/hep/hsfar/private/limits/histo/abcd/highmass_postunblinding_relaxed"
categories = [
    "mumu_OS", 
    "mumu_SS",
    "ee_OS", "ee_SS",
    "mue_OS", "mue_SS",
    "emu_OS", "emu_SS",
]


for topology in ["boosted", "resolved"]:
    for year in YEARS:
        cb = ch.CombineHarvester()
        proc_signal = ch.Process()
        proc_signal.set_process("fake_signal")
        proc_signal.set_era(year)
        proc_signal.set_signal(True)
        text = ""
        for category_name in categories:
            for region in ["A","B","C", "D"]:
                obs = ch.Observation()

                obs_hist = get_hist(
                    os.path.join(hist_path,"{}{}.root".format(year, topology)),
                        "{}/{}".format(category_name, "data"+region)
                )

                obs_hist.SetDirectory(0)
                obs.set_shape(obs_hist, True)
                name = category_name+"_"+region
                obs.set_bin(name)
                obs.set_era(year)
                cb.InsertObservation(obs)

                hist_signal = obs_hist.Clone("signal")
                hist_signal.Scale(1e-10)
                proc_signal.set_bin(name)
                proc_signal.set_shape(hist_signal, True)

                for ibin in range(3):
                    proc = ch.Process()
                    process_name = "bkg_{}_bin{}".format(name, ibin+1)

                    syst_name = "rate_bkg_{}_bin{}_{}_{}".format(name, ibin+1, topology, year)

                    proc.set_process(process_name)
                    proc.set_bin(name)
                    proc.set_era(year)

                    hist = ROOT.TH1F("bkgHist_{}_bin{}".format(name,ibin+1), "", 3, 0.5, 3.5)
                    hist.SetBinContent(ibin+1, 1)
                    hist.SetBinError(ibin+1, 1e-9)
                    hist.SetDirectory(0)

                    proc.set_shape(hist, True)
                    cb.InsertProcess(proc)

                    if region != "D":
                        cb.cp().process([process_name]).bin([name]).AddSyst(cb, syst_name, "rateParam", ch.SystMap("era")([year], 1.))

                        param = cb.GetParameter(syst_name)
                        obs_hist = get_hist(
                            os.path.join(hist_path,"{}{}.root".format(year, topology)),
                                "{}/{}".format(category_name, "data"+region)
                        )

                        content = obs_hist.GetBinContent(ibin+1)
                        err = math.sqrt(obs_hist.GetBinContent(ibin+1))

                        param.set_val(max(0.01, content))
                        param.set_range(0.01, max(4, content+5.*err))
            if region == "D":
                for ibin in range(3):
                    process_name = "bkg_{}_bin{}".format(name, ibin+1)
                    syst_name = "rate_bkg_{}_bin{}_{}_{}".format(name, ibin+1, topology, year)
                    syst_name_A = syst_name.replace("_D", "_A")
                    syst_name_B = syst_name.replace("_D", "_B")
                    syst_name_C = syst_name.replace("_D", "_C")
                    cb.cp().process([process_name]).bin([name]).AddSyst(cb, syst_name, "rateParam",
                        ch.SystMap("era")([year],(
                        "TMath::Max(@0,0.01)*TMath::Max(@1,0.01)/TMath::Max(@2,0.01)",
                        syst_name_B+","+syst_name_C+","+syst_name_A
                            )
                        )
                ) 
                    # 15 % non-closure uncert
                    uncertainty_name = "uncertainty_{}_{}_{}_{}".format(topology, year, ibin+1 , category_name) 
                    #uncertainty_name = "uncertainty_{}_{}_{}".format(topology, year, ibin+1 ) 
                    #uncertainty_name = "uncertainty_{}_{}".format(topology, year) 
                    if topology == "boosted" : 
                       cb.cp().process([process_name]).bin([name]).AddSyst(cb, uncertainty_name, "lnN", ch.SystMap("era")([year], 1.1))
		       if "mumu_OS" in category_name:
		              uncertainty_name_OSMUMU = "unc_boosted_{}_{}_{}".format(year,category_name, ibin+1)
		              cb.cp().process([process_name]).bin([category_name]).AddSyst(cb, uncertainty_name_OSMUMU, "lnN", ch.SystMap("era")([year], 1.25))
                    else :
                       if ibin == 3 :
                          cb.cp().process([process_name]).bin([name]).AddSyst(cb, uncertainty_name, "lnN", ch.SystMap("era")([year], 1.1))
                       else :  
                          cb.cp().process([process_name]).bin([name]).AddSyst(cb, uncertainty_name, "lnN", ch.SystMap("era")([year], 1.15))
                    # for determining uncertainy
                    #cb.cp().process([process_name]).bin([name]).AddSyst(cb, "uncertainty"+topology+year, "rateParam", ch.SystMap("era")([year], 1.))
                    #param = cb.GetParameter("uncertainty"+topology+year)
                    #param.set_range(0, 2)

        cb.InsertProcess(proc_signal)
                

        cb.PrintAll()
        f = ROOT.TFile.Open(os.path.join("closure", year+topology+".root"), "RECREATE")
        cb.cp().WriteDatacard(
            os.path.join("closure", year+topology+".txt"),
            f
        )

        f.Close()
