import CombineHarvester.CombineTools.ch as ch
import ROOT
import os
import errno
import argparse
import math
from multiprocessing import Pool
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

YEARS = ["2016", "2017", "2018"]
COUPLINGS = [1, 2, 7, 12, 47, 52] 
NBINS = 48
ZERO_BIN_RATE = 0.001
PSEUDO_DATA = 0
NWORKERS = 16


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

def make_sge_script(procs):
    submit_file = open("runCombine.sh","w")
    submit_file.write('''#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=00:30:00
#$ -e log/log.$TASK_ID.err
#$ -o log/log.$TASK_ID.out
#$ -t 1-'''+str(4*n_job)+'''
hostname
date
source ~/.cms_setup.sh
eval `scramv1 runtime -sh`
''')

    submit_file.write("JOBS=(\n")

    for proc in procs:
        for COUPLING in COUPLINGS:
            for YEAR in YEARS:
                if status_dict[YEAR][COUPLING][proc]:
                    path = os.path.join('$PWD/cards/{}/coupling_{}/{}'.format(YEAR, COUPLING, proc))
                    submit_file.write(" \"")
                    submit_file.write('''combineTool.py -M AsymptoticLimits --run expected --cminPreScan --cminPreFit 1 --rAbsAcc 0.000001 -d %s/out.txt --there -n HNL --mass %i''' % (path, COUPLING))
                    submit_file.write("\"")
                    submit_file.write("\n")

    
            path_2016 = os.path.expandvars(os.path.join('$PWD/cards/{}/coupling_{}/{}/'.format(2016, COUPLING, proc)))
            path_2017 = os.path.expandvars(os.path.join('$PWD/cards/{}/coupling_{}/{}/'.format(2017, COUPLING, proc)))
            path_2018 = os.path.expandvars(os.path.join('$PWD/cards/{}/coupling_{}/{}/'.format(2018, COUPLING, proc)))
            path_combined = os.path.join('cards/{}/coupling_{}/{}/'.format("combined", COUPLING, proc))

            status_2016 = status_dict["2016"][COUPLING][proc]
            status_2017 = status_dict["2017"][COUPLING][proc]
            status_2018 = status_dict["2018"][COUPLING][proc]

            if status_2016 or status_2017 or status_2018:
                if not os.path.exists(path_combined):
                    os.makedirs(path_combined)

                combine_string = ""
                if status_2016:
                    combine_string += path_2016+"out.txt "
                if status_2017:
                    combine_string += path_2017+"out.txt "
                if status_2018:
                    combine_string += path_2018+"out.txt "

                submit_file.write(" \"")
                submit_file.write("combineCards.py "+combine_string+" >> " +path_combined+"out.txt ")
                submit_file.write('''&& combineTool.py -M AsymptoticLimits --run expected --cminPreScan --cminPreFit 1 --rAbsAcc 0.000001 -d %sout.txt --there -n HNL --mass %i''' % (path_combined, COUPLING))
                submit_file.write("\"")
                submit_file.write("\n")


    submit_file.write(")")
    submit_file.write("\n")
    submit_file.write("echo ${JOBS[$SGE_TASK_ID-1]}")
    submit_file.write("\n")
    submit_file.write('''eval "${JOBS[$SGE_TASK_ID-1]}"''')
    submit_file.close()

def worker(proc):
    print("Started " + proc + ", "+str(hnl_sample_list.index(proc))+"/"+str(n_job/len(COUPLINGS)))
    for YEAR in YEARS:
        for COUPLING in COUPLINGS:
            path = os.path.join('cards/{}/coupling_{}/{}'.format(YEAR, COUPLING, proc))
            make_datacard(category_pairs, category_pairs_signal, proc, path, coupling=COUPLING, year=YEAR)

    print("Finished " + str(hnl_sample_list.index(proc))+"/"+str(n_job/len(COUPLINGS)))

# make a datacard for a single HNL mass/coupling scenario
def make_datacard(cats, cats_signal, signal_name, output_path, coupling=12, year="2016"):
    cb = ch.CombineHarvester()
    bkgs_mc = []
    bkgs_abcd = ["wjets", "dyjets", "qcd", "vgamma", "topbkg"]
    signal = ["HNL"]

    cb.AddProcesses(era=[year], procs=bkgs_mc, bin=cats, signal=False)
    cb.AddProcesses(era=[year], procs=signal, bin=cats, signal=True) 

    systematics_uncorrelated = ["pu", "unclEn", "jesTotal", "jer", "trigger", "tight_muon_iso", "tight_muon_id", "tight_electron_id", "tight_electron_reco", "loose_electron_reco", "displaced_track", "scale", "pdf"]
    systematics_correlated = []

    lumi_uncertainty = {"2016": 1.025, "2017": 1.023, "2018": 1.025}

    for syst in systematics_correlated:
        cb.cp().signals().AddSyst( cb, syst, "shape", ch.SystMap()(1.0) )
    
    for syst in systematics_uncorrelated:
        cb.cp().signals().AddSyst( cb, syst+"$ERA", "shape", ch.SystMap()(1.0) )

    cb.cp().AddSyst(cb, "lumi_$ERA", "lnN", ch.SystMap("era")([year], lumi_uncertainty[year]))

    cb.cp().signals().ExtractShapes(
            "{}/{}_{}.root".format(hist_path, signal_name, year),
            "$BIN/$PROCESS_coupling_{}".format(coupling),
            "$BIN/$PROCESS_coupling_{}_$SYSTEMATIC".format(coupling)
            )

    cb.cp().backgrounds().ExtractShapes(
            "{}/{}.root".format(hist_path, year),
            "$BIN/$PROCESS",
            "$BIN/$PROCESS_$SYSTEMATIC"
            )
    
    for _, category_name in cats:
        obs = ch.Observation()
        # pseudo-data = sum of MC
        if PSEUDO_DATA:
            for i, bkg in enumerate(bkgs_mc+bkgs_abcd):
                bkg_obs_hist = get_hist(
                    os.path.join(hist_path,"{}.root".format(year)),
                        "{}/{}".format(category_name, bkg)
                )
                if i == 0:
                    obs_sum_hist = bkg_obs_hist.Clone(category_name)
                else:
                    obs_sum_hist.Add(bkg_obs_hist)
        else: #data
            obs_sum_hist = get_hist(
                os.path.join(hist_path,"{}.root".format(year)),
                    "{}/{}".format(category_name, "data")
            )

        obs_sum_hist.SetDirectory(0)
        obs.set_shape(obs_sum_hist, True)
        obs.set_bin(category_name)
        obs.set_era(year)
        cb.InsertObservation(obs)

    bkg_processes = []

    for _, category_name in cats_signal:

        nbins = 6
        bin_min = 0.5
        bin_max = nbins+0.5


        for region in ["A","B","C", "D"]:
            name = category_name.replace("_D", "_"+region)
            for ibin in range(nbins):
                # Dummy histogram per bin , scale by rate parameters
                proc = ch.Process()
                process_name = "bkg_{}_bin{}".format(name, ibin+1)
                syst_name = "rate_bkg_{}_bin{}_{}".format(name, ibin+1, year)

                proc.set_process(process_name)
                proc.set_bin(name)
                proc.set_era(year)

                hist = ROOT.TH1F("bkgHist_{}_bin{}".format(name,ibin+1), "", nbins, bin_min, bin_max)
                hist.SetBinContent(ibin+1, 1)
                hist.SetBinError(ibin+1, 1e-9)
                hist.SetDirectory(0)

                proc.set_shape(hist, True)
                cb.InsertProcess(proc)

                if "_D" in process_name:
                    bkg_processes.append(process_name)

                else:
                    cb.cp().process([process_name]).bin([name]).AddSyst(cb, syst_name, "rateParam", ch.SystMap("era")([year], 1.))

                    param = cb.GetParameter(syst_name)

                    # Set ABCD rate parameter values from MC or data:
                    if PSEUDO_DATA:
                        for i, bkg in enumerate(bkgs_abcd):
                            bkg_hist = get_hist(
                                os.path.join(hist_path,"{}.root".format(year)),
                                    "{}/{}".format(name, bkg)
                            )
                            if i == 0:
                                bkg_hist_sum = bkg_hist.Clone(name)
                            else:
                                bkg_hist_sum.Add(bkg_hist)
                    else:
                        bkg_hist_sum = get_hist(
                            os.path.join(hist_path,"{}.root".format(year)),
                                "{}/{}".format(name, "data")
                        )

                    content = bkg_hist_sum.GetBinContent(ibin+1)
                    err = math.sqrt(bkg_hist_sum.GetBinContent(ibin+1))

                    param.set_val(max(0.01, content))
                    param.set_range(0.01, max(4, content+5.*err))
 
        # ABCD D = B*C/A
        for ibin in range(nbins):
            process_name = "bkg_{}_bin{}".format(category_name, ibin+1)
            syst_name = "rate_bkg_{}_bin{}_{}".format(category_name, ibin+1, year)
            syst_name_A = syst_name.replace("_D", "_A")
            syst_name_B = syst_name.replace("_D", "_B")
            syst_name_C = syst_name.replace("_D", "_C")
            cb.cp().process([process_name]).bin([category_name]).AddSyst(cb, syst_name, "rateParam",
                ch.SystMap("era")([year],(
                "TMath::Max(@0,0.01)*TMath::Max(@1,0.01)/TMath::Max(@2,0.01)",
                syst_name_B+","+syst_name_C+","+syst_name_A
                    )
                )
            ) 
            if ibin in [0, 1, 2]:
                uncertainty_name = "unc_boosted_{}".format(year)
            else:
                uncertainty_name = "unc_resolved_{}".format(year)
            cb.cp().process([process_name]).bin([category_name]).AddSyst(cb, uncertainty_name, "lnN", ch.SystMap("era")([year], 1.1))

    #if abs(cb.cp().signals().GetRate() - ZERO_BIN_RATE*NBINS)/ZERO_BIN_RATE*NBINS < 0.1:
        #print("Zero signal in histogram for the signal sample. Check your cuts! Skipping process: "+ signal_name +", coupling: "+ #str(coupling)+", year: "+str(year))
        #return False

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    f = ROOT.TFile.Open(os.path.join(output_path, "out.root"), "RECREATE")
    cb.cp().WriteDatacard(
        os.path.join(output_path, "out.txt"),
        f
    )
    f.Close()
    return True

parser = argparse.ArgumentParser()
parser.add_argument("--path", default="/vols/cms/vc1117/AN-19-207/histo/limits/hists_merged")
args = parser.parse_args()
hist_path = args.path

categories = [
    "mumu_OS", "mumu_SS",
    "ee_OS", "ee_SS",
    "mue_OS", "mue_SS",
    "emu_OS", "emu_SS",
]

hnl_sample_list = []
n_job=0
for proc in os.listdir(hist_path):
    if "HNL" not in proc:
        continue
    if "_pt20_" not in proc and proc.replace("_all_", "_pt20_") in os.listdir(hist_path):
        print("Skipping "+proc+", higher stat. sample exists")
        continue
    if "_2016" in proc:
        hnl_sample_list.append(proc.replace("_2016.root", ""))
        n_job += len(COUPLINGS)


# Make cards in parallel 
# Delete the old directory
try:
    os.mkdir("cards")
except OSError as e:
    if e.errno == errno.EEXIST:
        raise OSError('Directory "cards" already created.')


# Count the number of jobs
n_categories = len(categories)
category_pairs = []
category_pairs_signal = [] # So far signal only in region D (to save CPU time)

for index1, category_name in enumerate(categories):
    for index2, region in enumerate(["A", "B", "C", "D"]):
        index = index1*4 + index2 
        name = category_name+"_"+region
        pair=(index, name)
        if region == "D":
            category_pairs_signal.append(pair)
        category_pairs.append(pair)
pool = Pool(NWORKERS)
pool.map(worker, hnl_sample_list)


status_dict = {}
for year in YEARS:
    status_dict[year] = {}
    for coupling in COUPLINGS:
        status_dict[year][coupling] = {}
        for proc in hnl_sample_list:
            if os.path.exists(os.path.join('cards', year, 'coupling_'+str(coupling), proc)):
                status_dict[year][coupling][proc] = True
            else:
                status_dict[year][coupling][proc] = False
# Make sge submission script
make_sge_script(hnl_sample_list)
