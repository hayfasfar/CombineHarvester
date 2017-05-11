#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/regex.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/MSSMFull2016/interface/HttSystematics_MSSMRun2.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"
#include "TF1.h"

using namespace std;
using boost::starts_with;
namespace po = boost::program_options;

template <typename T>
void To1Bin(T* proc)
{
    std::unique_ptr<TH1> originalHist = proc->ClonedScaledShape();
    TH1F *hist = new TH1F("hist","hist",1,0,1);
    double err = 0;
    double rate =
        originalHist->IntegralAndError(0, originalHist->GetNbinsX() + 1, err);
    hist->SetDirectory(0);
    hist->SetBinContent(1, rate);
    hist->SetBinError(1, err);
    proc->set_shape(*hist, true);  // True means adjust the process rate to the
                                   // integral of the hist
}

bool BinIsControlRegion(ch::Object const* obj)
{
    return (boost::regex_search(obj->bin(),boost::regex{"_cr$"}) || (obj->channel() == std::string("zmm")));
}

// Useful to have the inverse sometimes too
bool BinIsNotControlRegion(ch::Object const* obj)
{
    return !BinIsControlRegion(obj);
}



int main(int argc, char** argv) {
  // First define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  string SM125= "";
  string mass = "mA";
  string output_folder = "mssm_run2";
  // TODO: option to pick up cards from different dirs depending on channel?
  // ^ Something like line 90?
  string input_folder_em="DESY/";
  string input_folder_et="Imperial/";
  string input_folder_mt="Imperial/";
  string input_folder_tt="Imperial/";
  string input_folder_zmm="KIT/";
  string postfix="";
  bool auto_rebin = false;
  bool manual_rebin = false;
  bool real_data = false;
  int control_region = 0;
  bool check_neg_bins = false;
  bool poisson_bbb = false;
  bool do_w_weighting = true;
  bool zmm_fit = true;
  bool do_jetfakes = false;
  string chan;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
    ("mass,m", po::value<string>(&mass)->default_value(mass))
    ("input_folder_em", po::value<string>(&input_folder_em)->default_value("DESY"))
    ("input_folder_et", po::value<string>(&input_folder_et)->default_value("Imperial"))
    ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value("Imperial"))
    ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value("Imperial"))
    ("input_folder_zmm", po::value<string>(&input_folder_zmm)->default_value("KIT"))
    ("postfix", po::value<string>(&postfix)->default_value(""))
    ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(false))
    ("real_data", po::value<bool>(&real_data)->default_value(false))
    ("manual_rebin", po::value<bool>(&manual_rebin)->default_value(false))
    ("output_folder", po::value<string>(&output_folder)->default_value("mssm_run2"))
    ("SM125,h", po::value<string>(&SM125)->default_value(SM125))
    ("control_region", po::value<int>(&control_region)->default_value(0))
    ("zmm_fit", po::value<bool>(&zmm_fit)->default_value(true))
    ("jetfakes", po::value<bool>(&do_jetfakes)->default_value(false))
    ("channel", po::value<string>(&chan)->default_value("all"))
    ("check_neg_bins", po::value<bool>(&check_neg_bins)->default_value(false))
    ("poisson_bbb", po::value<bool>(&poisson_bbb)->default_value(false))
    ("w_weighting", po::value<bool>(&do_w_weighting)->default_value(false));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  if (do_jetfakes && control_region){
    std::cerr << "\n" << "ERROR: jetfakes and control_region flags cannot be true at the same time. Use --control-region=0 if you want to use the jet fake estimates." << "\n";
    exit (EXIT_FAILURE);
  }

  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  std::map<string, string> input_dir;
  input_dir["em"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMFull2016/shapes/"+input_folder_em+"/";
  input_dir["mt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMFull2016/shapes/"+input_folder_mt+"/";
  input_dir["et"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMFull2016/shapes/"+input_folder_et+"/";
  input_dir["tt"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMFull2016/shapes/"+input_folder_tt+"/";
  input_dir["zmm"]  = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/MSSMFull2016/shapes/"+input_folder_zmm+"/";

  VString chns;
  if ( chan.find("mt") != std::string::npos ) chns.push_back("mt");
  if ( chan.find("et") != std::string::npos ) chns.push_back("et");
  if ( chan.find("em") != std::string::npos ) chns.push_back("em");
  if ( chan.find("tt") != std::string::npos ) chns.push_back("tt");
  if ( chan=="all" ) chns = {"mt","et","tt"/*,"em"*/};

  if (zmm_fit) chns.push_back("zmm");

  RooRealVar mA(mass.c_str(), mass.c_str(), 90., 3200.);
  mA.setConstant(true);
  RooRealVar mH("mH", "mH", 90., 3200.);
  RooRealVar mh("mh", "mh", 90., 3200.);

  map<string, VString> bkg_procs;
  if (do_jetfakes){
    bkg_procs["et"] = {"ZL", "TTT","VVT","ZTT","jetFakes"};
    bkg_procs["mt"] = {"ZL", "TTT","VVT","ZTT","jetFakes"};
    bkg_procs["tt"] = {"ZL", "TTT","VVT","ZTT","jetFakes"};
  }else{
    bkg_procs["et"] = {"W", "QCD", "ZL", "ZJ", "TTT","TTJ", "VVT","VVJ","ZTT"};
    bkg_procs["mt"] = {"W", "QCD", "ZL", "ZJ", "TTT","TTJ", "VVT","VVJ","ZTT"};
    bkg_procs["tt"] = {"W", "QCD", "ZL", "ZJ", "TTT","TTJ", "VVT","VVJ","ZTT"};
  }
  bkg_procs["em"] = {"W", "QCD", "ZLL", "TT", "VV", "ZTT"};
  bkg_procs["zmm"] = {"W", "QCD", "ZL", "ZJ", "TT", "VV", "ZTT"};

  VString SM_procs = {"ggH_SM125", "qqH_SM125", "ZH_SM125", "WminusH_SM125","WplusH_SM125"};

  //Example - could fill this map with hardcoded binning for different
  //categories if manual_rebin is turned on
  map<string, vector<double> > binning;
  binning["et_nobtag_tight"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,3900};
  binning["et_nobtag_loosemt"] = {0,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,700,900,1100,3900};
  binning["et_btag_loosemt"] = {0,60,80,100,120,140,160,180,200,250,300,350,400,500,3900};
  binning["et_btag_tight"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,3900};
  binning["mt_nobtag_tight"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,3900};
  binning["mt_nobtag_loosemt"] = {0,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,1100,3900};
  binning["mt_btag_tight"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,3900};
  binning["mt_btag_loosemt"] = {0,60,80,100,120,140,160,180,200,250,300,350,400,500,3900};
  binning["tt_nobtag"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,900,1100,3900};
  binning["tt_btag"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,900,3900};
  binning["em_nobtag_lowPzeta"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  binning["em_nobtag_mediumPzeta"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  binning["em_nobtag_highPzeta"] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,900,4000};
  binning["em_btag_lowPzeta"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};
  binning["em_btag_mediumPzeta"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};
  binning["em_btag_highPzeta"] = {0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,4000};
/*  binning["zmm_nobtag"] = {60,70,80,90,100,110,120};
  binning["zmm_btag"] = {60,70,80,90,100,110,120};*/

  // Create an empty CombineHarvester instance that will hold all of the
  // datacard configuration and histograms etc.
  ch::CombineHarvester cb;
  // Uncomment this next line to see a *lot* of debug information
  // cb.SetVerbosity(3);

  // Here we will just define two categories for an 8TeV analysis. Each entry in
  // the vector below specifies a bin name and corresponding bin_id.
  //
  map<string,Categories> cats;
  cats["et_13TeV"] = {
    {8, "et_nobtag_tight"},
    {9, "et_btag_tight"},
    {10, "et_nobtag_loosemt"},
    {11, "et_btag_loosemt"},
    };

  cats["em_13TeV"] = {
    {8, "em_nobtag_lowPzeta"},
    {9, "em_btag_lowPzeta"},
    {10, "em_nobtag_mediumPzeta"},
    {11, "em_btag_mediumPzeta"},
    {12, "em_nobtag_highPzeta"},
    {13, "em_btag_highPzeta"}
    };

  cats["tt_13TeV"] = {
    {8, "tt_nobtag"},
    {9, "tt_btag"}
    };

  cats["mt_13TeV"] = {
    {8, "mt_nobtag_tight"},
    {9, "mt_btag_tight"},
    {10, "mt_nobtag_loosemt"},
    {11, "mt_btag_loosemt"},
    };


  cats["zmm_13TeV"] = {
    {8, "zmm_nobtag"},
    {9, "zmm_btag"}
    };

  if (control_region > 0){
      // for each channel use the categories >= 10 for the control regions
      // the control regions are ordered in triples (10,11,12),(13,14,15)...
      for (auto chn : chns){
        // for em or tt do nothing
        if (ch::contains({"em", "tt", "zmm"}, chn)) {
          std::cout << " - Skipping extra control regions for channel " << chn << "\n";
          continue;
        }
        // if ( chn == "em") continue;
        Categories queue;
        int binid = 14;
        for (auto cat:cats[chn+"_13TeV"]){
          if(binid != 23 && binid!=20){ //Fairly shit way of dropping no b-tag loose mT wjets CR
            queue.push_back(make_pair(binid,cat.second+"_wjets_cr"));
          }
          queue.push_back(make_pair(binid+1,cat.second+"_qcd_cr"));
          if(binid != 23 && binid!=20){ //Fairly shit way of dropping no b-tag loose mT wjets ss CR
            queue.push_back(make_pair(binid+2,cat.second+"_wjets_ss_cr"));
          }
          binid += 3;
        }
        cats[chn+"_13TeV"].insert(cats[chn+"_13TeV"].end(),queue.begin(),queue.end());
      }
  }

  vector<string> masses = {"90","100","110","120","130","140","160","180", "200", "250", "350", "400", "450", "500"/*, "600"*/, "700", "800", "900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200"};

  map<string, VString> signal_types = {
    {"ggH", {"ggh_htautau", "ggH_Htautau", "ggA_Atautau"}},
    {"bbH", {"bbh_htautau", "bbH_Htautau", "bbA_Atautau"}}
  };
  if(mass=="MH"){
    signal_types = {
      {"ggH", {"ggH"}},
      {"bbH", {"bbH"}}
    };
  }
    vector<string> sig_procs = {"ggH","bbH"};
  for(auto chn : chns){
    cb.AddObservations({"*"}, {"htt"}, {"13TeV"}, {chn}, cats[chn+"_13TeV"]);

    cb.AddProcesses({"*"}, {"htt"}, {"13TeV"}, {chn}, bkg_procs[chn], cats[chn+"_13TeV"], false);

    cb.AddProcesses(masses, {"htt"}, {"13TeV"}, {chn}, signal_types["ggH"], cats[chn+"_13TeV"], true);
    cb.AddProcesses(masses, {"htt"}, {"13TeV"}, {chn}, signal_types["bbH"], cats[chn+"_13TeV"], true);
    if(SM125==string("bkg_SM125") && chn!="zmm") cb.AddProcesses({"*"}, {"htt"}, {"13TeV"}, {chn}, SM_procs, cats[chn+"_13TeV"], false);
    if(SM125==string("signal_SM125") && chn!="zmm") cb.AddProcesses({"*"}, {"htt"}, {"13TeV"}, {chn}, SM_procs, cats[chn+"_13TeV"], true);
    }
  if ((control_region > 0) || zmm_fit){
      // Since we now account for QCD in the high mT region we only
      // need to filter signal processes
      cb.FilterAll([](ch::Object const* obj) {
              return (BinIsControlRegion(obj) && obj->signal());
              });
  }




  ch::AddMSSMRun2Systematics(cb, control_region, zmm_fit);
  //! [part7]
  for (string chn:chns){
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + "htt_"+chn+".inputs-mssm-13TeV"+postfix+".root",
        "$BIN/$PROCESS",
        "$BIN/$PROCESS_$SYSTEMATIC");
    if(SM125==string("signal_SM125")) cb.cp().channel({chn}).process(SM_procs).ExtractShapes(
         input_dir[chn] + "htt_"+chn+".inputs-mssm-13TeV"+postfix+".root",
         "$BIN/$PROCESS",
         "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(signal_types["ggH"]).ExtractShapes(
        input_dir[chn] + "htt_"+chn+".inputs-mssm-13TeV"+postfix+".root",
        "$BIN/ggH$MASS",
        "$BIN/ggH$MASS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(signal_types["bbH"]).ExtractShapes(
        input_dir[chn] + "htt_"+chn+".inputs-mssm-13TeV"+postfix+".root",
        "$BIN/bbH$MASS",
        "$BIN/bbH$MASS_$SYSTEMATIC");
  }


 //Now delete processes with 0 yield
 cb.FilterProcs([&](ch::Process *p) {
  bool null_yield = !(p->rate() > 0. || BinIsControlRegion(p));
  if (null_yield){
     std::cout << "[Null yield] Removing process with null yield: \n ";
     std::cout << ch::Process::PrintHeader << *p << "\n"; 
     cb.FilterSysts([&](ch::Systematic *s){
       bool remove_syst = (MatchingProcess(*p,*s));
       return remove_syst;
    });
  }
  return null_yield;
 });

    //Scaling QCD in em btag by 1.45/2.2 
   /*cb.cp().process({"QCD"}).channel({"em"}).bin_id({9}).ForEachProc([&](ch::Process *proc){
         proc->set_rate(proc->rate()*(1.45/2.2));
   });*/

  if(SM125!=string("")) {
     cb.cp().process(SM_procs).ForEachProc([&](ch::Process * proc) {
       proc->set_rate(proc->rate()*0.0627);
      });
     cb.cp().process({"ggH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*48.58);
     });
     cb.cp().process({"qqH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*3.782);
     });
     cb.cp().process({"WplusH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.8400);
     });
     cb.cp().process({"WminusH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.5328);
     });
     cb.cp().process({"ZH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.8839);
     });
/*     cb.cp().process({"TTH_SM125"}).ForEachProc([&](ch::Process *proc){
       proc->set_rate(proc->rate()*0.5071);
     });*/
  }

  // And convert any shapes in the CRs to lnN:
  // Convert all shapes to lnN at this stage
  cb.cp().FilterSysts(BinIsNotControlRegion).syst_type({"shape"}).ForEachSyst([](ch::Systematic *sys) {
    sys->set_type("lnN");
  });

   //Replacing observation with the sum of the backgrounds (asimov) - nice to ensure blinding
    auto bins = cb.cp().bin_set();
    // For control region bins data (should) = sum of bkgs already
    // useful to be able to check this, so don't do the replacement
    // for these
  if(!real_data){
      for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
          std::cout << " - Replacing data with asimov in bin " << b << "\n";
          cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
            obs->set_shape(cb.cp().bin({b}).backgrounds().GetShape(), true);
          });
        }
  }


  // Merge to one bin for control region bins
  cb.cp().FilterAll(BinIsNotControlRegion).ForEachProc(To1Bin<ch::Process>);
  cb.cp().FilterAll(BinIsNotControlRegion).ForEachObs(To1Bin<ch::Observation>);

  // Rebinning
  // --------------------
  // Keep track of shapes before and after rebinning for comparison
  // and for checking the effect of negative bin contents
  std::map<std::string, TH1F> before_rebin;
  std::map<std::string, TH1F> after_rebin;
  std::map<std::string, TH1F> after_rebin_neg;
  if (check_neg_bins) {
    for (auto b : bins) {
      before_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
    }
  }


  auto rebin = ch::AutoRebin()
    .SetBinThreshold(0.)
    .SetBinUncertFraction(0.9)
    .SetRebinMode(1)
    .SetPerformRebin(true)
    .SetVerbosity(1);
  if(auto_rebin) rebin.Rebin(cb, cb);

  if(manual_rebin) {
    for(auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
      std::cout << "Rebinning by hand for bin: " << b <<  std::endl;
      cb.cp().bin({b}).VariableRebin(binning[b]);
    }
  }

  if (check_neg_bins) {
    for (auto b : bins) {
      after_rebin[b] = cb.cp().bin({b}).backgrounds().GetShape();
      // std::cout << "Bin: " << b << " (before)\n";
      // before_rebin[b].Print("range");
      // std::cout << "Bin: " << b << " (after)\n";
      // after_rebin[b].Print("range");
      // Build a sum-of-bkgs TH1 that doesn't truncate the negative yields
      // like the CH GetShape does
      for (auto p : cb.cp().bin({b}).backgrounds().process_set()) {
        TH1F proc_hist;
        cb.cp().bin({b}).process({p}).ForEachProc([&](ch::Process *proc) {
          proc_hist = proc->ShapeAsTH1F();
          proc_hist.Scale(proc->no_norm_rate());
          for (int i = 1; i <= proc_hist.GetNbinsX(); ++i) {
            if (proc_hist.GetBinContent(i) < 0.) {
              std::cout << p << " bin " << i << ": " << proc_hist.GetBinContent(i) << "\n";
            }
          }
        });
        if (after_rebin_neg.count(b)) {
          after_rebin_neg[b].Add(&proc_hist);
        } else {
          after_rebin_neg[b] = proc_hist;
        }
      }
      std::cout << "Bin: " << b << "\n";
      for (int i = 1; i <= after_rebin[b].GetNbinsX(); ++i) {
        double offset = after_rebin[b].GetBinContent(i) - after_rebin_neg[b].GetBinContent(i);
        double offset_by_yield = offset / after_rebin[b].GetBinContent(i);
        double offset_by_err = offset / after_rebin[b].GetBinError(i);
        printf("%-2i offset %-10.4f tot %-10.4f err %-10.4f off/tot %-10.4f off/err %-10.4f\n", i , offset, after_rebin[b].GetBinContent(i), after_rebin[b].GetBinError(i), offset_by_yield, offset_by_err);
      }
    }
  }

  // Uncomment this to inject 1 obs event in the last bin of every signal-region
  // category
  // if(!real_data){
  //     for (auto b : cb.cp().FilterAll(BinIsControlRegion).bin_set()) {
  //       std::cout << " - Adjusting data in bin " << b << "\n";
  //         cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
  //           TH1F new_obs = cb.cp().bin({b}).GetObservedShape();
  //           new_obs.SetBinContent(new_obs.GetNbinsX(), 1.);
  //           new_obs.Print("range");
  //           obs->set_shape(new_obs, true);
  //         });
  //       }
  // }

  // At this point we can fix the negative bins
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      std::cout << "[Negative bins] Fixing negative bins for " << p->bin()
                << "," << p->process() << "\n";
      // std::cout << "[Negative bins] Before:\n";
      // p->shape()->Print("range");
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      p->set_shape(std::move(newhist), false);
      // std::cout << "[Negative bins] After:\n";
      // p->shape()->Print("range");
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos) return;
    if (ch::HasNegativeBins(s->shape_u()) || ch::HasNegativeBins(s->shape_d())) {
      std::cout << "[Negative bins] Fixing negative bins for syst" << s->bin()
                << "," << s->process() << "," << s->name() << "\n";
      // std::cout << "[Negative bins] Before:\n";
      // s->shape_u()->Print("range");
      // s->shape_d()->Print("range");
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      // Set the new shape but do not change the rate, we want the rate to still
      // reflect the total integral of the events
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
      // std::cout << "[Negative bins] After:\n";
      // s->shape_u()->Print("range");
      // s->shape_d()->Print("range");
    }
  });

  cout << "Generating bbb uncertainties...";
  auto bbb = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(true)
    .SetMergeSaturatedBins(false)
    .SetPoissonErrors(poisson_bbb);
  for (auto chn : chns) {
    std::cout << " - Doing bbb for channel " << chn << "\n";
    bbb.MergeAndAdd(cb.cp().channel({chn}).process({"ZTT", "QCD", "W", "ZJ", "ZL", "TT", "VV", "Ztt", "ttbar", "EWK", "Fakes", "ZMM", "TTT","TTJ","VVT","VVJ", "WJets", "Dibosons"}).FilterAll([](ch::Object const* obj) {
                return BinIsControlRegion(obj);
                }), cb);
  }
  // And now do bbb for the control region with a slightly different config:
  auto bbb_ctl = ch::BinByBinFactory()
    .SetPattern("CMS_$ANALYSIS_$BIN_$ERA_$PROCESS_bin_$#")
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(false)  // contrary to signal region, bbb *should* change yield here
    .SetVerbosity(1);
  // Will merge but only for non W and QCD processes, to be on the safe side
  bbb_ctl.MergeBinErrors(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion));
  bbb_ctl.AddBinByBin(cb.cp().process({"QCD", "W"}, false).FilterProcs(BinIsNotControlRegion), cb);
  cout << " done\n";

  //Switch JES over to lnN:
  /*cb.cp().syst_name({"CMS_scale_j_13TeV"}).ForEachSyst([](ch::Systematic *sys) { sys->set_type("lnN");});
  cb.cp().syst_name({"CMS_scale_b_13TeV"}).ForEachSyst([](ch::Systematic *sys) { sys->set_type("lnN");});
  cb.cp().syst_name({"CMS_fake_b_13TeV"}).ForEachSyst([](ch::Systematic *sys) { sys->set_type("lnN");});*/

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  // which is commonly used in the htt analyses
  ch::SetStandardBinNames(cb);
  //! [part8]


  //! [part9]
  // First we generate a set of bin names:
  RooWorkspace ws("htt", "htt");

  TFile demo("htt_mssm_demo.root", "RECREATE");

  bool do_morphing = true;
  map<string, RooAbsReal *> mass_var = {
    {"ggh_htautau", &mh}, {"ggH_Htautau", &mH}, {"ggA_Atautau", &mA},
    {"bbh_htautau", &mh}, {"bbH_Htautau", &mH}, {"bbA_Atautau", &mA}
  };
  if(mass=="MH"){
    mass_var = {
      {"ggH", &mA},
      {"bbH", &mA}
    };
  }
  if (do_morphing) {
    auto bins = cb.bin_set();
    for (auto b : bins) {
      auto procs = cb.cp().bin({b}).process(ch::JoinStr({signal_types["ggH"], signal_types["bbH"]})).process_set();
      for (auto p : procs) {
        ch::BuildRooMorphing(ws, cb, b, p, *(mass_var[p]),
                             "norm", true, false, false, &demo);
      }
    }
  }
  demo.Close();
  cb.AddWorkspace(ws);
  cb.cp().process(ch::JoinStr({signal_types["ggH"], signal_types["bbH"]})).ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");


 //Write out datacards. Naming convention important for rest of workflow. We
 //make one directory per chn-cat, one per chn and cmb. In this code we only
 //store the individual datacards for each directory to be combined later, but
 //note that it's also possible to write out the full combined card with CH
  string output_prefix = "output/";
  if(output_folder.compare(0,1,"/") == 0) output_prefix="";
  ch::CardWriter writer(output_prefix + output_folder + "/$TAG/$BIN.txt",
                        output_prefix + output_folder + "/$TAG/htt_input.root");
  // We're not using mass as an identifier - which we need to tell the CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});
  writer.SetVerbosity(1);

  writer.WriteCards("cmb", cb);
  for (auto chn : chns) {
    if(chn == std::string("zmm"))
    {
        continue;
    }
    // per-channel
    writer.WriteCards(chn, cb.cp().channel({chn, "zmm"}));
    // And per-channel-category
    writer.WriteCards("htt_"+chn+"_8_13TeV", cb.cp().channel({chn,"zmm"}).attr({"tight","high"},"mtsel").attr({"nobtag"},"cat"));
    writer.WriteCards("htt_"+chn+"_9_13TeV", cb.cp().channel({chn,"zmm"}).attr({"tight","high"},"mtsel").attr({"btag"},"cat"));
    if(chn != std::string("tt")){
      writer.WriteCards("htt_"+chn+"_10_13TeV", cb.cp().channel({chn,"zmm"}).attr({"loose","high"},"mtsel").attr({"nobtag"},"cat"));
      writer.WriteCards("htt_"+chn+"_11_13TeV", cb.cp().channel({chn,"zmm"}).attr({"loose","high"},"mtsel").attr({"btag"},"cat"));
    }
    if(chn == std::string("em")){
        writer.WriteCards("htt_"+chn+"_12_13TeV", cb.cp().channel({chn,"zmm"}).bin_id({12}));
        writer.WriteCards("htt_"+chn+"_13_13TeV", cb.cp().channel({chn,"zmm"}).bin_id({13}));
    }
  }
  // For btag/nobtag areas want to include control regions. This will
  // work even if the extra categories aren't there.
  writer.WriteCards("htt_cmb_btag_13TeV", cb.cp().attr({"btag"},"cat"));
  writer.WriteCards("htt_cmb_nobtag_13TeV", cb.cp().attr({"nobtag"},"cat"));

  cb.PrintAll();
  cout << " done\n";
}
