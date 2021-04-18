import matplotlib
matplotlib.use('Agg') 
import style
import ROOT
ROOT.gROOT.SetBatch(1)
import json
from array import array
from scipy import interpolate
import numpy as np
import pandas as pd
import os

def get_mu(theory, exp):
    if np.isnan(theory) or np.isnan(exp):
        mu = 4.
    else:
        mu = exp-theory
    mu = np.clip(mu, -6., 4.)
    return mu

def parse_lookup_table(f, lookup_table):

    # parse lookup table
    proc = f.replace(year+"_", "").replace("limits_", "").replace(".json", "")
    lu_infos = lookup_table[proc]['weights'][str(int(scenario))]
    xsec = lu_infos['xsec']['nominal']*K_FACTOR
    coupling = lu_infos['couplings']['Ve']+lu_infos['couplings']['Vmu']+lu_infos['couplings']['Vtau']
    if coupling not in (2, 12, 67):
        coupling = coupling/2
    coupling = coupling ** 2
    mass = lookup_table[proc]['mass']

    return mass, coupling, xsec

def interpolate_point(coupling_points, mu_points):
    for i, (coupling, mu) in enumerate(zip(coupling_points, mu_points)):
        if i < len(mu_points)-1 and mu > 0 and mu_points[i+1] < 0:
            return (coupling+coupling_points[i+1])/2
    return 

def parse_limit_json(f, scenario=12):

    # limit json aggreggate 
    with open(os.path.join(json_path, f)) as json_file:
        xsec_dict = json.load(json_file)
    if str(scenario) not in xsec_dict.keys():
        raise ValueError
    xsec_dict = xsec_dict[str(scenario)]
    # means something failed with combine (shouldn't happend!)
    if "exp0" not in xsec_dict or "exp+1" not in xsec_dict or "exp+2" not in xsec_dict or "exp-1" not in xsec_dict or "exp-2" not in xsec_dict:
        return None

    return xsec_dict

def get_graph(filename, histname):
    root_file = ROOT.TFile(filename)
    hist = root_file.Get(histname)
    root_file.Close()
    return hist

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

json_path = "jsons_old"

K_FACTOR = 1.1

lumi = {"2016": 35.9, "2017": 41.5, "2018": 59.7, "combined": 137.1}
years = ["2016", "2017", "2018", "combined"]

coupling_dict = {}
coupling_dict[2.0] = ["ee", "U_{e} : U_{#mu} : U_{#tau} = 1 : 0 : 0"]
coupling_dict[7.0] = ["emu", "U_{e} : U_{#mu} : U_{#tau} = 1 : 1 : 0"]
coupling_dict[12.0] = ["mumu", "U_{e} : U_{#mu} : U_{#tau} = 0 : 1 : 0"]
coupling_dict[47.0] = ["etau", "U_{e} : U_{#mu} : U_{#tau} = 1 : 0 : 1"]
coupling_dict[52.0] = ["mutau", "U_{e} : U_{#mu} : U_{#tau} = 0 : 1 : 1"]
#coupling_dict[67.0] = ["tautau", "U_{e} : U_{#mu} : U_{#tau} = 0 : 0 : 1"]

n_bins = 200

mass_range = np.geomspace(1., 24., num=n_bins)
log_coupling_range = np.linspace(-8, -1., num=n_bins)
coupling_range = np.power(10, log_coupling_range)

with open("/vols/cms/LLP/gridpackLookupTable.json") as lookup_table_file:
    lookup_table = json.load(lookup_table_file)

for year in years:
    print(year)
    for scenario in coupling_dict.keys():
        print("Analyzing coupling scenario: "+str(scenario))
        coupling_text = coupling_dict[scenario][0]
        coupling_title = coupling_dict[scenario][1]

        for hnl_type in ["majorana", "dirac"]:
            print(hnl_type)
            # arrays to store mass, V2 and sigma/sigma_th values
            masses = []
            couplings = []
            sigma_dict = {}
            for exp_var in ["exp0", "exp+1", "exp+2", "exp-1", "exp-2", "theory"]:
                sigma_dict[exp_var] = []

            for f in os.listdir(json_path):
                if year not in f or ".json" not in f or hnl_type not in f:
                    continue

                xsec_dict = parse_limit_json(f, scenario)
                mass, coupling, xsec = parse_lookup_table(f, lookup_table)
                sigma_dict["theory"].append(xsec)

                for exp_var in ["exp0", "exp+1", "exp+2", "exp-1", "exp-2"]:
                    sigma_dict[exp_var].append(xsec_dict[exp_var])

                masses.append(mass)
                couplings.append(coupling)

            df = pd.DataFrame(list(zip(masses, couplings, sigma_dict["theory"], sigma_dict["exp0"], sigma_dict["exp+1"],
             sigma_dict["exp+2"], sigma_dict["exp-1"], sigma_dict["exp-2"])), 
                        columns =['mass', 'coupling', 'theory', 'exp0', 'exp+1', 'exp+2', 'exp-1', 'exp-2'])
            npoints = len(df)
            mass_coupling_pair = np.array([df['mass'], np.log10(df['coupling'])]).T
            log_theory_points = np.log10(np.array(df['theory']))
            log_expected_points = np.log10(np.array(df['exp0']))
            log_expected_points_plus = np.log10(np.array(df['exp+1']))
            log_expected_points_minus = np.log10(np.array(df['exp-1']))

            n_bins = 200
            grid_x, grid_y = np.meshgrid(mass_range, log_coupling_range, indexing='ij')
            fit_method = 'cubic'

            results_theory = interpolate.griddata(mass_coupling_pair, log_theory_points, (grid_x, grid_y), method=fit_method)
            results = interpolate.griddata(mass_coupling_pair, log_expected_points, (grid_x, grid_y), method=fit_method)
            results_plus = interpolate.griddata(mass_coupling_pair, log_expected_points_plus, (grid_x, grid_y), method=fit_method)
            results_minus = interpolate.griddata(mass_coupling_pair, log_expected_points_minus, (grid_x, grid_y), method=fit_method)

            hist_mu = ROOT.TH2D("mu"+hnl_type+str(scenario), "mu", n_bins-1, mass_range, n_bins-1, coupling_range)
            crossing_points = []
            masses_at_crossing_points = []
            errors_up = []
            errors_down = []

            for i in range(n_bins):
                mass = mass_range[i]
                coupling_values = []
                mu_values = []
                mu_plus_values = []
                mu_minus_values = []

                for j in range(n_bins):
                    coupling_values.append(log_coupling_range[j])
                    exp = results[i, j]
                    exp_plus = results_plus[i, j]
                    exp_minus = results_minus[i, j]
                    theory = results_theory[i, j]

                    mu_values.append(get_mu(theory, exp))
                    mu_plus_values.append(get_mu(theory, exp_plus))
                    mu_minus_values.append(get_mu(theory, exp_minus))

                    hist_mu.SetBinContent(i+1, j+1, np.power(10, get_mu(theory, exp)))
                crossing_point = interpolate_point(coupling_values, mu_values)
                crossing_point_plus = interpolate_point(coupling_values, mu_plus_values) 
                crossing_point_minus = interpolate_point(coupling_values, mu_minus_values)

                if crossing_point and crossing_point_plus and crossing_point_minus:
                    crossing_points.append(np.power(10, crossing_point))
                    errors_down.append(np.power(10, crossing_point_plus)-np.power(10,crossing_point))
                    errors_up.append(np.power(10, crossing_point)-np.power(10, crossing_point_minus))
                    masses_at_crossing_points.append(mass)

            x_errors = np.zeros(len(crossing_points))  
            graph = ROOT.TGraphAsymmErrors(len(crossing_points), array('d', masses_at_crossing_points), array('d', crossing_points), x_errors, x_errors, array('d', errors_up), array('d', errors_down))

            cv = style.makeCanvas()
            cv.SetRightMargin(0.2)
            cv.Draw("")
            cv.SetLogy()
            cv.SetLogx()
            cv.SetLogz()

            hist_mu.GetXaxis().SetTitle("m_{N} (GeV)")
            hist_mu.GetYaxis().SetTitle("|V_{lN}|^{2}")
            hist_mu.GetZaxis().SetTitle("#sigma/#sigma_{th}")
            hist_mu.Draw("COLZ")
            points_graph = ROOT.TGraph(npoints, array('d', mass_coupling_pair.T[0]), array('d',np.power(10, mass_coupling_pair.T[1])))
            points_graph.SetMarkerStyle(33)
            points_graph.SetMarkerSize(1)

            graph.SetLineColor(ROOT.kBlack)
            graph.SetMarkerColor(ROOT.kBlack)
            graph.SetFillColorAlpha(ROOT.kBlack, 0.3)
            graph.SetLineWidth(2)

            graph.Draw("SAMEC3")
            points_graph.Draw("P SAME")
            cv.SaveAs("limits/interpolation_{}_coupling_{}_{}.pdf".format(hnl_type, scenario, year))
            if hnl_type == 'majorana':
                graph_majorana = graph.Clone("majorana")
            elif hnl_type == 'dirac':
                graph_dirac = graph.Clone("dirac")
               
        cv = style.makeCanvas()
        cv.Draw("")
        cv.SetLogy()
        cv.SetLogx()
        cv.SetLogz()
        graph_majorana.Draw("AC3")
        graph_majorana.GetXaxis().SetTitle("m_{N} (GeV)")
        graph_majorana.GetYaxis().SetTitle("|V_{lN}|^{2}")

        graph_majorana.GetXaxis().SetLimits(1, 20.)
        graph_majorana.SetMinimum(1e-7)
        graph_majorana.SetMaximum(1e3)
        graph_majorana.GetXaxis().SetMoreLogLabels()

        graph_majorana.SetLineColor(ROOT.kAzure)
        graph_majorana.SetMarkerColor(ROOT.kAzure)
        graph_majorana.SetFillColorAlpha(ROOT.kAzure, 0.3)

        graph_dirac.SetLineColor(ROOT.kOrange)
        graph_dirac.SetMarkerColor(ROOT.kOrange)
        graph_dirac.SetFillColorAlpha(ROOT.kOrange, 0.3)

        graph_dirac.Draw("SAMEC3")
        graph_majorana.Draw("SAME C3")
        points_graph.Draw("SAME P")

        leg = style.makeLegend(0.16, 0.47, 0.5, 0.67)
        leg.AddEntry(graph_majorana, "Majorana", "lpf")
        leg.AddEntry(graph_dirac, "Dirac", "lpf")

        # Table 3: displaced HNL, Vmu , Dirac
        # Table 4: displaced HNL, Vmu, Majorana
        # Table 6: Prompt HNL, Vmu Dirac+Majorana
        # Table 5: prompt HNL, Ve Dirac+Majorana

        if scenario == 12:
            hist_displaced_dirac = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 3/Graph1D_y1")
            hist_displaced_majorana = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 4/Graph1D_y1")
            hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 6/Graph1D_y1")
            hist_displaced_dirac.SetLineColor(ROOT.kRed)
            hist_displaced_majorana.SetLineColor(ROOT.kBlue)
            hist_prompt.SetLineStyle(2)
            hist_displaced_dirac.SetLineWidth(2)
            hist_displaced_majorana.SetLineWidth(2)
            hist_prompt.SetLineWidth(2)

            hist_displaced_dirac.Draw("SAME")
            hist_displaced_majorana.Draw("SAME")
            hist_prompt.Draw("SAME")

            leg.AddEntry(hist_displaced_dirac, "ATLAS 3l, displaced LNC")
            leg.AddEntry(hist_displaced_majorana, "ATLAS 3l, displaced LNV")
            leg.AddEntry(hist_prompt, "ATLAS 3l, prompt")

        elif scenario == 2:
            hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 5/Graph1D_y1")
            hist_prompt.SetLineStyle(2)
            hist_prompt.SetLineWidth(2)
            hist_prompt.Draw("SAME")
            leg.AddEntry(hist_prompt, "ATLAS 3l, prompt")

        leg.Draw("SAME")
        style.makeText(0.18, 0.7, 0.2, 0.7, coupling_title)
        style.makeCMSText(0.18, 0.87, additionalText="Simulation Preliminary")
        style.makeLumiText(0.18, 0.8, year=year, lumi=lumi[year])
        cv.SaveAs("limits/limit_coupling_{}_{}.pdf".format(scenario, year))