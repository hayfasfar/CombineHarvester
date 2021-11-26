import matplotlib
matplotlib.use('Agg') 
import style
import ROOT
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetPalette(ROOT.kDarkRainBow)
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
    lu_infos = lookup_table[proc.replace('pt20', 'all')]['weights'][str(int(scenario))]
    xsec = lu_infos['xsec']['nominal']*K_FACTOR
    coupling = lu_infos['couplings']['Ve']+lu_infos['couplings']['Vmu']+lu_infos['couplings']['Vtau']
    if coupling not in (2, 12, 67):
        coupling = coupling/2
    coupling = coupling ** 2
    mass = lookup_table[proc.replace('pt20', 'all')]['mass']

    return mass, coupling, xsec

def interpolate_point(coupling_points, mu_points):
    for i, (coupling, mu) in enumerate(zip(coupling_points, mu_points)):
        if i < len(mu_points)-1 and mu > 0 and mu_points[i+1] < 0:
            return (coupling+coupling_points[i+1])/2
    return min(coupling_points)

def parse_limit_json(f, scenario=12):
    # limit json aggreggate 
    with open(os.path.join(json_path, f)) as json_file:
        xsec_dict = json.load(json_file)
    if str(scenario) not in xsec_dict.keys():
        return None
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

json_path = "jsons"

K_FACTOR = 1.1

blinded = True

limit_var_list = ["exp0", "exp+1", "exp+2", "exp-1", "exp-2"]

if not blinded:
    limit_var_list += ["obs"]

lumi = {"2016": 35.9, "2017": 41.5, "2018": 59.7, "combined": 137.1}
years = ["2016", "2017", "2018", "combined"]

coupling_dict = {}
coupling_dict[1.0] = ["emutau", "V_{e} : V_{#mu} : V_{#tau} = 1 : 1 : 1"]
coupling_dict[2.0] = ["ee", "V_{e} : V_{#mu} : V_{#tau} = 1 : 0 : 0"]
coupling_dict[7.0] = ["emu", "V_{e} : V_{#mu} : V_{#tau} = 1 : 1 : 0"]
coupling_dict[12.0] = ["mumu", "V_{e} : V_{#mu} : V_{#tau} = 0 : 1 : 0"]
coupling_dict[47.0] = ["etau", "V_{e} : V_{#mu} : V_{#tau} = 1 : 0 : 1"]
coupling_dict[52.0] = ["mutau", "V_{e} : V_{#mu} : V_{#tau} = 0 : 1 : 1"]

n_bins = 200

mass_range = np.geomspace(1., 19., num=n_bins)
log_coupling_range = np.linspace(-7, 0., num=n_bins)
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
            for exp_var in limit_var_list+["theory"]:
                sigma_dict[exp_var] = []

            for f in os.listdir(json_path):
                if year not in f or ".json" not in f or hnl_type not in f:
                    continue
            
                xsec_dict = parse_limit_json(f, scenario)
                if xsec_dict is None:
                    continue
                mass, coupling, xsec = parse_lookup_table(f, lookup_table)

                # if "ctau1p0e-05" in f: 
                #      continue
                # if "ctau1p0e-03_massHNL10p0" in f:
                #     continue
                # if "ctau1p0e-03_massHNL12p0" in f:
                #     continue
                # if "ctau1p0e-02_massHNL16p0" in f:
                #     continue
                if "_pt20_" not in f and f.replace("_all_", "_pt20_") in os.listdir(json_path):
                    print("Skipping "+f+", higher stat. sample exists")
                    continue

                sigma_dict["theory"].append(xsec)

                

                for exp_var in limit_var_list:
                    sigma_dict[exp_var].append(xsec_dict[exp_var])

                masses.append(mass)
                couplings.append(coupling)

            df = pd.DataFrame.from_dict(sigma_dict)
            df['mass'] = masses
            df['coupling'] = couplings
    
            npoints = len(df)
            mass_coupling_pair = np.array([df['mass'], np.log10(df['coupling'])]).T
            log_theory_points = np.log10(np.array(df['theory']))
            log_expected_points = np.log10(np.array(df['exp0']))
            log_expected_points_plus = np.log10(np.array(df['exp+1']))
            log_expected_points_minus = np.log10(np.array(df['exp-1']))
            log_expected_points_plus_two = np.log10(np.array(df['exp+2']))
            log_expected_points_minus_two = np.log10(np.array(df['exp-2']))

            if not blinded:
                log_observed_points = np.log10(np.array(df['obs']))

            df.to_csv("csv/limit_{}_coupling_{}_{}.csv".format(hnl_type, scenario, year))
            n_bins = 200
            grid_x, grid_y = np.meshgrid(mass_range, log_coupling_range, indexing='ij')
            fit_method = 'cubic'

            results_theory = interpolate.griddata(mass_coupling_pair, log_theory_points, (grid_x, grid_y), method=fit_method)
            results = interpolate.griddata(mass_coupling_pair, log_expected_points, (grid_x, grid_y), method=fit_method)
            results_plus = interpolate.griddata(mass_coupling_pair, log_expected_points_plus, (grid_x, grid_y), method=fit_method)
            results_minus = interpolate.griddata(mass_coupling_pair, log_expected_points_minus, (grid_x, grid_y), method=fit_method)
            results_plus_two = interpolate.griddata(mass_coupling_pair, log_expected_points_plus_two, (grid_x, grid_y), method=fit_method)
            results_minus_two = interpolate.griddata(mass_coupling_pair, log_expected_points_minus_two, (grid_x, grid_y), method=fit_method)

            if not blinded:
                results_obs = interpolate.griddata(mass_coupling_pair, log_observed_points, (grid_x, grid_y), method=fit_method)


            hist_mu = ROOT.TH2D("mu"+hnl_type+str(scenario), "mu", n_bins-1, mass_range, n_bins-1, coupling_range)
            crossing_points = []
            masses_at_crossing_points = []
            errors_up = []
            errors_down = []
            errors_up_two = []
            errors_down_two = []

            yellow = ROOT.TGraph(2*n_bins)    # yellow band
            green = ROOT.TGraph(2*n_bins)     # green band
            median = ROOT.TGraph(n_bins)      # median line 
            if not blinded:
                observed = ROOT.TGraph(n_bins)    

            for i in range(n_bins):
                mass = mass_range[i]
                coupling_values = []
                mu_values = []
                
                if not blinded:
                    mu_obs_values = []

                mu_plus_values = []
                mu_minus_values = []
                mu_plus_two_values = []
                mu_minus_two_values = []

                for j in range(n_bins):
                    coupling_values.append(log_coupling_range[j])
                    exp = results[i, j]

                    if not blinded:
                        obs = results_obs[i, j]

                    exp_plus = results_plus[i, j]
                    exp_minus = results_minus[i, j]
                    exp_plus_two = results_plus_two[i, j]
                    exp_minus_two = results_minus_two[i, j]
                    theory = results_theory[i, j]

                    mu_values.append(get_mu(theory, exp))

                    if not blinded:
                        mu_obs_values.append(get_mu(theory, obs))

                    mu_plus_values.append(get_mu(theory, exp_plus))
                    mu_minus_values.append(get_mu(theory, exp_minus))
                    mu_plus_two_values.append(get_mu(theory, exp_plus_two))
                    mu_minus_two_values.append(get_mu(theory, exp_minus_two))

                    hist_mu.SetBinContent(i+1, j+1, np.power(10, get_mu(theory, exp)))

    
                crossing_point = np.power(10, interpolate_point(coupling_values, mu_values))

                if not blinded:
                    crossing_point_observed = np.power(10, interpolate_point(coupling_values, mu_obs_values))

                crossing_point_plus = np.power(10, interpolate_point(coupling_values, mu_plus_values))
                crossing_point_minus = np.power(10, interpolate_point(coupling_values, mu_minus_values))
                crossing_point_plus_two = np.power(10, interpolate_point(coupling_values, mu_plus_two_values)) 
                crossing_point_minus_two = np.power(10, interpolate_point(coupling_values, mu_minus_two_values))
                yellow.SetPoint(    i,    mass, crossing_point_plus_two ) # + 2 sigma
                green.SetPoint(     i,    mass, crossing_point_plus ) # + 1 sigma
                median.SetPoint(    i,    mass, crossing_point ) # median

                if not blinded:
                    observed.SetPoint(    i,    mass, crossing_point_observed ) # observed

                green.SetPoint(  2*n_bins-1-i, mass, crossing_point_minus ) # - 1 sigma
                yellow.SetPoint( 2*n_bins-1-i, mass, crossing_point_minus_two ) # - 2 sigma

            points_graph = ROOT.TGraph(npoints, array('d', mass_coupling_pair.T[0]), array('d',np.power(10, mass_coupling_pair.T[1])))
            points_graph.SetMarkerStyle(33)
            points_graph.SetMarkerSize(1)
            points_graph.SetMarkerColor(0)

            cv = style.makeCanvas()
            cv.SetRightMargin(0.1) # 0.1
            cv.SetLogy()
            cv.SetLogx()
            cv.SetLogz()
            cv.Draw("")

            hist_mu.GetXaxis().SetTitle("m_{N} [GeV]")
            hist_mu.GetYaxis().SetTitle("|V_{lN}|^{2}")
            #hist_mu.GetZaxis().SetTitle("#sigma/#sigma_{th}")
            hist_mu.Draw("COLZ")
            hist_mu.SetMaximum(1e3)

            points_graph.Draw("P SAME")

            median.SetLineWidth(2)

            if not blinded:
                observed.SetLineColor(ROOT.kRed)
                observed.SetLineWidth(2)

            median.Draw('Lsame')

            if not blinded:
                observed.Draw('Lsame')

            median.SetTitle("")

            leg = style.makeLegend(0.16, 0.5, 0.5, 0.77)
            leg.AddEntry(median, "expected",'L')

            if not blinded:
                leg.AddEntry(observed, "observed",'L')

            leg.Draw()

            text = style.makeText(0.18, 0.8, 0.2, 0.80, coupling_title+", "+hnl_type.capitalize()+ " HNL")
            text.SetTextFont(63)
            text.SetTextSize(31)
            #Text(0.18, 0.89, additionalText="Simulation Preliminary")
            style.makeLumiText(0.64, 0.95, year=year, lumi=lumi[year])

            cv.SaveAs("limits/interpolation_{}_coupling_{}_{}.pdf".format(hnl_type, scenario, year))
            cv.SaveAs("limits/interpolation_{}_coupling_{}_{}.png".format(hnl_type, scenario, year))

            hist_mu.Draw("AXIS")
            leg = style.makeLegend(0.16, 0.5, 0.5, 0.77)

            if not blinded:
                leg.AddEntry(observed, "observed",'L')
            leg.AddEntry(median, "expected",'L')

            yellow.SetFillColor(ROOT.kOrange)
            yellow.SetLineColor(ROOT.kOrange)
            yellow.SetFillStyle(1001)
        
            green.SetFillColor(ROOT.kGreen+1)
            green.SetLineColor(ROOT.kGreen+1)
            green.SetFillStyle(1001)

            yellow.Draw('Fsame')
            green.Draw('Fsame')
            median.Draw("Lsame")

            if not blinded:
                observed.Draw("Lsame")

            leg.AddEntry(green, "#pm 1 std. deviation",'f')
            leg.AddEntry(yellow,"#pm 2 std. deviation",'f')

            leg.Draw("same")

            # ATLAS Results:
            # Table 3: displaced HNL, Vmu , Dirac
            # Table 4: displaced HNL, Vmu, Majorana
            # Table 6: Prompt HNL, Vmu Dirac+Majorana
            # Table 5: prompt HNL, Ve Dirac+Majorana

            if scenario == 12:
                hist_displaced_dirac = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 3/Graph1D_y1")
                hist_displaced_majorana = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 4/Graph1D_y1")
                hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 6/Graph1D_y1")
                hist_displaced_dirac.SetLineColor(ROOT.kBlue)
                hist_displaced_majorana.SetLineColor(ROOT.kBlue)
                hist_prompt.SetLineColor(ROOT.kBlue)
                hist_prompt.SetLineStyle(3)
                hist_displaced_dirac.SetLineWidth(2)
                hist_displaced_majorana.SetLineWidth(2)
                hist_prompt.SetLineWidth(3)

                if hnl_type == "dirac":
                    hist_displaced_dirac.Draw("SAME")
                    leg.AddEntry(hist_displaced_dirac, "ATLAS displaced")

                elif hnl_type == "majorana":
                    hist_displaced_majorana.Draw("SAME")
                    leg.AddEntry(hist_displaced_majorana, "ATLAS displaced")
                    hist_prompt.Draw("SAME")
                    leg.AddEntry(hist_prompt, "ATLAS prompt")

                mass_values_lhcb = np.asarray([5., 10., 15., 20., 30., 50.])
                x_errors_lhcb = np.zeros(6)
                if hnl_type == "majorana":
                    coupling_values_lhcb = np.asarray([
                        0.0004912281,
                        0.0002422232,
                        0.000210303,
                        0.0002470208,
                        0.0003078862,
                        0.0008244407]
                    )
                    y_errors_lhcb = np.asarray([
                        1.534516e-05,
                        2.356009e-06,
                        2.852664e-06,
                        2.554021e-06,
                        2.57904e-06,
                        5.739265e-06   
                    ])
                else:
                    coupling_values_lhcb = np.asarray([
                        0.001240287,
                            0.0007974421,
                            0.001519619,
                            0.001442267,
                            0.002283984,
                            0.009250932
                    ])
                    y_errors_lhcb = np.asarray([
                        1.306788e-05,
                        1.036377e-05,
                        1.792985e-05,
                        5.016125e-05,
                        0.0003602391,
                        0.0002851714
                    ])
                print(mass_values_lhcb, coupling_values_lhcb)
                graph_lhcb = ROOT.TGraphErrors(6,mass_values_lhcb,coupling_values_lhcb, x_errors_lhcb,y_errors_lhcb)
                graph_lhcb.SetLineColor(ROOT.kRed)
                graph_lhcb.SetLineWidth(3)
                graph_lhcb.SetMarkerSize(0)
                graph_lhcb.SetLineStyle(6)
                graph_lhcb.Draw("SAME L")
                leg.AddEntry(graph_lhcb, "LHCb")

            elif scenario == 2 and hnl_type == "majorana":
                hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 5/Graph1D_y1")
                hist_prompt.SetLineStyle(3)
                hist_prompt.SetLineWidth(3)
                hist_prompt.Draw("SAME")
                leg.AddEntry(hist_prompt, "ATLAS 3l, prompt")
                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")
            else:
                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")           

            text = style.makeText(0.18, 0.8, 0.2, 0.80, coupling_title+", "+hnl_type.capitalize())
            text.SetTextFont(63)
            text.SetTextSize(31)
            style.makeCMSText(0.18, 0.89, additionalText="Simulation Preliminary")
            style.makeLumiText(0.64, 0.95, year=year, lumi=lumi[year])

            cv.SaveAs("limits/{}_coupling_{}_{}.pdf".format(hnl_type, scenario, year))
            cv.SaveAs("limits/{}_coupling_{}_{}.png".format(hnl_type, scenario, year))