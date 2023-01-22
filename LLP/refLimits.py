# ATLAS Results:
# Table 3: displaced HNL, Vmu , Dirac
# Table 4: displaced HNL, Vmu, Majorana
# Table 6: Prompt HNL, Vmu Dirac+Majorana
# Table 5: prompt HNL, Ve Dirac+Majorana

def get_graph(filename, histname):
    root_file = ROOT.TFile(filename)
    graph = root_file.Get(histname)
    root_file.Close()
    return graph


ref_ATLAS_mumu_dirac_displaced = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 3/Graph1D_y1")
ref_ATLAS_mumu_majorana_displaced = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 4/Graph1D_y1")
ref_ATLAS_mumu_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 6/Graph1D_y1")
ref_ATLAS_ee_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 5/Graph1D_y1")


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
                graph_lhcb = ROOT.TGraphErrors(6,mass_values_lhcb,coupling_values_lhcb, x_errors_lhcb,y_errors_lhcb)
                graph_lhcb.SetLineColor(ROOT.kRed)
                graph_lhcb.SetLineWidth(3)
                graph_lhcb.SetMarkerSize(0)
                graph_lhcb.SetLineStyle(6)
                graph_lhcb.Draw("SAME L")
                leg.AddEntry(graph_lhcb, "LHCb")

                if hnl_type == "majorana":
                    mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 15., 15., 14., 11., 9., 7., 6.])
                    coupling_values_exo_20_009 = np.asarray([1e-4, 2e-5, 5e-6, 2e-6, 1e-6, 5e-7, 4e-7, 3e-7, 3e-7, 3.5e-7, 6e-7, 1e-6, 2e-6, 6e-6, 4e-5, 1.5e-4, 1e-3, 2.5e-3])
                elif hnl_type == "dirac":
                    mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 16.5, 16.5, 14., 11., 9., 7., 6.])
                    coupling_values_exo_20_009 = np.asarray([1.7e-4, 2e-5, 7e-6, 2.7e-6, 1.5e-6, 7e-7, 5e-7, 4e-7, 3e-7, 3e-7, 4e-7, 1e-6, 2.5e-6, 1.5e-5, 1.5e-4, 4e-4, 2e-3, 5e-3])

                graph_exo_20_009 = ROOT.TGraph(len(mass_values_exo_20_009),mass_values_exo_20_009,coupling_values_exo_20_009)
                graph_exo_20_009.SetLineColor(ROOT.kViolet)
                graph_exo_20_009.SetLineWidth(3)
                graph_exo_20_009.SetMarkerSize(0)
                graph_exo_20_009.SetLineStyle(5)
                graph_exo_20_009.Draw("SAME L")
                leg.AddEntry(graph_exo_20_009, "EXO-20-009")
                if hnl_type == "dirac":
                    leg.AddEntry(median, "", "")

            if scenario == 2:
                if hnl_type == "majorana":
                    hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 5/Graph1D_y1")
                    hist_prompt.SetLineStyle(3)
                    hist_prompt.SetLineWidth(3)
                    hist_prompt.Draw("SAME")

                    mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 13., 13., 11., 9., 7., 6.])
                    coupling_values_exo_20_009 = np.asarray([1e-3, 1e-4, 2.6e-5, 1.2e-5, 6.3e-6, 1.6e-6, 1.2e-6, 9.3e-7, 8.5e-7, 8.8e-7, 1e-6, 6.5e-6, 3.1e-5, 1.4e-4, 8e-4, 2.3e-3])

                elif hnl_type == "dirac":
                    mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 14., 11., 9., 7., 6.])
                    coupling_values_exo_20_009 = np.asarray([1e-3, 1e-4, 3.2e-5, 1.5e-5, 7.4e-6, 2e-6, 1.5e-6, 1e-6, 9e-7, 8.7e-7, 1e-6, 1e-5, 8.7e-5, 3.5e-4, 2e-3, 5e-3])

                leg.AddEntry(hist_prompt, "ATLAS 3l, prompt")

                graph_exo_20_009 = ROOT.TGraph(len(mass_values_exo_20_009),mass_values_exo_20_009,coupling_values_exo_20_009)
                graph_exo_20_009.SetLineColor(ROOT.kViolet)
                graph_exo_20_009.SetLineWidth(3)
                graph_exo_20_009.SetMarkerSize(0)
                graph_exo_20_009.SetLineStyle(5)
                graph_exo_20_009.Draw("SAME L")

