import ROOT
import numpy as np
import style 

rootObj = []

def get_graph(filename, histname):
    root_file = ROOT.TFile(filename)
    graph = root_file.Get(histname)
    rootObj.append(graph)
    root_file.Close()
    return graph
    
    
color_ATLAS = style.newColorHLS(0.60,0.55,0.94)
color_LHCb = style.newColorHLS(0.10,0.51,0.98)#style.newColorHLS(0.80,0.53,0.98)
color_CMS = style.newColorHLS(0.02,0.52,0.97)


def getRefs(hnl_type,scenario):
    entries = []

    if scenario == 12:
        hist_displaced_dirac = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 3/Graph1D_y1")
        hist_displaced_majorana = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 4/Graph1D_y1")
        hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 6/Graph1D_y1")
        hist_displaced_dirac.SetLineColor(color_ATLAS.GetNumber())
        hist_displaced_majorana.SetLineColor(color_ATLAS.GetNumber())
        hist_prompt.SetLineColor(color_ATLAS.GetNumber())
        hist_prompt.SetLineStyle(3)
        hist_displaced_dirac.SetLineWidth(2)
        hist_displaced_majorana.SetLineWidth(2)
        hist_prompt.SetLineWidth(3)

        if hnl_type == "dirac":
            #hist_displaced_dirac.Draw("SAME")
            entries.append([
                [hist_displaced_dirac, "ATLAS, displaced, LNC"],
                [None, "JHEP 10 (2019) 265"]
            ])

        elif hnl_type == "majorana":
            hist_displaced_majorana.Draw("SAME")
            entries.append([
                [hist_displaced_majorana, "ATLAS, displaced, LNV"],
                [hist_prompt, "ATLAS, prompt, LNV"],
                [None, "JHEP 10 (2019) 265"]
            ])
            #hist_prompt.Draw("SAME")

        
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
            graph_lhcb = ROOT.TGraphErrors(6,mass_values_lhcb,coupling_values_lhcb, x_errors_lhcb,y_errors_lhcb)
            entries.append([
                [graph_lhcb, "LHCb, displaced, LNV"],
                [None, "EPJC 81 (2021) 248"]
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
            entries.append([
                [graph_lhcb, "LHCb, displaced, LNC"],
                [None, "EPJC 81 (2021) 248"]
            ])
        
        graph_lhcb.SetLineColor(color_LHCb.GetNumber())
        graph_lhcb.SetLineWidth(3)
        graph_lhcb.SetMarkerSize(0)
        graph_lhcb.SetLineStyle(6)
        #graph_lhcb.Draw("SAME L")
        
        if hnl_type == "majorana":
            mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 15., 15., 14., 11., 9., 7., 6.])
            coupling_values_exo_20_009 = np.asarray([1e-4, 2e-5, 5e-6, 2e-6, 1e-6, 5e-7, 4e-7, 3e-7, 3e-7, 3.5e-7, 6e-7, 1e-6, 2e-6, 6e-6, 4e-5, 1.5e-4, 1e-3, 2.5e-3])
        elif hnl_type == "dirac":
            mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 16.5, 16.5, 14., 11., 9., 7., 6.])
            coupling_values_exo_20_009 = np.asarray([1.7e-4, 2e-5, 7e-6, 2.7e-6, 1.5e-6, 7e-7, 5e-7, 4e-7, 3e-7, 3e-7, 4e-7, 1e-6, 2.5e-6, 1.5e-5, 1.5e-4, 4e-4, 2e-3, 5e-3])

        graph_exo_20_009 = ROOT.TGraph(len(mass_values_exo_20_009),mass_values_exo_20_009,coupling_values_exo_20_009)
        graph_exo_20_009.SetLineColor(color_CMS.GetNumber())
        graph_exo_20_009.SetLineWidth(3)
        graph_exo_20_009.SetMarkerSize(0)
        graph_exo_20_009.SetLineStyle(5)
        #graph_exo_20_009.Draw("SAME L")
        entries.append([
            [graph_exo_20_009, "CMS, displaced, 3l"],
            [None,"JHEP 07 (2022) 081"]
        ])
        
        

    if scenario == 2:
        if hnl_type == "majorana":
            hist_prompt = get_graph("hepdata/HEPData-ins1736526-v1.root", "Table 5/Graph1D_y1")
            hist_prompt.SetLineStyle(3)
            hist_prompt.SetLineWidth(3)
            hist_prompt.SetLineColor(color_ATLAS.GetNumber())
            #hist_prompt.Draw("SAME")
            
            entries.append([
                [hist_prompt, "ATLAS, prompt, LNV"],
                [None, "JHEP 10 (2019) 265"]
            ])

            mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 13., 13., 11., 9., 7., 6.])
            coupling_values_exo_20_009 = np.asarray([1e-3, 1e-4, 2.6e-5, 1.2e-5, 6.3e-6, 1.6e-6, 1.2e-6, 9.3e-7, 8.5e-7, 8.8e-7, 1e-6, 6.5e-6, 3.1e-5, 1.4e-4, 8e-4, 2.3e-3])

        elif hnl_type == "dirac":
            mass_values_exo_20_009 = np.asarray([1., 2., 3., 4., 5., 7., 8., 9., 11., 12., 14., 14., 11., 9., 7., 6.])
            coupling_values_exo_20_009 = np.asarray([1e-3, 1e-4, 3.2e-5, 1.5e-5, 7.4e-6, 2e-6, 1.5e-6, 1e-6, 9e-7, 8.7e-7, 1e-6, 1e-5, 8.7e-5, 3.5e-4, 2e-3, 5e-3])

        

        graph_exo_20_009 = ROOT.TGraph(len(mass_values_exo_20_009),mass_values_exo_20_009,coupling_values_exo_20_009)
        graph_exo_20_009.SetLineColor(color_CMS.GetNumber())
        graph_exo_20_009.SetLineWidth(3)
        graph_exo_20_009.SetMarkerSize(0)
        graph_exo_20_009.SetLineStyle(5)
        #graph_exo_20_009.Draw("SAME L")
        entries.append([
            [graph_exo_20_009, "CMS, displaced, 3l"],
            [None, "JHEP 07 (2022) 081"]
        ])


    return entries


         
