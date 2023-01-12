import matplotlib
matplotlib.use('Agg') 
import style
import ROOT
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetPalette(ROOT.kDarkRainBow)
import json
from array import array
import scipy
import scipy.spatial
from scipy import interpolate
import numpy as np
import math
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

json_path = "/vols/cms/hsfar/jsons"

K_FACTOR = 1.1

blinded = False

limit_var_list = ["exp0", "exp+1", "exp+2", "exp-1", "exp-2"]

if not blinded:
    limit_var_list += ["obs"]

lumi = {"2016": 36, "2017": 42, "2018": 60, "combined": 138}
#years = ["2016", "2017", "2018", "combined"]

years = ["2016"]

coupling_dict = {}
#coupling_dict[1.0] = ["emutau", "V_{Ne} : V_{N#mu} : V_{N#tau} = 1 : 1 : 1"]
#coupling_dict[2.0] = ["ee", "V_{Ne} : V_{N#mu} : V_{N#tau} = 1 : 0 : 0"]
#coupling_dict[7.0] = ["emu", "V_{Ne} : V_{N#mu} : V_{N#tau} = 1 : 1 : 0"]
coupling_dict[12.0] = ["mumu", "V_{Ne} : V_{N#mu} : V_{N#tau} = 0 : 1 : 0"]
#coupling_dict[47.0] = ["etau", "V_{Ne} : V_{N#mu} : V_{N#tau} = 1 : 0 : 1"]
#coupling_dict[52.0] = ["mutau", "V_{Ne} : V_{N#mu} : V_{N#tau} = 0 : 1 : 1"]

n_bins = 200

mass_range = np.geomspace(1., 20., num=n_bins)
log_coupling_range = np.linspace(-7, 0., num=n_bins)
coupling_range = np.power(10, log_coupling_range)


def smoothPoints(points,values,addTrigs=False,splineFit=False):
    delaunay = scipy.spatial.Delaunay(points)
    
    newPoints = []
    newValues = []
    for i,trigIdx in enumerate(delaunay.simplices):
        smooth1 = 0.5*(values[trigIdx[0]]+values[trigIdx[1]])
        smooth2 = 0.5*(values[trigIdx[0]]+values[trigIdx[2]])
        smooth3 = 0.5*(values[trigIdx[1]]+values[trigIdx[2]])
        smooth = 0.3333*(values[trigIdx[0]]+values[trigIdx[1]]+values[trigIdx[2]])
        
        midpoint = 0.333*(points[trigIdx[0]]+points[trigIdx[1]]+points[trigIdx[2]])
        
        newPoints.append(midpoint)
        newValues.append(smooth)
        
    newPoints = np.stack(newPoints,axis=0)
    newValues = np.array(newValues)
    
    smoothValues = np.zeros_like(values)
    
    newDelaunay = scipy.spatial.Delaunay(newPoints)
    for ipoint in range(points.shape[0]):
        simplexIndex = newDelaunay.find_simplex(points[ipoint])
        pointIndices = newDelaunay.simplices[simplexIndex]
        
        smoothValues[ipoint] = 0.5*values[ipoint]+0.5*(0.333*(newValues[pointIndices[0]]+newValues[pointIndices[1]]+newValues[pointIndices[2]]))
        
    allPoints = np.concatenate([newPoints,points],axis=0)
    allValues = np.concatenate([newValues,smoothValues],axis=0)
    
    
    if splineFit:
        tck = interpolate.bisplrep(allPoints[:,0], allPoints[:,1], allValues, kx=2, ky=2, s=0.01*allPoints.shape[0], eps=1e-8)
        for ipoint in range(newPoints.shape[0]):
            newValues[ipoint] = interpolate.bisplev(newPoints[ipoint,0],newPoints[ipoint,1], tck)
        for ipoint in range(points.shape[0]):
            smoothValues[ipoint] = interpolate.bisplev(points[ipoint,0],points[ipoint,1], tck)
            
        
    if addTrigs:
        return np.concatenate([newPoints,points],axis=0),np.concatenate([newValues,smoothValues],axis=0)
    else:
        return points,smoothValues
        

    

def interpolatedFct(points,values,transformations=[]):
    transformedPoints = []
    if len(transformations)!=points.shape[1]:
        raise Exception("ERROR: transformation (%i) need to be of same length as point dims (%i)"%(
            len(transformations),
            points.shape[1]
        ))
    
    for i,transformation in enumerate(transformations):
        transformedPoints.append(
            transformation(points[:,i])
        )
    transformedPoints = np.stack(transformedPoints,axis=1)
    
    transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=True,splineFit=False)
    transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=True,splineFit=False)
    transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=True,splineFit=True)
    #transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=False,splineFit=True)
    print (points.shape,transformedPoints.shape)
    #transformedPoints,values = smoothPoints(transformedPoints,values)
    #transformedPoints,values = smoothPoints(transformedPoints,values)
    
    delaunay = scipy.spatial.Delaunay(transformedPoints)
    
    
    
    
    
    def getValue(point):
        transformedPoint = []
        for i,transformation in enumerate(transformations):
            transformedPoint.append(
                transformation(point[i])
            )
        transformedPoint = np.array(transformedPoint)
        simplexIndex = delaunay.find_simplex(transformedPoint)
        pointIndices = delaunay.simplices[simplexIndex]
        simplexPoints = transformedPoints[pointIndices]
        
        minSimplexPoints = np.amin(simplexPoints,axis=0)
        maxSimplexPoints = np.amax(simplexPoints,axis=0)
        distance = np.fabs(maxSimplexPoints-minSimplexPoints)+1e-6
        
        
        weight = np.prod(1.0-np.abs(simplexPoints-transformedPoint)/distance,axis=1)
        return 0.9*np.sum(values[pointIndices]*weight)/np.sum(weight)+0.1*np.mean(values[pointIndices])
     
    return getValue
'''
def interpolatedFct(points,values,transformations=[]):
    transformedPoints = []
    if len(transformations)!=points.shape[1]:
        raise Exception("ERROR: transformation (%i) need to be of same length as point dims (%i)"%(
            len(transformations),
            points.shape[1]
        ))
    
    for i,transformation in enumerate(transformations):
        transformedPoints.append(
            transformation(points[:,i])
        )
    transformedPoints = np.stack(transformedPoints,axis=1)
    f = interpolate.interp2d(transformedPoints[:,0], transformedPoints[:,1], values, kind='cubic')
    def getValue(point):
        transformedPoint = []
        for i,transformation in enumerate(transformations):
            transformedPoint.append(
                transformation(point[i])
            )
        transformedPoint = np.array(transformedPoint)
        return f(transformedPoint[0],transformedPoint[1])
    return getValue
'''

with open("/vols/cms/LLP/gridpackLookupTable.json") as lookup_table_file:
    lookup_table = json.load(lookup_table_file)

for year in years:
    print(year)
    for scenario in coupling_dict.keys():
        print("Analyzing coupling scenario: "+str(scenario))
        coupling_text = coupling_dict[scenario][0]
        coupling_title = coupling_dict[scenario][1]

        for hnl_type in ["dirac"]:#,"majorana"]:
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
                
                

                #print mass, coupling, xsec

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
            
            
            results_theory = np.zeros((n_bins,n_bins),dtype=np.float32)
            results = np.zeros((n_bins,n_bins),dtype=np.float32)
            results_plus = np.zeros((n_bins,n_bins),dtype=np.float32)
            results_minus = np.zeros((n_bins,n_bins),dtype=np.float32)
            results_plus_two = np.zeros((n_bins,n_bins),dtype=np.float32)
            results_minus_two = np.zeros((n_bins,n_bins),dtype=np.float32)
            
            results_theory_fct = interpolatedFct(mass_coupling_pair,log_theory_points,transformations=[lambda x: np.log(x), lambda x: x])
            results_fct = interpolatedFct(mass_coupling_pair,log_expected_points,transformations=[lambda x: np.log(x), lambda x: x])
            results_plus_fct = interpolatedFct(mass_coupling_pair,log_expected_points_plus,transformations=[lambda x: np.log(x), lambda x: x])
            results_minus_fct = interpolatedFct(mass_coupling_pair,log_expected_points_minus,transformations=[lambda x: np.log(x), lambda x: x])
            results_plus_two_fct = interpolatedFct(mass_coupling_pair,log_expected_points_plus_two,transformations=[lambda x: np.log(x), lambda x: x])
            results_minus_two_fct = interpolatedFct(mass_coupling_pair,log_expected_points_minus_two,transformations=[lambda x: np.log(x), lambda x: x])
            
            for i in range(n_bins):
                for j in range(n_bins):
                    results_theory[i,j] = results_theory_fct([mass_range[i],log_coupling_range[j]])
                    results[i,j] = results_fct([mass_range[i],log_coupling_range[j]])
                    results_plus[i,j] = results_plus_fct([mass_range[i],log_coupling_range[j]])
                    results_minus[i,j] = results_minus_fct([mass_range[i],log_coupling_range[j]])
                    results_plus_two[i,j] = results_plus_two_fct([mass_range[i],log_coupling_range[j]])
                    results_minus_two[i,j] = results_minus_two_fct([mass_range[i],log_coupling_range[j]])
            
            '''
            grid_x, grid_y = np.meshgrid(mass_range, log_coupling_range, indexing='ij')
            fit_method = 'cubic'
            
            results_theory = interpolate.griddata(mass_coupling_pair, log_theory_points, (grid_x, grid_y), method=fit_method,rescale=True)
            results = interpolate.griddata(mass_coupling_pair, log_expected_points, (grid_x, grid_y), method=fit_method,rescale=True)
            results_plus = interpolate.griddata(mass_coupling_pair, log_expected_points_plus, (grid_x, grid_y), method=fit_method,rescale=True)
            results_minus = interpolate.griddata(mass_coupling_pair, log_expected_points_minus, (grid_x, grid_y), method=fit_method,rescale=True)
            results_plus_two = interpolate.griddata(mass_coupling_pair, log_expected_points_plus_two, (grid_x, grid_y), method=fit_method,rescale=True)
            results_minus_two = interpolate.griddata(mass_coupling_pair, log_expected_points_minus_two, (grid_x, grid_y), method=fit_method,rescale=True)
            
            for i in range(n_bins):
                for j in range(n_bins):
                    if not np.isfinite(results[i,j]).all():
                    
                        print i,j,results[i,j]
            '''
            
            if not blinded:
                
                results_obs = np.zeros((n_bins,n_bins),dtype=np.float32)
                results_obs_fct = interpolatedFct(mass_coupling_pair,log_observed_points,transformations=[lambda x: np.log(x), lambda x: x])
                for i in range(n_bins):
                    for j in range(n_bins):
                        results_obs[i,j] = results_obs_fct([mass_range[i],log_coupling_range[j]])
                
                #results_obs = interpolate.griddata(mass_coupling_pair, log_observed_points, (grid_x, grid_y), method=fit_method,rescale=True)


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

            hist_mu.GetXaxis().SetTitle("m#lower[0.3]{#scale[0.7]{N}} (GeV)")
            hist_mu.GetYaxis().SetTitle("|V#lower[0.3]{#scale[0.7]{#font[12]{l}N}}|#lower[-0.7]{#scale[0.7]{2}}")
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

            leg = style.makeLegend(0.7, 0.7, 0.89, 0.85)
            leg.AddEntry(median, "expected",'L')

            if not blinded:
                leg.AddEntry(observed, "observed",'L')

            leg.Draw()

            text = style.makeText(0.18, 0.8, 0.2, 0.80, coupling_title+", "+hnl_type.capitalize()+ " HNL")
            text.SetTextFont(63)
            text.SetTextSize(31)
            #Text(0.18, 0.89, additionalText="Simulation Preliminary")
            #style.makeLumiText(0.64, 0.95, year=year, lumi=lumi[year])
            style.makeLumiText(0.64, 0.95, year=year, lumi=lumi[year])

            cv.SaveAs("limits/interpolation_{}_coupling_{}_{}.pdf".format(hnl_type, scenario, year))
            cv.SaveAs("limits/interpolation_{}_coupling_{}_{}.png".format(hnl_type, scenario, year))

            hist_mu.Draw("AXIS")
            #leg = style.makeLegend(0.16, 0.4, 0.5, 0.77)
            leg = style.makeLegend(0.58, 0.49, 0.88, 0.88)
            leg.SetTextSize(24)

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
                
            rootObj = []
            for ipoint in range(mass_coupling_pair.shape[0]):
                m = ROOT.TMarker(mass_coupling_pair[ipoint,0],10**mass_coupling_pair[ipoint,1],20)
                rootObj.append(m)
                m.SetMarkerSize(0.8)
                m.SetMarkerColor(ROOT.kBlack)
                m.Draw()
            
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
                leg.AddEntry(graph_exo_20_009, "EXO-20-009")

                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")
                if hnl_type == "dirac":
                    leg.AddEntry(median, "", "")

            else:
                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")
                leg.AddEntry(median, "", "")           
                leg.AddEntry(median, "", "")           

            #text = style.makeText(0.18, 0.8, 0.2, 0.80, coupling_title+", "+hnl_type.capitalize())
            text = style.makeText(0.18, 0.76, 0.2, 0.80, coupling_title)
            text2 = style.makeText(0.18, 0.7, 0.2, 0.70, hnl_type.capitalize())
            text.SetTextFont(63)
            text.SetTextSize(31)

            text2.SetTextFont(63)
            text2.SetTextSize(31)
            style.makeCMSText(0.18, 0.89)#, additionalText="Simulation Preliminary")
            #style.makeLumiText(0.9, 0.95, year=year, lumi=lumi[year])
            style.makeLumiText(0.9, 0.95, year="13 TeV", lumi=lumi[year])

            cv.SaveAs("limits/{}_coupling_{}_{}.pdf".format(hnl_type, scenario, year))
            cv.SaveAs("limits/{}_coupling_{}_{}.png".format(hnl_type, scenario, year))
