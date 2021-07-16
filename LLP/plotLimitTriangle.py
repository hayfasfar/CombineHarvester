import numpy
import json
import scipy
import scipy.spatial
import ROOT
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import ctypes
import itertools 

typeHNL = "Dirac"
ctau = "1"

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptDate(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFile(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.gStyle.SetLabelFont(43,"XYZ")
ROOT.gStyle.SetLabelSize(28,"XYZ")
ROOT.gStyle.SetTitleFont(43,"XYZ")
ROOT.gStyle.SetTitleSize(31,"XYZ")

colors = []
    
def newColorRGB(red,green,blue):
    newColorRGB.colorindex+=1
    color=ROOT.TColor(newColorRGB.colorindex,red,green,blue)
    colors.append(color)
    return color
    
newColorRGB.colorindex=100

    
def HLS2RGB(hue,light,sat):
    r, g, b = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
    ROOT.TColor.HLS2RGB(
        int(round(hue*255.)),
        int(round(light*255.)),
        int(round(sat*255.)),
        r,g,b
    )
    return r.value/255.,g.value/255.,b.value/255.
    
def newColorHLS(hue,light,sat):
    r,g,b = HLS2RGB(hue,light,sat)
    return newColorRGB(r,g,b)
    


colorList = [
    [0.,newColorHLS(0.8, 0.6,0.95)],
    [0.,newColorHLS(0.7, 0.61,0.95)],
    [0.,newColorHLS(0.6, 0.63,0.95)],
    [0.,newColorHLS(0.4, 0.65,0.9)],
    [0.,newColorHLS(0.15, 0.68,0.9)],
    [0.,newColorHLS(0.0, 0.72,0.9)],
]

'''
lumiMin = min(map(lambda x:x[1].GetLight(),colorList))
lumiMax = max(map(lambda x:x[1].GetLight(),colorList))

for color in colorList:
    color[0] = ((color[1].GetLight()-lumiMin)/(lumiMax-lumiMin))
'''
#stops = numpy.array(map(lambda x:x[0],colorList))
stops = numpy.linspace(0,1,len(colorList))
red   = numpy.array(map(lambda x:x[1].GetRed(),colorList))
green = numpy.array(map(lambda x:x[1].GetGreen(),colorList))
blue  = numpy.array(map(lambda x:x[1].GetBlue(),colorList))


start=ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 200)
print start
ROOT.gStyle.SetNumberContours(200)



def generate(n=11):
    index = 0
    for i,l3 in enumerate(numpy.linspace(0,1,n)):
        for l2 in numpy.linspace(0,max(0,1-l3),n-i):
            l1 = max(0,1 - l2 - l3)
            yield index,l1,l2,l3
            index+=1

def makePlot(res, ctau="1", typeHNL="Dirac"):
    def applyAxisStyle(axis):
        axis.SetLabelFont(43)
        axis.SetLabelSize(29)
        axis.SetTitleFont(43)
        axis.SetTitleSize(33)
        axis.CenterTitle(True)
        
    def getX(ve,vmu,vtau):
        x = 0*ve+1*vmu+0.5*vtau
        x*=0.75
        x+=0.1
        return x
        
    def getY(ve,vmu,vtau):
        y = 0*ve+0*vmu+1*vtau
        y*=0.7
        y+=0.15
        return y

    axis1 = ROOT.TGaxis(
        getX(0,1,0),getY(0,1,0),
        getX(1,0,0),getY(1,0,0),
        0.0,1,
        512,"+U"
    )
    
    ve = "#lbar#kern[-0.7]{ }V#lower[0.3]{#scale[0.8]{eN}}#kern[-0.7]{ }#lbar"
    vmu = "#lbar#kern[-0.7]{ }V#lower[0.3]{#scale[0.8]{#muN}}#kern[-0.7]{ }#lbar"
    vtau = "#lbar#kern[-0.7]{ }V#lower[0.3]{#scale[0.8]{#tauN}}#kern[-0.7]{ }#lbar"
    
    rootObj = []
    
    cv = ROOT.TCanvas("cv","",800,750)
    cv.SetLeftMargin(0.14)
    cv.SetRightMargin(0.04)
    cv.SetBottomMargin(0.12)
    cv.SetTopMargin(0.08)
    cv.Range(0,0,1,1)
    cv.RangeAxis(0,0,1,1)
        
    maxMass = numpy.amax(res[:,3])
    minMass = numpy.amin(res[:,3])
    massRange = maxMass-minMass
    maxMass+=massRange*0.05
    minMass-=massRange*0.05
    
    def getCol(mass):
        return start+int(199.*(mass-minMass)/(maxMass-minMass))
    
        
    for i in range(res.shape[0]):
        
    	marker = ROOT.TMarker(
    	    getX(res[i,0],res[i,1],res[i,2]),
    	    getY(res[i,0],res[i,1],res[i,2]),
    	    20
	    )
        rootObj.append(marker)
        marker.SetMarkerColor(getCol(res[i,3]))
        marker.SetMarkerSize(2)
        marker.Draw()
        #print i,getX(res[i,0],res[i,1],res[i,2]),getY(res[i,0],res[i,1],res[i,2]),x
        
    
    applyAxisStyle(axis1)
    axis1.SetTickSize(0.025)
    axis1.Draw()

    axis1L = ROOT.TGaxis(
        
        getX(0,1,0),getY(0,1,0),
        getX(1,0,0),getY(1,0,0),
        0.0,1,
        512,"-S"
    )
    applyAxisStyle(axis1L)
    axis1L.SetTickSize(0.0)
    axis1L.SetTitleOffset(1.5)
    axis1L.SetTitle("f#lower[0.3]{#scale[0.8]{e#mu}}")##times|V#lower[0.3]{#scale[0.8]{e}}|#lower[-0.7]{#scale[0.8]{2}}+(1-f#lower[0.3]{#scale[0.8]{e#mu}})#times|V#lower[0.3]{#scale[0.8]{#mu}}|#lower[-0.7]{#scale[0.8]{2}}")
    #axis1L.SetTitle(ve+" #lower[0.15]{#scale[1.4]{/}} #lower[0.]{#scale[1.2]{(}}#kern[-0.5]{ }"+ve+" + "+vmu+"#kern[-0.5]{ }#lower[0.]{#scale[1.2]{)}}")##times|V#lower[0.3]{#scale[0.8]{e}}|#lower[-0.7]{#scale[0.8]{2}}+(1-f#lower[0.3]{#scale[0.8]{e#mu}})#times|V#lower[0.3]{#scale[0.8]{#mu}}|#lower[-0.7]{#scale[0.8]{2}}")
    axis1L.SetLabelOffset(0.01)
    axis1L.Draw()


    axis2 = ROOT.TGaxis(  
        getX(0,0,1),getY(0,0,1),
        getX(0,1,0),getY(0,1,0),
        0.0,1,
        512,"+U"
    )
    applyAxisStyle(axis2)
    axis2.SetTickSize(0.025)
    axis2.Draw()

    axis2L = ROOT.TGaxis(
        
        getX(0,0,1),getY(0,0,1),
        getX(0,1,0),getY(0,1,0),
        0.0,1,
        512,"+S="
    )
    applyAxisStyle(axis2L)
    axis2L.SetTickSize(0.0)
    axis2L.SetTitleOffset(1.75)
    axis2L.SetLabelOffset(0.07)
    axis2L.SetTitle("f#lower[0.3]{#scale[0.8]{#mu#tau}}")##times|V#lower[0.3]{#scale[0.8]{#mu}}|#lower[-0.7]{#scale[0.8]{2}}+(1-f#lower[0.3]{#scale[0.8]{#mu#tau}})#times|V#lower[0.3]{#scale[0.8]{#tau}}|#lower[-0.7]{#scale[0.8]{2}}")
    axis2L.Draw()


    axis3 = ROOT.TGaxis(
        getX(1,0,0),getY(1,0,0),
        getX(0,0,1),getY(0,0,1),
        0.0,1,
        512,"+U"
    )
    applyAxisStyle(axis3)
    axis3.SetTickSize(0.025)
    axis3.Draw()

    axis3L = ROOT.TGaxis(
        getX(0,0,1),getY(0,0,1),
        getX(1,0,0),getY(1,0,0),
        
        0.0,1,
        512,"-S="
    )
    applyAxisStyle(axis3L)
    axis3L.SetTickSize(0.0)
    axis3L.SetTitleOffset(1.75)
    axis3L.SetLabelOffset(0.07)
    axis3L.SetTitle("f#lower[0.3]{#scale[0.8]{e#tau}}")##times|V#lower[0.3]{#scale[0.8]{e}}|#lower[-0.7]{#scale[0.8]{2}}+(1-f#lower[0.3]{#scale[0.8]{e#tau}})#times|V#lower[0.3]{#scale[0.8]{#tau}}|#lower[-0.7]{#scale[0.8]{2}}")
    axis3L.Draw()
    
    for i,mass in enumerate(numpy.linspace(minMass,maxMass,100)):
        y1 = 0.4+i*0.5/100.
        y2 = 0.4+(i+1)*0.5/100.
        box = ROOT.TBox(0.89,y1,0.93,y2)
        rootObj.append(box)
        box.SetFillColor(getCol(mass))
        box.Draw("F")
    
    axisColor = ROOT.TGaxis(
        0.93,0.4,
        0.93,0.9,
        minMass,maxMass,
        512,"+U"
    )
    applyAxisStyle(axis3)
    axisColor.SetTickSize(0.025)
    axisColor.Draw()

    axisColorL = ROOT.TGaxis(
        0.93,0.4,
        0.93,0.9,
        minMass,maxMass,
        512,"-S="
    )
    applyAxisStyle(axisColorL)
    axisColorL.SetTickSize(0.0)
    axisColorL.SetTitleOffset(0.5)
    axisColorL.SetLabelOffset(0.05)
    axisColorL.SetTitle("min. m#lower[0.3]{#scale[0.8]{HNL}} (GeV)")##times|V#lower[0.3]{#scale[0.8]{e}}|#lower[-0.7]{#scale[0.8]{2}}+(1-f#lower[0.3]{#scale[0.8]{e#tau}})#times|V#lower[0.3]{#scale[0.8]{#tau}}|#lower[-0.7]{#scale[0.8]{2}}")
    axisColorL.Draw()
    
    
    for v in numpy.linspace(0.1,0.9,9):
    
        l1 = ROOT.TLine(
            getX(1,0,v),getY(1,0,v),
            getX(0,1-v,v),getY(0,1-v,v)
        )
        rootObj.append(l1)
        l1.SetLineColor(ROOT.kGray+1)
        l1.SetLineStyle(2)
        l1.Draw("Same")
        
        l2 = ROOT.TLine(
            getX(v,1-v,0),getY(v,1-v,0),
            getX(v,1-v,v),getY(v,1-v,v)
        )
        rootObj.append(l2)
        l2.SetLineColor(ROOT.kGray+1)
        l2.SetLineStyle(2)
        l2.Draw("Same")
        
        l3 = ROOT.TLine(
            getX(v,1-v,0),getY(v,1-v,0),
            getX(v,0,1-v),getY(v,0,1-v)
        )
        rootObj.append(l3)
        l3.SetLineColor(ROOT.kGray+1)
        l3.SetLineStyle(2)
        l3.Draw("Same")
        
    pText = ROOT.TPaveText(0.03,0.95,0.03,0.95,"NDC")
    pText.SetBorderSize(0)
    pText.SetFillStyle(0)
    pText.SetTextFont(63)
    pText.SetTextAlign(11)
    pText.SetTextSize(40)
    pText.AddText("CMS")
    pText.Draw()
    
    pText2 = ROOT.TPaveText(0.14,0.95,0.14,0.95,"NDC")
    pText2.SetBorderSize(0)
    pText2.SetFillStyle(0)
    pText2.SetTextFont(53)
    pText2.SetTextAlign(11)
    pText2.SetTextSize(40)
    pText2.AddText("Preliminary")
    pText2.Draw()
    
    pText3 = ROOT.TPaveText(0.98,0.95,0.98,0.95,"NDC")
    pText3.SetBorderSize(0)
    pText3.SetFillStyle(0)
    pText3.SetTextFont(43)
    pText3.SetTextAlign(31)
    pText3.SetTextSize(40)
    pText3.AddText("137.1 fb#lower[-0.7]{#scale[0.7]{-1}}")
    pText3.Draw()
    
    pText4 = ROOT.TPaveText(0.03,0.75,0.03,0.88,"NDC")
    pText4.SetBorderSize(0)
    pText4.SetFillStyle(0)
    pText4.SetTextFont(43)
    pText4.SetTextAlign(11)
    pText4.SetTextSize(35)
    pText4.AddText("{} HNL".format(typeHNL))
    pText4.AddText("c#tau = {}#kern[-0.4]{{ }}mm".format(ctau))
    pText4.Draw()
            
    cv.Print("limit_{}_{}.pdf".format(typeHNL, ctau))

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
    transformedPoints = numpy.stack(transformedPoints,axis=1)
    
    delaunay = scipy.spatial.Delaunay(transformedPoints)
    def getValue(point):
        transformedPoint = []
        for i,transformation in enumerate(transformations):
            transformedPoint.append(
                transformation(point[i])
            )
        transformedPoint = numpy.array(transformedPoint)
        simplexIndex = delaunay.find_simplex(transformedPoint)
        pointIndices = delaunay.simplices[simplexIndex]
        simplexPoints = transformedPoints[pointIndices]
        minSimplexPoints = numpy.amin(simplexPoints,axis=0)
        maxSimplexPoints = numpy.amax(simplexPoints,axis=0)
        distance = numpy.fabs(maxSimplexPoints-minSimplexPoints)+1e-6
        weight = numpy.prod(1.0-numpy.fabs(simplexPoints-transformedPoint)/distance,axis=1)
        return 0.9*numpy.sum(values[pointIndices]*weight)/numpy.sum(weight)+0.1*numpy.mean(values[pointIndices])
     
    return getValue

gridpacksTable = json.load(open("/vols/cms/LLP/gridpackLookupTable.json"))

limitsArrays = {
    "ctau":[],
    "mass":[],
    "barycentric":{
        "Ve":[],
        "Vmu":[],
        "Vtau":[]
    },
    "couplings":{
        "Ve":[],
        "Vmu":[],
        "Vtau":[]
    },
    "limit":[]
}

xsecArrays = {
    "ctau":[],
    "mass":[],
    "barycentric":{
        "Ve":[],
        "Vmu":[],
        "Vtau":[]
    },
    "couplings":{
        "Ve":[],
        "Vmu":[],
        "Vtau":[]
    },
    "xsec":[]
}

for sample in sorted(gridpacksTable.keys()):
    if sample.find("dirac")<0 and typeHNL == "Dirac":
        continue
    elif sample.find("majorana")<0 and typeHNL == "Majorana":
        continue
    print("Adding sample {}".format(sample))
    limitFile = json.load(open("jsons/limits_combined_"+sample+".json"))
    #print sample
    for weight in sorted(gridpacksTable[sample]["weights"].keys()):
        xsecArrays['ctau'].append(gridpacksTable[sample]['ctau'])
        xsecArrays['mass'].append(gridpacksTable[sample]['mass'])
        xsecArrays['barycentric']['Ve'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Ve"])
        xsecArrays['barycentric']['Vmu'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Vmu"])
        xsecArrays['barycentric']['Vtau'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Vtau"])
        xsecArrays['couplings']['Ve'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Ve"])
        xsecArrays['couplings']['Vmu'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Vmu"])
        xsecArrays['couplings']['Vtau'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Vtau"])
        xsecArrays['xsec'].append(gridpacksTable[sample]["weights"]["10"]["xsec"]["nominal"])

        weightKey = "%.1f"%float(weight)
        if not limitFile.has_key(weightKey):
            continue
        if not limitFile[weightKey].has_key("exp0"):
            continue
        #print "\t",weight,gridpacksTable[sample]["weights"]["10"]["xsec"]["nominal"],
        limitsArrays['ctau'].append(gridpacksTable[sample]['ctau'])
        limitsArrays['mass'].append(gridpacksTable[sample]['mass'])
        limitsArrays['barycentric']['Ve'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Ve"])
        limitsArrays['barycentric']['Vmu'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Vmu"])
        limitsArrays['barycentric']['Vtau'].append(gridpacksTable[sample]["weights"][weight]["barycentric"]["Vtau"])
        limitsArrays['couplings']['Ve'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Ve"])
        limitsArrays['couplings']['Vmu'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Vmu"])
        limitsArrays['couplings']['Vtau'].append(gridpacksTable[sample]["weights"][weight]["couplings"]["Vtau"])
        limitsArrays['limit'].append(limitFile[weightKey]["exp0"])
        
'''
tree = scipy.spatial.KDTree(numpy.stack([
    xsecArrays['ctau'],
    xsecArrays['mass'],
    xsecArrays['barycentric']['Ve'],
    xsecArrays['barycentric']['Vmu'],
    #xsecArrays['xsec']
],axis=1))
distances,indices = tree.query(numpy.array([6.,11.,0.35,0.5]), k=5)
for i in indices:
    print xsecArrays['ctau'][i],xsecArrays['mass'][i],xsecArrays['barycentric']['Ve'][i],xsecArrays['barycentric']['Vmu'][i],xsecArrays['barycentric']['Vtau'][i]
'''

xsecFct = interpolatedFct(
    numpy.stack([
        xsecArrays['ctau'],
        xsecArrays['mass'],
        xsecArrays['barycentric']['Ve'],
        xsecArrays['barycentric']['Vmu']
    ],axis=1),
    numpy.array(xsecArrays['xsec']),
    transformations = [
        lambda x: numpy.log10(x),
        lambda x: numpy.log10(x),
        lambda x: x,
        lambda x: x,
    ]
)
limitFct = interpolatedFct(
    numpy.stack([
        limitsArrays['ctau'],
        limitsArrays['mass'],
        limitsArrays['barycentric']['Ve'],
        limitsArrays['barycentric']['Vmu']
    ],axis=1),
    numpy.array(limitsArrays['limit']),
    transformations = [
        lambda x: numpy.log10(x),
        lambda x: numpy.log10(x),
        lambda x: x,
        lambda x: x,
    ]
)

def minMass(coords,ctau=1.0):
    
    l1 = coords[0]
    l2 = coords[1]
    excluded = []
    for mass in numpy.linspace(1,20,100):
        point = numpy.array([ctau,mass,l1,l2])
        if xsecFct(point)<limitFct(point):
            excluded.append(mass)
    print coords,min(excluded),max(excluded)
    return min(excluded)#,max(excluded),
    
'''
def maxCtau(coords):
    mass = 10.
    
    l1 = coords[0]
    l2 = coords[1]
    excluded = []
    for ctau in numpy.logspace(-2,2,100):
        point = numpy.array([ctau,mass,l1,l2])
        xsec = xsecFct(point)
        lim = limitFct(point)
        if xsec<lim:
            print xsec,lim,ctau,mass 
            excluded.append(ctau)
    print coords,min(excluded),max(excluded)
    return min(excluded)#,max(excluded),
'''
    
res = []
for _,l1,l2,l3 in generate(41):
    res.append([l1,l2,l3,minMass([l1,l2,l3],ctau = float(ctau))])
   
makePlot(numpy.array(res), ctau=ctau, typeHNL=typeHNL)



