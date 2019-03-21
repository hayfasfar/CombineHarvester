import os
import sys
import numpy
import math
import ROOT
import random
import json
import scipy
import scipy.interpolate

# ROOT.gStyle.SetHistFillStyle(0)
# ROOT.gStyle.SetLegoInnerR(Float_t rad = 0.5)
# ROOT.gStyle.SetNumberContours(Int_t number = 20)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetEndErrorSize(2)
#ROOT.gStyle.SetErrorMarker(20)
ROOT.gStyle.SetErrorX(0.)

ROOT.gStyle.SetMarkerStyle(20)
#ROOT.gStyle.SetMarkerStyle(20)

#For the fit/function:
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetFitFormat("5.4g")
ROOT.gStyle.SetFuncColor(2)
ROOT.gStyle.SetFuncStyle(1)
ROOT.gStyle.SetFuncWidth(1)

#For the date:
ROOT.gStyle.SetOptDate(0)
# ROOT.gStyle.SetDateX(Float_t x = 0.01)
# ROOT.gStyle.SetDateY(Float_t y = 0.01)

# For the statistics box:
ROOT.gStyle.SetOptFile(0)
ROOT.gStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetStatFont(42)
ROOT.gStyle.SetStatFontSize(0.025)
ROOT.gStyle.SetStatTextColor(1)
ROOT.gStyle.SetStatFormat("6.4g")
ROOT.gStyle.SetStatBorderSize(1)
ROOT.gStyle.SetStatH(0.1)
ROOT.gStyle.SetStatW(0.15)

ROOT.gStyle.SetHatchesSpacing(1)
ROOT.gStyle.SetHatchesLineWidth(2)

# ROOT.gStyle.SetStaROOT.TStyle(Style_t style = 1001)
# ROOT.gStyle.SetStatX(Float_t x = 0)
# ROOT.gStyle.SetStatY(Float_t y = 0)


#ROOT.gROOT.ForceStyle(True)
#end modified

# For the Global title:

ROOT.gStyle.SetOptTitle(0)



# ROOT.gStyle.SetTitleH(0) # Set the height of the title box
# ROOT.gStyle.SetTitleW(0) # Set the width of the title box
#ROOT.gStyle.SetTitleX(0.35) # Set the position of the title box
#ROOT.gStyle.SetTitleY(0.986) # Set the position of the title box
# ROOT.gStyle.SetTitleStyle(Style_t style = 1001)
#ROOT.gStyle.SetTitleBorderSize(0)

# For the axis titles:
ROOT.gStyle.SetTitleColor(1, "XYZ")
ROOT.gStyle.SetTitleFont(43, "XYZ")
ROOT.gStyle.SetTitleSize(35, "XYZ")
# ROOT.gStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
# ROOT.gStyle.SetTitleYSize(Float_t size = 0.02)
ROOT.gStyle.SetTitleXOffset(1.2)
#ROOT.gStyle.SetTitleYOffset(1.2)
ROOT.gStyle.SetTitleOffset(1.45, "Y") # Another way to set the Offset
ROOT.gStyle.SetTitleOffset(1.38, "Z")
# For the axis labels:

ROOT.gStyle.SetLabelColor(1, "XYZ")
ROOT.gStyle.SetLabelFont(43, "XYZ")
ROOT.gStyle.SetLabelOffset(0.0077, "XYZ")
ROOT.gStyle.SetLabelSize(29, "XYZ")
#ROOT.gStyle.SetLabelSize(0.04, "XYZ")

# For the axis:

ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.025, "Y")
ROOT.gStyle.SetTickLength(0.025, "X")


ROOT.gStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
ROOT.gStyle.SetPadTickY(1)

# Change for log plots:
ROOT.gStyle.SetOptLogx(0)
ROOT.gStyle.SetOptLogy(0)
ROOT.gStyle.SetOptLogz(0)

#ROOT.gStyle.SetPalette(1) #(1,0)

# another top group addition

# Postscript options:
#ROOT.gStyle.SetPaperSize(20., 20.)
#ROOT.gStyle.SetPaperSize(ROOT.TStyle.kA4)
#ROOT.gStyle.SetPaperSize(27., 29.7)
#ROOT.gStyle.SetPaperSize(27., 29.7)
ROOT.gStyle.SetPaperSize(8.5*1.4,7.0*1.4)
ROOT.TGaxis.SetMaxDigits(3)
ROOT.gStyle.SetLineScalePS(2)

# ROOT.gStyle.SetLineStyleString(Int_t i, const char* text)
# ROOT.gStyle.SetHeaderPS(const char* header)
# ROOT.gStyle.SetTitlePS(const char* pstitle)
#ROOT.gStyle.SetColorModelPS(1)

# ROOT.gStyle.SetBarOffset(Float_t baroff = 0.5)
# ROOT.gStyle.SetBarWidth(Float_t barwidth = 0.5)
# ROOT.gStyle.SetPaintTextFormat(const char* format = "g")
# ROOT.gStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
# ROOT.gStyle.SetTimeOffset(Double_t toffset)
# ROOT.gStyle.SetHistMinimumZero(kTRUE)

ROOT.gStyle.SetPaintTextFormat("3.0f")

NRGBs = 7;
NCont = 255;

stops = numpy.array( [0.0,0.15, 0.4, 0.55,0.65, 0.88, 1.00] )
red  = numpy.array( [0.6,0.00, 0.00,0.3, 0.87, 1.00, 0.51] )
green = numpy.array( [0.0,0.00, 0.81,0.8, 1.00, 0.20, 0.00] )
blue = numpy.array( [0.8,0.51, 1.00,0.1, 0.12, 0.00, 0.00] )

colWheelDark = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)

for i in range(NRGBs):
    red[i]=min(1,red[i]*1.1+0.25)
    green[i]=min(1,green[i]*1.1+0.25)
    blue[i]=min(1,blue[i]*1.1+0.25)

colWheel = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
ROOT.gStyle.SetNumberContours(NCont)
'''

NRGBs = 3;
NCont = 8;

stops = numpy.array( [0.0, 0.5, 1.00] )
red  = numpy.array( [0.05,  0.1,  0.95] )
green = numpy.array( [0.35, 0.95,  0.75] )
blue = numpy.array( [0.95,  0.75, 0.05] )

colWheelDark = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)


colWheel = ROOT.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)
ROOT.gStyle.SetNumberContours(NCont)

'''
'''
stops = numpy.array( [0.00, 0.20, 0.4, 0.6, 0.8, 1.00] )
red  = numpy.array( [163/255., 123/255.,  96/255., 145/255., 221/255., 253/255.] )
green = numpy.array( [33/255.,  77/255., 180/255., 242/255., 223/255., 217/255.] )
blue = numpy.array([254./255., 251/255., 237/255., 180/255., 165/255., 209/255.] )
colWheel = ROOT.TColor.CreateGradientColorTable(6, stops, red, green, blue, 255)
ROOT.gStyle.SetNumberContours(255)
'''
ROOT.gRandom.SetSeed(123)

colors=[]
def hex2rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16)/255.0 for i in range(0, lv, lv // 3))

def newColor(red,green,blue):
    newColor.colorindex+=1
    color=ROOT.TColor(newColor.colorindex,red,green,blue)
    colors.append(color)
    return color
    
newColor.colorindex=301

def getDarkerColor(color):
    darkerColor=newColor(color.GetRed()*0.6,color.GetGreen()*0.6,color.GetBlue()*0.6)
    return darkerColor

with open('eventyields.json',) as f:
    genweights = json.load(f)

ctauValues = ["0p001","0p01","0p1","1","10","100","1000","10000"]
#ctauValues = ["10"]
ctauLabels = {
    "0p001":"1#kern[-0.5]{ }#mum",
    "0p01":"10#kern[-0.5]{ }#mum",
    "0p1":"100#kern[-0.5]{ }#mum",
    "1":"1#kern[-0.5]{ }mm",
    "10":"10#kern[-0.5]{ }mm",
    "100":"100#kern[-0.5]{ }mm",
    "1000":"1#kern[-0.5]{ }m",
    "10000":"10#kern[-0.5]{ }m",
}
ctauValueMap = {
    "0p001":0,
    "0p01":1,
    "0p1":2,
    "1":3,
    "10":4,
    "100":5,
    "1000":6,
    "10000":7,
}

massesDict = {}
for ctau in ctauValues:
    ctauSampleName = "SMS-T1qqqq_ctau-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"%ctau
    if not massesDict.has_key(ctau):
        massesDict[ctau] = {}
    for signalSample in [ctauSampleName,ctauSampleName+"_extra"]:
        for llpMass in genweights[signalSample].keys():
            for lspMass in genweights[signalSample][llpMass].keys():
                if genweights[signalSample][llpMass][lspMass]["sum"]<10:
                    continue
                if int(llpMass)==int(lspMass):
                    continue
                if not massesDict[ctau].has_key(llpMass):
                    massesDict[ctau][llpMass] = []
                if not lspMass in massesDict[ctau][llpMass]:
                    massesDict[ctau][llpMass].append(lspMass)

def getTheoryXsecFct(filePath):
    f = open(filePath)
    llpvalues = []
    xsecvalues = []
    xsecvaluesUp = []
    xsecvaluesDown = []
    for l in f:
        if len(l)==0:
            continue
        splitted = l.split(",")
        if len(splitted)!=3:
            print "cannot parse xsec line: ",l
            continue
        llpvalues.append(float(splitted[0]))
        xsec = float(splitted[1])
        xsecRelUnc = float(splitted[2])
        xsecvalues.append(math.log(xsec))
        xsecvaluesUp.append(math.log(xsec*(1+xsecRelUnc)))
        xsecvaluesDown.append(math.log(xsec*(1-xsecRelUnc)))
        
    llpvalues = numpy.array(llpvalues,dtype=numpy.float32)
    xsecvalues = numpy.array(xsecvalues,dtype=numpy.float32)
    xsecvaluesUp = numpy.array(xsecvaluesUp,dtype=numpy.float32)
    xsecvaluesDown = numpy.array(xsecvaluesDown,dtype=numpy.float32)
    
    tckXsec = scipy.interpolate.splrep(llpvalues,xsecvalues,s=1e-3)
    tckXsecUp = scipy.interpolate.splrep(llpvalues,xsecvaluesUp,s=1e-3)
    tckXsecDown = scipy.interpolate.splrep(llpvalues,xsecvaluesDown,s=1e-3)

    def getValue(llpMass):
        xsec = math.exp(scipy.interpolate.splev(llpMass,tckXsec))
        xsecUp = math.exp(scipy.interpolate.splev(llpMass,tckXsecUp))
        xsecDown = math.exp(scipy.interpolate.splev(llpMass,tckXsecDown))
        return xsec,xsecUp,xsecDown
    return getValue
    
    
def interpolatedHist(limitFct,binningX,binningY):
    newHist = ROOT.TH2F(
        "interpolated"+str(random.random()),
        "",
        len(binningX)-1,
        binningX,
        len(binningY)-1,
        binningY
    )
    newHist.SetDirectory(0)
    
    for ibin in range(newHist.GetNbinsX()):
        llpMass = newHist.GetXaxis().GetBinCenter(ibin+1)*1000
        for jbin in range(newHist.GetNbinsY()):
            lspMass = newHist.GetYaxis().GetBinCenter(jbin+1)*1000
            if (llpMass-lspMass)<100:
                continue
            newHist.SetBinContent(ibin+1,jbin+1,limitFct(llpMass,lspMass))
    
    return newHist 

def interpolatedLimitFct(result,kind="median"):
    xvalues = []
    yvalues = []
    zvalues = []
    for llpMass in sorted(result.keys()):
        for lspMass in sorted(result[llpMass].keys()):
            loglimit = math.log(result[llpMass][lspMass][kind])
            xvalues.append(int(llpMass))
            yvalues.append(math.sqrt(100+int(llpMass)-int(lspMass)))
            zvalues.append(loglimit)
    xvalues = numpy.array(xvalues,dtype=numpy.float32)
    yvalues = numpy.array(yvalues,dtype=numpy.float32)
    zvalues = numpy.array(zvalues,dtype=numpy.float32)
   
    tck = scipy.interpolate.bisplrep(
        xvalues,
        yvalues,
        zvalues,
        kx=3,
        ky=3,
        s=1e-3
    )
       
    def getValue(llpMass,lspMass):
        tckCopy = numpy.copy(tck)
        return numpy.exp(scipy.interpolate.bisplev(
            1.*int(llpMass),1.*math.sqrt(100+int(llpMass)-int(lspMass)),tckCopy
        ))
        
    n = 0
    meanDiff = 0.
    meanDiff2 = 0.
    maxDiff = 0.
    for llpMass in sorted(result.keys()):
        for lspMass in sorted(result[llpMass].keys()):
            n+=1
            limit = (result[llpMass][lspMass][kind])
            limitSmooth = getValue(1.*int(llpMass),1.*int(lspMass))
            diff = 1-limitSmooth/limit
            meanDiff += math.fabs(diff)
            meanDiff2 += math.fabs(diff)**2
            maxDiff = max(maxDiff, math.fabs(diff))
    print "rel. interpolation difference mean: %5.3f+-%.3f (max: %5.3f)"%(meanDiff/n,math.sqrt(meanDiff2/n-(meanDiff/n)**2),maxDiff)
 
        
    return getValue
    
def interpolatedLimitFct2(result,kind="median"):
    xvalues = []
    yvalues = []
    zvalues = []
    for llpMass in sorted(result.keys()):
        for lspMass in sorted(result[llpMass].keys()):
            loglimit = math.log(result[llpMass][lspMass][kind])
            xvalues.append(int(llpMass))
            yvalues.append((100+int(llpMass)-int(lspMass)))
            zvalues.append(loglimit)
    xvalues = numpy.array(xvalues,dtype=numpy.float32)
    yvalues = numpy.array(yvalues,dtype=numpy.float32)
    zvalues = numpy.array(zvalues,dtype=numpy.float32)
   
       
    def getValue(llpMass,lspMass):
        xval = 1.*int(llpMass)
        yval = (100+int(llpMass)-int(lspMass))
        
        distances = numpy.zeros(len(xvalues))
        values = numpy.zeros(len(xvalues))
        
        for i in range(len(xvalues)):
            distance = math.sqrt((xvalues[i]-xval)**2+(yvalues[i]-yval)**2)
            distances[i] = distance
            values[i] = zvalues[i]
        
        sortedIndices = numpy.argsort(distances)
        distances = distances[sortedIndices]
        values = values[sortedIndices]
        #print distances
        #print values
        wsum = 0.
        vsum = 0.
        
        wsumExclNearest = 0.
        vsumExclNearest = 0.
        
        for i in range(len(distances)):
            weight = math.exp(-distances[i]**1.7/((20+numpy.mean(distances[0:2]))**1.7))
            vsum += weight*values[i]
            wsum+=weight
            
            if i>=1:
                weightExclNearest = math.exp(-distances[i]**1.7/((20+numpy.mean(distances[1:2]))**1.7))
                wsumExclNearest += weightExclNearest*values[i]
                vsumExclNearest+=weightExclNearest
                
                
        result = math.exp(vsum/wsum)
        #resultExclNearest = math.exp(vsumExclNearest/wsumExclNearest)
        
        
        return result#(distances[1]*resultExclNearest+distances[0]*result)/(distances[1]+distances[0])
        
        '''
        nearestX = (numpy.abs(xvalues[1:-1]-xval)).argmin()+1
        nearestY = (numpy.abs(yvalues[1:-1]-yval)).argmin()+1
        
        
            
        if xvalues[nearestX]>xval:
            leftX = xvalues[nearestX-1]
            rightX = xvalues[nearestX]
            zleftX = zvalues[nearestX-1]
            zrightX = zvalues[nearestX]
        else:
            leftX = xvalues[nearestX]
            rightX = xvalues[nearestX+1]
            zleftX = zvalues[nearestX]
            zrightX = zvalues[nearestX+1]
            
        if yvalues[nearestY]>yval:
            leftY = yvalues[nearestY-1]
            rightY = yvalues[nearestY]
            zleftY = zvalues[nearestY-1]
            zrightY = zvalues[nearestY]
        else:
            leftY = yvalues[nearestY]
            rightY = yvalues[nearestY+1]
            zleftY = zvalues[nearestY]
            zrightY = zvalues[nearestY+1]
        
        xfrac = (xval-leftX)/(rightX-leftX)
        yfrac = (yval-leftY)/(rightY-leftY)
        
        print yval,leftY,rightY
        
        return 0.5*(xfrac*zleftX+(1-xfrac)*zrightX+zleftY*yfrac+(1-yfrac)*zrightY)
        '''
        
    n = 0
    meanDiff = 0.
    meanDiff2 = 0.
    maxDiff = 0.
    for llpMass in sorted(result.keys()):
        for lspMass in sorted(result[llpMass].keys()):
            n+=1
            limit = (result[llpMass][lspMass][kind])
            limitSmooth = getValue(1.*int(llpMass),1.*int(lspMass))
            diff = 1-limitSmooth/limit
            meanDiff += math.fabs(diff)
            meanDiff2 += math.fabs(diff)**2
            maxDiff = max(maxDiff, math.fabs(diff))
    print "rel. interpolation difference mean: %5.3f+-%.3f (max: %5.3f)"%(meanDiff/n,math.sqrt(meanDiff2/n-(meanDiff/n)**2),maxDiff)
 
        
    return getValue
    
    
def parseCombineJson(filePath):
    if not os.path.exists(filePath):
        return None    
    data = json.load(open(filePath))
    mapping = {
        "median": "exp0", 
        "+1":"exp+1", 
        "-1":"exp-1",
        "+2":"exp+2",
        "-2":"exp-2"
    }
    result = {}
    for k,v in mapping.iteritems():
        if not data.has_key("120.0"):
            return None
        if data["120.0"].has_key(v):
            result[k]=data["120.0"][v]
        
    return result

def parseCombineResult(filePath):
    rootFile = ROOT.TFile(filePath)
    limitTree = rootFile.Get("limit")
    result = {}
    #note: combine seems to be not very precise here
    mapping = {
        "median": 0.5, 
        "+1":0.840, 
        "-1":0.160,
        "+2":0.975,
        "-2":0.025
    }
    for entry in range(limitTree.GetEntries()):
        limitTree.GetEntry(entry)
        found = False
        for k,v in mapping.iteritems():
            if math.fabs(limitTree.quantileExpected-v)<0.0001:
                result[k] = limitTree.limit
                found = True
                break
        if not found and limitTree.quantileExpected>0:
            print "unknown quantile: ",limitTree.quantileExpected

    return result

basePath = "limits"

results = {}



for ctau in ctauValues:
    results[ctau] = {}
    for llpMass in sorted(massesDict[ctau].keys()):
        results[ctau][llpMass] = {}
        for lspMass in massesDict[ctau][llpMass]:
            signalProcess = "ctau%s_llp%s_lsp%s"%(ctau,llpMass,lspMass)
            combineJsonFile = os.path.join(basePath,"limits_%s.json"%signalProcess)
            result = parseCombineJson(combineJsonFile)
            if not result:
                print "ERROR: file '"+combineJsonFile+"' not found -> skip"
                continue
            if len(result.keys())!=5:
                print "WARNING: Not all quantiles found in file ",combineJsonFile
            if not result.has_key("median"):
                print "ERROR: Median found in file ",combineJsonFile," -> skip"
                continue
            results[ctau][llpMass][lspMass] = result
                
   
limitsU = {}
limitsC = {}
   
for ctau in results.keys():
    cv = ROOT.TCanvas("cv"+ctau,"",850,700)
    cv.SetLeftMargin(0.145)
    cv.SetRightMargin(0.195)
    cv.SetBottomMargin(0.14)
    cv.SetTopMargin(0.08)
    cv.SetLogz(1)
    xmin = 0.600
    ymin = 0
    xmax = 2.400
    ymax = 2.400
    zmax = 0.35
    zmin = 0.0005

    axis = ROOT.TH2F("axis"+ctau,";m#lower[0.2]{#scale[0.8]{#tilde{g}}} (TeV); m#lower[0.2]{#scale[0.8]{#tilde{#chi}#lower[-0.5]{#scale[0.65]{0}}#kern[-1.2]{#lower[0.6]{#scale[0.65]{1}}}}} (TeV)",
        int(round((xmax-xmin)/0.050))+1,numpy.linspace(xmin-0.025,xmax+0.025,int(round((xmax-xmin)/0.050))+2),
        int(round((ymax-ymin)/0.050))+1,numpy.linspace(ymin-0.025,ymax+0.025,int(round((ymax-ymin)/0.050))+2)
    )
    axis.Fill(-1,-1)
    '''
    for xbin in range(axis.GetNbinsX()):
        value = xmin+xbin*(xmax-xmin)/50
        if xbin%4==0:
            axis.GetXaxis().SetBinLabel(xbin+1,"%.0f"%axis.GetXaxis().GetBinCenter(xbin+1))
    '''
    axis.Draw("colz")    
    axis.GetZaxis().SetTitle("95% CL expected limit #sigma#lower[0.2]{#scale[0.8]{pp#rightarrow#tilde{g}#tilde{g}}} (pb)")
    axis.GetZaxis().SetRangeUser(zmin,zmax)
    
    axis.GetXaxis().SetNoExponent(True)
    #axis.GetXaxis().LabelsOption("v")
    axis.GetYaxis().SetNoExponent(True)
    axis.GetXaxis().SetNdivisions(510)
    axis.GetYaxis().SetNdivisions(510)
    
    
    
    limitHist = ROOT.TH2F("limitHist"+ctau,";m#lower[0.2]{#scale[0.8]{#tilde{g}}} (TeV); m#lower[0.2]{#scale[0.8]{#tilde{#chi}#lower[-0.5]{#scale[0.65]{0}}#kern[-1.2]{#lower[0.6]{#scale[0.65]{1}}}}} (TeV)",
        int(round((xmax-xmin)/0.050))+1,numpy.linspace(xmin-0.025,xmax+0.025,int(round((xmax-xmin)/0.050))+2),
        int(round((ymax-ymin)/0.050))+1,numpy.linspace(ymin-0.025,ymax+0.025,int(round((ymax-ymin)/0.050))+2)
    )
    boxes = []
    for llpMass in sorted(results[ctau].keys()):
        for lspMass in sorted(results[ctau][llpMass].keys()):
        
            if 0.001*int(llpMass)>xmax:
                continue
            if 0.001*int(lspMass)>ymax:
                continue
            limitHist.Fill(0.001*int(llpMass),0.001*int(lspMass),results[ctau][llpMass][lspMass]["median"])
            '''
            box= ROOT.TBox(0.001*int(llpMass)-0.025,0.001*int(lspMass)-0.025,0.001*int(llpMass)+0.025,0.001*int(lspMass)+0.025)
            box.SetFillColor(ROOT.kGray)
            box.SetLineWidth(0)
            box.SetFillStyle(1001)
            boxes.append(box)
            '''
            
            marker= ROOT.TMarker(0.001*int(llpMass),0.001*int(lspMass),20)
            marker.SetMarkerColor(ROOT.kWhite)
            marker.SetMarkerSize(1.4)
            boxes.append(marker)
            
    
    limitFct = interpolatedLimitFct2(
        results[ctau],
        kind="median"
    )
    
    limitFctUp = interpolatedLimitFct2(
        results[ctau],
        kind="+1"
    )
    limitFctDown = interpolatedLimitFct2(
        results[ctau],
        kind="-1"
    )
    
    #print limitFct(600,0)
    
    
    limitHistSmooth = interpolatedHist(
        limitFct,
        numpy.linspace(xmin-0.025,xmax+0.025,(xmax-xmin)/0.010),
        numpy.linspace(ymin-0.025,ymax+0.025,(ymax-ymin)/0.010)
    )
    
    limitHistSmooth.GetZaxis().SetRangeUser(zmin,zmax)
    limitHistSmooth.Draw("colSame")
    
    limitHist.GetZaxis().SetRangeUser(zmin,zmax)
    #limitHist.Draw("colSame")
    
    for box in boxes:
        box.Draw("L")
    
    
    poly = ROOT.TPolyLine(3,
        numpy.array([0.600-0.025,2.400+0.025,0.600-0.025],dtype=numpy.float32), 
        numpy.array([0.600-0.025,2.400+0.025,2.400+0.025],dtype=numpy.float32),
    )
    poly.SetFillColor(ROOT.kGray)
    poly.SetFillStyle(3445)
    poly.Draw("F")
    
    
    
    
    xsecFct = getTheoryXsecFct("theory_xsec.dat")
    
    
    llpMassExpMedian = []
    lspMassExpMedian = []
    
    llpMassExpUp = []
    lspMassExpUp = []
    
    llpMassExpDown = []
    lspMassExpDown = []
    
    for angle in numpy.linspace(0.0,math.pi/4,50):
        foundDown = False
        foundMedian = False
        foundUp = False
        for r in numpy.linspace(2.600,0.600,100):
            llpMass = r*math.cos(angle)
            lspMass = r*math.sin(angle)
            #print llpMass,lspMass
            if llpMass>(xmax) or lspMass>(ymax):
                continue
            if (llpMass-lspMass)<0.1:
                continue
            xsecTheo,xsecTheoUp,xsecTheoDown = xsecFct(llpMass*1000.)
            xsecLimit = limitFct(llpMass*1000.,lspMass*1000.)
            #print llpMass,lspMass,xsecLimit,xsecTheo
            xsecLimitUp = limitFctUp(llpMass*1000.,lspMass*1000.)
            xsecLimitDown = limitFctDown(llpMass*1000.,lspMass*1000.)
            if not foundDown and xsecLimitDown<xsecTheo:
                llpMassExpDown.append(llpMass)
                lspMassExpDown.append(lspMass)
                foundDown = True
            if not foundMedian and xsecLimit<xsecTheo:
                llpMassExpMedian.append(llpMass)
                lspMassExpMedian.append(lspMass)
                foundMedian = True
            if not foundUp and xsecLimitUp<xsecTheo:
                llpMassExpUp.append(llpMass)
                lspMassExpUp.append(lspMass)
                foundUp = True
            if foundDown and foundMedian and foundUp:
                break
    #print llpMassExpMedian
    #print lspMassExpMedian
                
    foundC=False
    foundU=False
    for r in numpy.linspace(3.000,0.600,200):
        llpMassC = r
        lspMassC = r-0.100
        llpMassU = r
        lspMassU = 0.100
        
        xsecTheoC,_,_ = xsecFct(llpMassC*1000.)
        xsecLimitC = limitFct(llpMassC*1000.,lspMassC*1000.)
        if not foundC and xsecLimitC<xsecTheoC:
            limitsC[ctau]=llpMassC
            foundC = True
        
        xsecTheoU,_,_ = xsecFct(llpMassU*1000.)
        xsecLimitU = limitFct(llpMassU*1000.,lspMassU*1000.)
        if not foundU and xsecLimitU<xsecTheoU:
            limitsU[ctau]=llpMassU
            foundU = True
            
        if foundC and foundU:
            break
        
    if not foundC:
        limitsC[ctau]=0
    if not foundU:
        limitsU[ctau]=0
        
    llpMassExpMedian = numpy.array(llpMassExpMedian)
    lspMassExpMedian = numpy.array(lspMassExpMedian)
    
    llpMassExpUp = numpy.array(llpMassExpUp)
    lspMassExpUp = numpy.array(lspMassExpUp)
    
    llpMassExpDown = numpy.array(llpMassExpDown)
    lspMassExpDown = numpy.array(lspMassExpDown)
    
    polyExpMedian = ROOT.TPolyLine(
        len(llpMassExpMedian),
        llpMassExpMedian,
        lspMassExpMedian
    )
    polyExpMedian.SetLineColor(ROOT.kBlack)
    polyExpMedian.SetLineWidth(3)
    polyExpMedian.Draw("L")
    
    polyExpDown = ROOT.TPolyLine(
        len(llpMassExpDown),
        llpMassExpDown,
        lspMassExpDown
    )
    polyExpDown.SetLineColor(ROOT.kBlack)
    polyExpDown.SetLineWidth(2)
    polyExpDown.SetLineStyle(2)
    polyExpDown.Draw("L")
    
    polyExpUp = ROOT.TPolyLine(
        len(llpMassExpUp),
        llpMassExpUp,
        lspMassExpUp
    )
    polyExpUp.SetLineColor(ROOT.kBlack)
    polyExpUp.SetLineWidth(2)
    polyExpUp.SetLineStyle(2)
    polyExpUp.Draw("L")
    
    ROOT.gPad.RedrawAxis()
    
    pTextCMS = ROOT.TPaveText(cv.GetLeftMargin(),1-cv.GetTopMargin()*0.4,cv.GetLeftMargin(),1-cv.GetTopMargin()*0.4,"NDC")
    pTextCMS.AddText("CMS")
    pTextCMS.SetTextFont(63)
    pTextCMS.SetTextSize(32)
    pTextCMS.SetTextAlign(13)
    pTextCMS.Draw("Same")
    
    pTextPreliminary = ROOT.TPaveText(cv.GetLeftMargin()+0.08,1-cv.GetTopMargin()*0.4,cv.GetLeftMargin()+0.08,1-cv.GetTopMargin()*0.4,"NDC")
    pTextPreliminary.AddText("Preliminary")
    pTextPreliminary.SetTextFont(53)
    pTextPreliminary.SetTextSize(32)
    pTextPreliminary.SetTextAlign(13)
    pTextPreliminary.Draw("Same")
    
    pTextCTau = ROOT.TPaveText(1-cv.GetRightMargin(),1-cv.GetTopMargin()*0.4,1-cv.GetRightMargin(),1-cv.GetTopMargin()*0.4,"NDC")
    pTextCTau.AddText("c#tau#kern[-0.5]{ }=#kern[-0.8]{ }%s"%ctauLabels[ctau])
    pTextCTau.SetTextFont(43)
    pTextCTau.SetTextSize(32)
    pTextCTau.SetTextAlign(33)
    pTextCTau.Draw("Same")
    
    cv.Print("limits_ctau%s.pdf"%ctau)
    cv.Print("limits_ctau%s.png"%ctau)
    
xvalues = numpy.logspace(-6,1,num=8)
#xvalues = numpy.power(10,xvalues)
yvaluesC = numpy.zeros(8)
yvaluesU = numpy.zeros(8)

for ctau in limitsU.keys():
    yvaluesC[ctauValueMap[ctau]]=limitsC[ctau]
    yvaluesU[ctauValueMap[ctau]]=limitsU[ctau]
    
    print ctau,limitsC[ctau],limitsU[ctau]
    
cv = ROOT.TCanvas("summary","",800,650)
cv.SetGridx(True)
cv.SetGridy(True)
cv.SetLogx(1)
cv.SetMargin(0.155,0.04,0.15,0.09)
ROOT.gStyle.SetGridColor(ROOT.kBlack)
ROOT.gStyle.SetGridStyle(2)
ROOT.gStyle.SetGridWidth(1)
axis = ROOT.TH2F("axis",";c#tau (m); 95% CL lower limit on m#lower[0.2]{#scale[0.8]{#tilde{g}}} (TeV)",50,10**-6.2,10**1.2,50,0.5,2.6)
axis.Draw("AXIS")
graphC = ROOT.TGraph(8,xvalues,yvaluesC)
graphC.SetLineColor(ROOT.kOrange+7)
graphC.SetLineWidth(3)
graphC.SetLineStyle(2)
graphC.SetMarkerStyle(20)
graphC.SetMarkerSize(1.5)
graphC.SetMarkerColor(ROOT.kOrange+7)
graphC.Draw("PL")
graphU = ROOT.TGraph(8,xvalues,yvaluesU)
graphU.SetLineColor(ROOT.kRed+1)
graphU.SetLineWidth(2)
graphU.SetLineStyle(1)
graphU.SetMarkerStyle(20)
graphU.SetMarkerSize(1.5)
graphU.SetMarkerColor(ROOT.kRed+1)
graphU.Draw("PL")

pTextCMS = ROOT.TPaveText(cv.GetLeftMargin(),1-cv.GetTopMargin()+0.055,cv.GetLeftMargin(),1-cv.GetTopMargin()+0.055,"NDC")
pTextCMS.AddText("CMS")
pTextCMS.SetTextFont(63)
pTextCMS.SetTextSize(32)
pTextCMS.SetTextAlign(13)
pTextCMS.Draw("Same")

pTextPreliminary = ROOT.TPaveText(cv.GetLeftMargin()+0.09,1-cv.GetTopMargin()+0.055,cv.GetLeftMargin()+0.09,1-cv.GetTopMargin()+0.05,"NDC")
pTextPreliminary.AddText("Preliminary")
pTextPreliminary.SetTextFont(53)
pTextPreliminary.SetTextSize(32)
pTextPreliminary.SetTextAlign(13)
pTextPreliminary.Draw("Same")

pInfo = ROOT.TPaveText(1-cv.GetRightMargin(),1-cv.GetTopMargin()+0.06,1-cv.GetRightMargin(),1-cv.GetTopMargin()+0.06,"NDC")
pInfo.AddText("p#kern[-0.6]{ }p#kern[-0.5]{ }#rightarrow#kern[-0.5]{ }#tilde{g}#kern[-0.6]{ }#tilde{g}, #tilde{g}#kern[-0.5]{ }#rightarrow#kern[-0.5]{ }q#kern[-0.6]{ }#bar{q}#kern[-0.6]{ }#tilde{#chi}#lower[-0.5]{#scale[0.65]{0}}#kern[-1.2]{#lower[0.6]{#scale[0.65]{1}}}")
pInfo.SetTextFont(43)
pInfo.SetTextSize(32)
pInfo.SetTextAlign(33)
pInfo.Draw("Same")

legend = ROOT.TLegend(1-cv.GetRightMargin()-0.03-0.45,cv.GetBottomMargin()+0.03+0.17,1-cv.GetRightMargin()-0.08,cv.GetBottomMargin()+0.05)
legend.SetFillColor(ROOT.kWhite)
legend.SetBorderSize(1)
legend.SetTextFont(43)
legend.SetTextSize(30)
legend.AddEntry(graphU,"#kern[-0.5]{ }m#lower[0.2]{#scale[0.8]{#tilde{#chi}#lower[-0.5]{#scale[0.65]{0}}#kern[-1.2]{#lower[0.6]{#scale[0.65]{1}}}}}#kern[-0.4]{ }=#kern[-0.5]{ }100#kern[-0.6]{ }GeV","PL")
legend.AddEntry(graphC,"m#lower[0.2]{#scale[0.8]{#tilde{g}}}#kern[-0.5]{ }-#kern[-0.5]{ }m#lower[0.2]{#scale[0.8]{#tilde{#chi}#lower[-0.5]{#scale[0.65]{0}}#kern[-1.2]{#lower[0.6]{#scale[0.65]{1}}}}}#kern[-0.4]{ }=#kern[-0.5]{ }100#kern[-0.6]{ }GeV","PL")
legend.Draw("Same")

cv.Print("summary.pdf")
cv.Print("summary.png")



