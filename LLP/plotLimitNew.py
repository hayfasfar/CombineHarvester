import ROOT
import json
from array import array
import scipy
import scipy.spatial
from scipy import interpolate
from sklearn.neighbors import NearestNeighbors
import numpy as np
import math
import pandas as pd
import os
import random

import style
import refLimits

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPaperSize(7.0*1.35,6.7*1.35)

    
def parse_lookup_table(f, lookup_table, couplingWeight):

    # parse lookup table
    proc = f.replace(year+"_", "").replace("limits_", "").replace(".json", "")
    lu_infos = lookup_table[proc.replace('pt20', 'all')]['weights'][str(int(couplingWeight))]
    xsec = lu_infos['xsec']['nominal']*K_FACTOR
    coupling = lu_infos['couplings']['Ve']+lu_infos['couplings']['Vmu']+lu_infos['couplings']['Vtau']
    if coupling not in (2, 12, 67):
        coupling = coupling/2
    coupling = coupling ** 2
    mass = lookup_table[proc.replace('pt20', 'all')]['mass']

    return mass, coupling, xsec

def parse_limit_json(f, couplingWeight="12.0"):
    # limit json aggreggate 
    with open(os.path.join(json_path, f)) as json_file:
        xsec_dict = json.load(json_file)
    if str(couplingWeight) not in xsec_dict.keys():
        return None
    xsec_dict = xsec_dict[str(couplingWeight)]
    # means something failed with combine (shouldn't happend!)
    if "exp0" not in xsec_dict or "exp+1" not in xsec_dict or "exp+2" not in xsec_dict or "exp-1" not in xsec_dict or "exp-2" not in xsec_dict:
        return None

    return xsec_dict


json_path = "/vols/cms/hsfar/jsons"

K_FACTOR = 1.10

blinded = False

limit_var_list = ["exp0", "exp+1", "exp+2", "exp-1", "exp-2"]

if not blinded:
    limit_var_list += ["obs"]

lumi = {"2016": "36.3", "2017": "41.5", "2018": "59.7", "combined": "138"}
#years = ["2016", "2017", "2018", "combined"]
years = ["combined"]
#years = ["2016"]

veSym = "V#lower[0.3]{#scale[0.7]{e#kern[-0.7]{ }N}}"
vmuSym = "V#lower[0.3]{#scale[0.7]{#mu#kern[-0.7]{ }N}}"
vtauSym = "V#lower[0.3]{#scale[0.7]{#tau#kern[-0.7]{ }N}}"
vlSym = "V#lower[0.3]{#font[12]{l}#kern[-0.7]{ }#scale[0.7]{N}}"

ve2Sym = "|"+veSym+"|#lower[-0.9]{#scale[0.7]{2}}"
vmu2Sym = "|"+vmuSym+"|#lower[-0.9]{#scale[0.7]{2}}"
vtau2Sym = "|"+vtauSym+"|#lower[-0.9]{#scale[0.7]{2}}"
vl2Sym = "|"+vlSym+"|#lower[-0.9]{#scale[0.7]{2}}"

mHNLSym = "m#lower[0.3]{#scale[0.7]{N}}"


limitGreen = style.newColorHLS(0.34,0.45,0.95)
limitYellow = style.newColorHLS(0.15,0.68,0.98)

def scenarioTitle(Ve,Vmu,Vtau):
    Ve = str(Ve)
    Vmu = str(Vmu)
    Vtau = str(Vtau)
    
    
    if Ve=="0":
        Ve = "#kern[-0.8]{ }"+Ve
    if Vmu=="0":
        Vmu = "#kern[-0.8]{ }"+Vmu
    if Vtau=="0":
        Vtau = "#kern[-0.8]{ }"+Vtau

    return veSym+"#kern[-0.5]{ }:#kern[-0.5]{ }"+vmuSym+"#kern[-0.4]{ }:#kern[-0.4]{ }"+vtauSym+"#kern[0.1]{ }=#kern[0.1]{ }"+str(Ve)+"#kern[-0.5]{ }:#kern[-0.5]{ }"+str(Vmu)+"#kern[-0.5]{ }:#kern[-0.5]{ }"+str(Vtau)

scenarios = {
    'emutau': {
        'title': scenarioTitle(1,1,1),
        'ylabel': ve2Sym+" = "+vmu2Sym+" = "+vtau2Sym,
        'couplingWeight':1.0,
        'couplingFct': lambda couplingSum2Value: (math.sqrt(couplingSum2Value)/3.)**2
    },
    'ee': {
        'title': scenarioTitle(1,0,0),
        'ylabel': ve2Sym,
        'couplingWeight':2.0,
        'couplingFct': lambda couplingSum2Value: couplingSum2Value
    },
    'etau': {
        'title': scenarioTitle(1,0,1),
        'ylabel': ve2Sym+" = "+vtau2Sym,
        'couplingWeight':47.0,
        'couplingFct': lambda couplingSum2Value: (math.sqrt(couplingSum2Value)/2.)**2
    },
    'tautau': {
        'title': scenarioTitle(0,0,1),
        'ylabel': vtau2Sym,
        'couplingWeight':67.0,
        'couplingFct': lambda couplingSum2Value: (math.sqrt(couplingSum2Value)/2.)**2
    },
    'emu': {
        'title': scenarioTitle(1,1,0),
        'ylabel': ve2Sym+" = "+vmu2Sym,
        'couplingWeight':7.0,
        'couplingFct': lambda couplingSum2Value: (math.sqrt(couplingSum2Value)/2.)**2
    },
    'mumu': {
        'title': scenarioTitle(0,1,0),
        'ylabel': vmu2Sym,
        'couplingWeight':12.0,
        'couplingFct': lambda couplingSum2Value: couplingSum2Value
    },
    'mutau': {
        'title': scenarioTitle(0,1,1),
        'ylabel': vmu2Sym+" = "+vtau2Sym,
        'couplingWeight':52.0,
        'couplingFct': lambda couplingSum2Value: (math.sqrt(couplingSum2Value)/2.)**2
    }
}

scenarios = {

     'ee': {
        'title': scenarioTitle(1,0,0),
        'ylabel': ve2Sym,
        'couplingWeight':2.0,
        'couplingFct': lambda couplingSum2Value: couplingSum2Value
    },
    'mumu': {
        'title': scenarioTitle(0,1,0),
        'ylabel': vmu2Sym,
        'couplingWeight':12.0,
        'couplingFct': lambda couplingSum2Value: couplingSum2Value
    },
}





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
        
        smoothValues[ipoint] = 0.75*values[ipoint]+0.25*(0.333*(newValues[pointIndices[0]]+newValues[pointIndices[1]]+newValues[pointIndices[2]]))
        
    allPoints = np.concatenate([newPoints,points],axis=0)
    allValues = np.concatenate([newValues,smoothValues],axis=0)
    
    
    if splineFit:
        tck = interpolate.bisplrep(allPoints[:,0], allPoints[:,1], allValues, kx=2, ky=2, s=0.01*allPoints.shape[0], eps=1e-8)
        for ipoint in range(newPoints.shape[0]):
            newValues[ipoint] = 0.75*newValues[ipoint]+0.25*interpolate.bisplev(newPoints[ipoint,0],newPoints[ipoint,1], tck)
        for ipoint in range(points.shape[0]):
            smoothValues[ipoint] = 0.75*smoothValues[ipoint]+0.25*interpolate.bisplev(points[ipoint,0],points[ipoint,1], tck)
            
        
    if addTrigs:
        return np.concatenate([newPoints,points],axis=0),np.concatenate([newValues,smoothValues],axis=0)
    else:
        return points,smoothValues
        

def fitTheoryXsec(massCoupling2Arr,logTheoryXsecs):
    #since cross section scales with |V|^2 one can just fit a 1D function to sigma(mass)/|V|^2 through all points
    normXsecPerMass = {}
    for i in range(massCoupling2Arr.shape[0]):
        mass = int(massCoupling2Arr[i,0]*10)
        coupling2Value = massCoupling2Arr[i,1]
        xsec = math.exp(logTheoryXsecs[i])
        #print mass,coupl,xsec
        if mass not in normXsecPerMass.keys():
            normXsecPerMass[mass] = []
        normXsecPerMass[mass].append(xsec/coupling2Value)
        
    normLogXsecArr = []
    massArr = []
    for mass,normXsec in normXsecPerMass.items():
        normXsec = np.array(normXsec)
        normLogXsecArr.append(np.log(np.mean(normXsec)))
        massArr.append(mass/10.)
        
    massArr = np.array(massArr)
    normLogXsecArr = np.array(normLogXsecArr)
    
    idx = np.argsort(massArr)
    massArr = massArr[idx]
    normLogXsecArr = normLogXsecArr[idx]
        
    #print massArr
    #print normLogXsecArr
        
    fct = scipy.interpolate.interp1d(massArr, normLogXsecArr, kind='cubic')
    
    def getValue(point):
        mass = np.clip(point[0],1.,20.)
        coupling2Value = point[1]
        normLogXsec = fct(mass)
        return math.log(math.exp(normLogXsec)*coupling2Value)

    return getValue

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
    transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=False,splineFit=True)
    transformedPoints,values = smoothPoints(transformedPoints,values,addTrigs=False,splineFit=True)
    
    print (points.shape,transformedPoints.shape)
    
    delaunay = scipy.spatial.Delaunay(transformedPoints)
    
    kneighbors = NearestNeighbors(n_neighbors=24, algorithm='ball_tree').fit(transformedPoints)

    def getValue(point):
        transformedPoint = []
        for i,transformation in enumerate(transformations):
            transformedPoint.append(
                transformation(point[i])
            )
        transformedPoint = np.array(transformedPoint)
        
        distances, pointIndices = kneighbors.kneighbors(np.expand_dims(transformedPoint,axis=0))
        
        pointIndices = pointIndices[0]
        
        '''
        simplexIndex = delaunay.find_simplex(transformedPoint)
        pointIndices = delaunay.simplices[simplexIndex]
        '''
        
        
        simplexPoints = transformedPoints[pointIndices]
        simplexValues = values[pointIndices]
        
        minSimplexPoints = np.amin(simplexPoints,axis=0)
        maxSimplexPoints = np.amax(simplexPoints,axis=0)
        distance = np.fabs(maxSimplexPoints-minSimplexPoints)+1e-6
        
        
        
        #closestIdx = np.argmin(np.sum(np.square(simplexPoints-transformedPoint),axis=1))
        #return simplexValues[closestIdx]
        #print closestIdx
        
        #weight = np.sum(distances[0]/distance,axis=1)
        #weight = np.prod(1.0-np.abs(simplexPoints-transformedPoint)/distance,axis=1)
        
        #weight = 1./(1e-12+distances[0])
        #weight = 1./(0.1*distances[0,-1]+distances[0])
        
        weight = np.square(np.mean(distances[0])/(np.mean(distances[0])+distances[0]))
        #print weight/np.sum(weight)
        
        #print simplexValues.shape,weight.shape,(simplexValues*weight).shape
        
        #return 0.9*np.sum(values[pointIndices]*weight)/np.sum(weight)+0.1*np.mean(values[pointIndices])
        return np.sum(simplexValues*weight)/np.sum(weight)
        '''
        if distances[0,0]>0.05:
            return -3.0
        else:
            return values[pointIndices[0]]
        '''
        
    return getValue

'''
def getContour(graph,level=1.0):
    graph.GetHistogram().SetContour(1,array('d',[1.0]))
    graph.Draw("cont list")
    contourList = []
    contLevel = graph.GetContourList(level)
    if contLevel and contLevel.GetSize()>0:
        # print (graph.GetName())
        # print (contLevel.GetSize())
        for i in range(0,contLevel.GetSize()):
            contourTemp = contLevel[i].Clone()
            if type(contourTemp) != int: contourTemp.SetName(graph.GetName()+"_r1Contour_{0}".format(i))
            else: continue
            contourTemp.SetLineWidth(3)
            contourList.append(contourTemp)
        contourList = sorted(contourList,key = lambda contour: contour.GetN())[::-1]
        return contourList
    else:
        return None
    
def traceContour(logXsecTheoryFct, logXsecLimitFct):
    logMassArr = np.linspace(math.log10(1.),math.log10(20.),20)
    logCoupling2Arr = np.linspace(-7,0,20)
    outputGraph = ROOT.TGraph2D(len(logMassArr)*len(logCoupling2Arr))
    
    for i in range(logMassArr.shape[0]):
        for j in range(logCoupling2Arr.shape[0]):
            logXsecTheory = logXsecTheoryFct([10**logMassArr[i],10**logCoupling2Arr[j]])
            logXsecLimit = logXsecLimitFct([10**logMassArr[i],10**logCoupling2Arr[j]])
        
            outputGraph.SetPoint(
                i*logCoupling2Arr.shape[0]+j,
                logMassArr[i],
                logCoupling2Arr[j],
                math.exp(logXsecTheory)/math.exp(logXsecLimit)
            )
    outputGraph.SetNpx(400)
    outputGraph.SetNpy(400)
    
    interpHist = outputGraph.GetHistogram().Clone(outputGraph.GetName()+"interp")
    contours = getContour(outputGraph,1.)
    
    #if contours and len(contours)>=1:
    print contours
    if contours and len(contours)>=1:
        n = contours[0].GetN()
        graph = ROOT.TGraph(n)
        for i in range(n):
            graph.SetPoint(i,10**contours[0].GetX()[i],10**contours[0].GetY()[i])
        return graph
    return None

'''
def traceContour(logXsecTheoryFct, logXsecLimitFct):
    startingMass = 4.5
    startingCoupling2 = -1.
    
    def fct(logMass,logCoupling2):
        mass = 10**logMass
        coupling2 = 10**logCoupling2
        logXsecTheory = logXsecTheoryFct([mass,coupling2])
        logXsecLimit = logXsecLimitFct([mass,coupling2])
        return logXsecLimit-logXsecTheory
    
    for coupling2Value in np.logspace(-7,1,100):
        logXsecTheory = logXsecTheoryFct([startingMass,coupling2Value])
        logXsecLimit = logXsecLimitFct([startingMass,coupling2Value])
        if logXsecTheory>logXsecLimit:
            startingCoupling2 = coupling2Value
            break
    if startingCoupling2<0.:
        raise Exception("No crossing found")
        
    
    
    
    def findNextPoint(currentPoint, grad, stepSize, previousAngle=0.):
        print "find next",currentPoint,stepSize
        while (True):
            minValue = 1e12
            minAngle = 0.
            for angle in np.linspace(-0.6*math.pi,0.6*math.pi,100):
                value = fct(
                    currentPoint[0]-stepSize*(grad[1]*math.cos(angle)-grad[0]*math.sin(angle)),
                    currentPoint[1]+stepSize*(grad[0]*math.cos(angle)+grad[1]*math.sin(angle))
                )
                if math.fabs(value)<minValue:
                    minValue = math.fabs(value)
                    minAngle = angle
                    
            #beta = 0.5
            #minAngle = beta*minAngle+(1.-beta)*previousAngle
            
            print "\tbest angle",minAngle/math.pi*180.,minValue,stepSize
            
            newPoint = [
                    currentPoint[0]-stepSize*(grad[1]*math.cos(minAngle)-grad[0]*math.sin(minAngle)),
                    currentPoint[1]+stepSize*(grad[0]*math.cos(minAngle)+grad[1]*math.sin(minAngle))
                ]
            if 10**newPoint[0]>20. or 10**newPoint[0]<1.:
                return newPoint,stepSize,minAngle
                
                
            if stepSize>=7e-2 and math.fabs(minAngle)>15./180.*math.pi:
                stepSize = np.clip(stepSize*0.8,2e-2,1e-1)
            elif stepSize<=7e-2 and math.fabs(minAngle)<2./180.*math.pi:
                stepSize = np.clip(stepSize*1.2,2e-2,1e-1)
            else:
                return newPoint,stepSize,minAngle
    
    
    currentLogMass = np.log10(startingMass)
    currentLogCoupling2 = np.log10(startingCoupling2)
    
    
    previousGrad = [0.,0.]
    iteration = 0
    contourMass = []#[10**currentLogMass]
    contourCoupling2 = []#[10**currentLogCoupling2]
    
    currentStepSize = 7e-2
    currentAngle = 0.
    
    
    while (True):
        
        h = 1e-2
        
        
        diffMass = (fct(currentLogMass+h,currentLogCoupling2)-fct(currentLogMass,currentLogCoupling2))/h
        diffCoupling = (fct(currentLogMass,currentLogCoupling2+h)-fct(currentLogMass,currentLogCoupling2))/h
        
        norm = math.sqrt(diffMass**2+diffCoupling**2)
        
        diffMass = diffMass/norm
        diffCoupling = diffCoupling/norm
        
        beta = 0.7
        if iteration>0:
            diffMass = beta*diffMass+(1.-beta)*previousGrad[0]
            diffCoupling = beta*diffCoupling+(1.-beta)*previousGrad[1]
        
        previousGrad = [diffMass,diffCoupling]
        

        nextPoint, currentStepSize, currentAngle = findNextPoint(
            [currentLogMass,currentLogCoupling2],
            [-diffMass,-diffCoupling],
            currentStepSize,
            currentAngle
        )
        
        currentLogMass=0.1*currentLogMass+0.8*nextPoint[0]+0.1*(currentLogMass+currentStepSize*previousGrad[1])
        currentLogCoupling2=0.1*currentLogCoupling2+0.8*nextPoint[1]+0.1*(currentLogCoupling2-currentStepSize*previousGrad[0])
        
        #currentLogMass = currentLogMass+currentStepSize*previousGrad[1]
        #currentLogCoupling2 = currentLogCoupling2-currentStepSize*previousGrad[0]
        
        #currentLogMass=nextPoint[0]
        #currentLogCoupling2=nextPoint[1]
        
        
        contourMass.append(10**currentLogMass)
        contourCoupling2.append(10**currentLogCoupling2)
        
        if (10**currentLogMass<1. or 10**currentLogMass>20.):
            break

        if (iteration>100):
            break
        
        #print 10**currentLogMass,currentLogCoupling2
        
        iteration+=1
        
        print iteration,10**currentLogMass,currentLogCoupling2
    
    
    contourMass = contourMass[::-1]
    contourCoupling2 = contourCoupling2[::-1]
    
    currentLogMass = np.log10(contourMass[-1])
    currentLogCoupling2 = np.log10(contourCoupling2[-1])
    
    previousGrad = [0.,0.]
    iteration = 0
    
    currentStepSize = 7e-2
    currentAngle = 0.
    
    
    while (True):
        
        h = 1e-2
        
        
        diffMass = (fct(currentLogMass+h,currentLogCoupling2)-fct(currentLogMass,currentLogCoupling2))/h
        diffCoupling = (fct(currentLogMass,currentLogCoupling2+h)-fct(currentLogMass,currentLogCoupling2))/h
        
        norm = math.sqrt(diffMass**2+diffCoupling**2)
        
        diffMass = diffMass/norm
        diffCoupling = diffCoupling/norm
        
        beta = 0.7
        if iteration>0:
            diffMass = beta*diffMass+(1.-beta)*previousGrad[0]
            diffCoupling = beta*diffCoupling+(1.-beta)*previousGrad[1]
        
        previousGrad = [diffMass,diffCoupling]
        

        nextPoint, currentStepSize, currentAngle = findNextPoint(
            [currentLogMass,currentLogCoupling2],
            [diffMass,diffCoupling],
            currentStepSize,
            currentAngle
        )
        
        currentLogMass=0.1*currentLogMass+0.8*nextPoint[0]+0.1*(currentLogMass-currentStepSize*previousGrad[1])
        currentLogCoupling2=0.1*currentLogCoupling2+0.8*nextPoint[1]+0.1*(currentLogCoupling2+currentStepSize*previousGrad[0])
        
        #currentLogMass = currentLogMass-currentStepSize*previousGrad[1]
        #currentLogCoupling2 = currentLogCoupling2+currentStepSize*previousGrad[0]
        #currentLogMass=nextPoint[0]
        #currentLogCoupling2=nextPoint[1]
        
        contourMass.append(10**currentLogMass)
        contourCoupling2.append(10**currentLogCoupling2)
        
        if (10**currentLogMass<1. or 10**currentLogMass>20.):
            break
        
        if (iteration>100):
            break
        
        #print 10**currentLogMass,currentLogCoupling2
        
        iteration+=1
        
        print iteration,10**currentLogMass,currentLogCoupling2
        
    return np.clip(np.array(contourMass,dtype=np.float64),1.,20.),np.clip(np.array(contourCoupling2,dtype=np.float64),1e-7,1)
    #print startingCoupling2




with open("/vols/cms/LLP/gridpackLookupTable.json") as lookup_table_file:
    lookup_table = json.load(lookup_table_file)

for year in years:
    for scenarioName,scenarioCfg in scenarios.items():
        for hnl_type in ["majorana","dirac"]:#,"majorana"]:
            print("Plotting coupling scenario: "+str(scenarioName),hnl_type,year)


            masses = []
            coupling2Values = []
            logTheoryXsecs = []
            logLimitObs = []
            logLimitExp = []
            logLimitExp68Up = []
            logLimitExp68Down = []
            logLimitExp95Up = []
            logLimitExp95Down = []
            
            for f in os.listdir(json_path):
                if not f.endswith('.json'):
                    continue
                if f.find(hnl_type)<0:
                    continue
                if f.find(year)<0:
                    continue
                    
                if "_pt20_" not in f and f.replace("_all_", "_pt20_") in os.listdir(json_path):
                    print("Skipping "+f+", higher stat. sample exists")
                    continue
                    
                limits = parse_limit_json(f, scenarioCfg['couplingWeight'])
                if limits is None:
                    continue
                
                    
                mass, couplingSum2Value, theoryXsec = parse_lookup_table(f, lookup_table, scenarioCfg['couplingWeight'])
                coupling2Value = scenarioCfg['couplingFct'](couplingSum2Value)

                
                masses.append(mass)
                coupling2Values.append(coupling2Value)
                logTheoryXsecs.append(math.log(theoryXsec))
                logLimitExp.append(math.log(limits['exp0']))
                logLimitExp68Up.append(math.log(limits['exp+1']))
                logLimitExp68Down.append(math.log(limits['exp-1']))
                logLimitExp95Up.append(math.log(limits['exp+2']))
                logLimitExp95Down.append(math.log(limits['exp-2']))
                
                if not blinded:
                    logLimitObs.append(math.log(limits['obs']))
                
                #expand to 0.5 to remove edge artifacts
                if int(mass)==1:
                    masses.append(0.5)
                    coupling2Values.append(coupling2Value)
                    logTheoryXsecs.append(math.log(theoryXsec))
                    logLimitExp.append(math.log(limits['exp0']))
                    logLimitExp68Up.append(math.log(limits['exp+1']))
                    logLimitExp68Down.append(math.log(limits['exp-1']))
                    logLimitExp95Up.append(math.log(limits['exp+2']))
                    logLimitExp95Down.append(math.log(limits['exp-2']))
                    if not blinded:
                        logLimitObs.append(math.log(limits['obs']))
                
            
            masses = np.array(masses,dtype=np.float32)
            coupling2Values = np.array(coupling2Values,dtype=np.float32)
            massCoupling2Arr = np.stack([masses,coupling2Values],axis=1)
            logTheoryXsecs = np.array(logTheoryXsecs,dtype=np.float32)
            logLimitExp = np.array(logLimitExp,dtype=np.float32)
            logLimitExp68Up = np.array(logLimitExp68Up,dtype=np.float32)
            logLimitExp68Down = np.array(logLimitExp68Down,dtype=np.float32)
            logLimitExp95Up = np.array(logLimitExp95Up,dtype=np.float32)
            logLimitExp95Down = np.array(logLimitExp95Down,dtype=np.float32)
            
            if not blinded:
                logLimitObs = np.array(logLimitObs,dtype=np.float32)
            
            
            logTheoryXsecsFct = fitTheoryXsec(massCoupling2Arr,logTheoryXsecs)
            
            #logTheoryXsecsFct = interpolatedFct(massCoupling2Arr,logTheoryXsecs,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])

            logLimitExpFct = interpolatedFct(massCoupling2Arr,logLimitExp,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
            logLimitExp68UpFct = interpolatedFct(massCoupling2Arr,logLimitExp68Up,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
            logLimitExp68DownFct = interpolatedFct(massCoupling2Arr,logLimitExp68Down,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
            logLimitExp95UpFct = interpolatedFct(massCoupling2Arr,logLimitExp95Up,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
            logLimitExp95DownFct = interpolatedFct(massCoupling2Arr,logLimitExp95Down,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
            
            if not blinded:
                logLimitObsFct = interpolatedFct(massCoupling2Arr,logLimitObs,transformations=[lambda x: np.log10(x)/np.log10(20.), lambda x: np.log10(x)/7.])
                
            rootObj = []
            
            cv = style.makeCanvas("cvLimit"+str(random.random()),700,670)
            cv.SetPad(0.0, 0.0, 1.0, 1.0)
            cv.SetFillStyle(4000)
            cvxmin=0.145
            #cvxmax=0.8
            cvxmax=0.965
            cvymin=0.13
            cvymax=0.92
            cv.SetBorderMode(0)
            cv.SetGridx(False)
            cv.SetGridy(False)
            cv.SetFrameBorderMode(0)
            cv.SetFrameBorderSize(1)
            cv.SetFrameFillColor(0)
            cv.SetFrameFillStyle(0)
            cv.SetFrameLineColor(1)
            cv.SetFrameLineStyle(1)
            cv.SetFrameLineWidth(1)
            cv.SetLeftMargin(cvxmin)
            cv.SetRightMargin(1-cvxmax)
            cv.SetTopMargin(1-cvymax)
            cv.SetBottomMargin(cvxmin)
            cv.SetTitle("")
            cv.SetTickx(1) 
            cv.SetTicky(1)
            
            cv.SetLogx(1)
            cv.SetLogy(1)
            cv.SetLogz(1)
            
            xsec2DHist = ROOT.TH2F(
                "xsec2DHist"+str(random.random()),";"+mHNLSym+" (GeV); "+scenarioCfg['ylabel'],
                200,np.logspace(math.log10(1.),math.log10(20.),201),
                200,np.logspace(-7,0,201),
            )
            xsec2DHist.GetYaxis().SetTitleOffset(1.4)
            xsec2DHist.GetXaxis().SetLabelSize(0)
            xsec2DHist.GetXaxis().SetTickLength(0.0)
            '''
            for ibin in range(xsec2DHist.GetNbinsX()):
                massValue = xsec2DHist.GetXaxis().GetBinCenter(ibin+1)
                for jbin in range(xsec2DHist.GetNbinsY()):
                    coupling2Value = xsec2DHist.GetYaxis().GetBinCenter(jbin+1)
                    
                    if blinded:
                        logLimitValue = logLimitExpFct([massValue,coupling2Value])
                    else:
                        logLimitValue = logLimitObsFct([massValue,coupling2Value])
                    
                    #logLimitValue = logTheoryXsecsFct([massValue,coupling2Value])
                    #xsec2DHist.SetBinContent(ibin+1,jbin+1,math.exp(np.clip(logLimitValue,-10,10)))
                    
                    
                    logTheoryValue = logTheoryXsecsFct([massValue,coupling2Value])
                    xsec2DHist.SetBinContent(ibin+1,jbin+1,math.exp(np.clip(logLimitValue,-10,10))/math.exp(np.clip(logTheoryValue,-10,10)))
            
                    
            #xsec2DHist.GetXaxis().SetMoreLogLabels()
            #xsec2DHist.GetXaxis().SetTickLength(0.015/(1-cv.GetLeftMargin()-cv.GetRightMargin()))
            xsec2DHist.GetXaxis().SetTickLength(0.0)
            xsec2DHist.GetXaxis().SetLabelSize(0)
            
            xsec2DHist.GetYaxis().SetTickLength(0.015/(1-cv.GetTopMargin()-cv.GetBottomMargin()))
            
            if blinded:
                xsec2DHist.GetZaxis().SetTitle("Expected 95% CL upper limit (pb)")
            else:
                xsec2DHist.GetZaxis().SetTitle("Observed 95% CL upper limit (pb)")
            
            #xsec2DHist.GetZaxis().SetTitle("Theory HNL production (pb)")
            
            #xsec2DHist.GetZaxis().SetTitle("Signal strength (#sigma#lower[0.3]{#scale[0.7]{obs.}}/#sigma#lower[0.3]{#scale[0.7]{theo.}})")
            
            xsec2DHist.GetZaxis().SetTitleOffset(1.5)
            xsec2DHist.Draw("colz")
            '''
            xsec2DHist.Draw("AXIS")
            
            '''
            for ipoint in range(massCoupling2Arr.shape[0]):
                marker = ROOT.TMarker(massCoupling2Arr[ipoint,0],massCoupling2Arr[ipoint,1],20)
                marker.SetMarkerColor(ROOT.kWhite)
                marker.SetMarkerSize(1.)
                rootObj.append(marker)
                marker.Draw("Same")
            
                if massCoupling2Arr[ipoint,0]>20. or massCoupling2Arr[ipoint,0]<1.:
                    continue
                if massCoupling2Arr[ipoint,1]>1e0 or massCoupling2Arr[ipoint,1]<1e-7:   
                    continue
                markerText = ROOT.TText(
                    massCoupling2Arr[ipoint,0],
                    math.exp(math.log(massCoupling2Arr[ipoint,1])+0.1),
                    "%.2f"%(logTheoryXsecs[ipoint])
                )
                markerText = ROOT.TText(
                    massCoupling2Arr[ipoint,0],
                    math.exp(math.log(massCoupling2Arr[ipoint,1])+0.1),
                    "%.2f"%(logLimitObs[ipoint])
                )
                markerText = ROOT.TText(
                    massCoupling2Arr[ipoint,0],
                    math.exp(math.log(massCoupling2Arr[ipoint,1])+0.1),
                    "%.2f"%(math.log10(math.exp(logLimitObs[ipoint])/math.exp(logTheoryXsecs[ipoint]))
                ))
                
                
                markerText.SetTextAlign(21)
                markerText.SetTextFont(43)
                markerText.SetTextSize(20)
                rootObj.append(markerText)
                markerText.Draw("Same")
            '''
            
            expMassContourArr,expCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitExpFct)
            exp68UpMassContourArr,exp68UpCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitExp68UpFct)
            exp68DownMassContourArr,exp68DownCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitExp68DownFct)
            exp95UpMassContourArr,exp95UpCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitExp95UpFct)
            exp95DownMassContourArr,exp95DownCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitExp95DownFct)
            
            
            exp95MassContourArr = np.concatenate([exp95UpMassContourArr,np.flip(exp95DownMassContourArr)])
            exp95Coupling2ContourArr = np.concatenate([exp95UpCoupling2ContourArr,np.flip(exp95DownCoupling2ContourArr)])
            
            exp95GraphLimit = ROOT.TGraph(exp95MassContourArr.shape[0],exp95MassContourArr,exp95Coupling2ContourArr)
            rootObj.append(exp95GraphLimit)
            exp95GraphLimit.SetFillStyle(1001)
            exp95GraphLimit.SetLineWidth(0)
            exp95GraphLimit.SetFillColor(limitYellow.GetNumber())
            exp95GraphLimit.Draw("F")
            
            exp68MassContourArr = np.concatenate([exp68UpMassContourArr,np.flip(exp68DownMassContourArr)])
            exp68Coupling2ContourArr = np.concatenate([exp68UpCoupling2ContourArr,np.flip(exp68DownCoupling2ContourArr)])
            
            exp68GraphLimit = ROOT.TGraph(exp68MassContourArr.shape[0],exp68MassContourArr,exp68Coupling2ContourArr)
            rootObj.append(exp68GraphLimit)
            exp68GraphLimit.SetFillStyle(1001)
            exp68GraphLimit.SetLineWidth(0)
            exp68GraphLimit.SetFillColor(limitGreen.GetNumber())
            exp68GraphLimit.Draw("F")
            
            expGraphLimit = ROOT.TGraph(expMassContourArr.shape[0],expMassContourArr,expCoupling2ContourArr)
            rootObj.append(expGraphLimit)
            expGraphLimit.SetLineWidth(2)
            expGraphLimit.SetLineStyle(2)
            expGraphLimit.SetLineColor(ROOT.kBlack)
            expGraphLimit.Draw("L")
            
            
            
            if not blinded:
                obsMassContourArr,obsCoupling2ContourArr = traceContour(logTheoryXsecsFct,logLimitObsFct)
                obsGraphLimit = ROOT.TGraph(obsMassContourArr.shape[0],obsMassContourArr,obsCoupling2ContourArr)
                rootObj.append(obsGraphLimit)
                obsGraphLimit.SetLineWidth(2)
                obsGraphLimit.SetLineStyle(1)
                obsGraphLimit.SetLineColor(ROOT.kBlack)
                obsGraphLimit.Draw("L")
                  
            #print contourList
            
            #if contour:
            #    contour.Draw("L")
            
            
            
            
            pCMS=ROOT.TPaveText(cvxmin+0.035,cvymax-0.025,cvxmin+0.035,cvymax-0.025,"NDC")
            pCMS.SetFillColor(ROOT.kWhite)
            pCMS.SetBorderSize(0)
            pCMS.SetTextFont(63)
            pCMS.SetTextSize(27)
            pCMS.SetTextAlign(13)
            pCMS.AddText("CMS")
            pCMS.Draw("Same")
            
            pPreliminary=ROOT.TPaveText(cvxmin+0.125,cvymax-0.025,cvxmin+0.125,cvymax-0.025,"NDC")
            pPreliminary.SetFillColor(ROOT.kWhite)
            pPreliminary.SetBorderSize(0)
            pPreliminary.SetTextFont(53)
            pPreliminary.SetTextSize(27)
            pPreliminary.SetTextAlign(13)
            pPreliminary.AddText("Preliminary")
            
            
            '''
            couplingTitle = ROOT.TPaveText(cvxmin+0.025,cvymax-0.075,cvxmin+0.025,cvymax-0.075,"NDC")
            couplingTitle.SetBorderSize(0)
            couplingTitle.SetFillStyle(0)
            couplingTitle.SetTextFont(63)
            couplingTitle.SetTextSize(27)
            couplingTitle.SetTextAlign(13)
            couplingTitle.AddText(scenarioCfg['title'])
            couplingTitle.Draw("Same")
            '''
            
            couplingTitle = ROOT.TPaveText(cvxmin,0.94,cvxmin,0.94,"NDC")
            couplingTitle.SetBorderSize(0)
            couplingTitle.SetFillStyle(0)
            couplingTitle.SetTextFont(43)
            couplingTitle.SetTextSize(30)
            couplingTitle.SetTextAlign(11)
            couplingTitle.AddText(scenarioCfg['title'])
            couplingTitle.Draw("Same")

            pLumi=ROOT.TPaveText(cvxmax,0.94,cvxmax,0.94,"NDC")
            pLumi.SetFillColor(ROOT.kWhite)
            pLumi.SetBorderSize(0)
            pLumi.SetTextFont(43)
            pLumi.SetTextSize(30)
            pLumi.SetTextAlign(31)
            pLumi.AddText(lumi[year]+"#kern[-0.5]{ }fb#lower[-0.7]{#scale[0.7]{-1}} (13 TeV)")
            pLumi.Draw("Same")
            
            legendTitle = ROOT.TPaveText(cvxmin+0.035,cvymax-0.02-0.04,cvxmin+0.035,cvymax-0.02-0.04-2*0.037,"NDC")
            legendTitle.SetBorderSize(0)
            legendTitle.SetFillStyle(0)
            legendTitle.SetTextFont(43)
            legendTitle.SetTextSize(20)
            legendTitle.SetTextAlign(13)
            if hnl_type=="dirac":
                ttext = legendTitle.AddText("Dirac HNL production")
                ttext.SetTextFont(63)
            else:
                ttext = legendTitle.AddText("Majorana HNL production")
                ttext.SetTextFont(63)
            legendTitle.AddText("95% CL upper limits")
            legendTitle.Draw("Same")
            
            legend = ROOT.TLegend(cvxmin+0.035,cvymax-0.02-0.055-2*0.037,cvxmin+0.4,cvymax-0.02-0.055-2*0.037-0.037*4,"","NDC")
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(43)
            legend.SetTextSize(20)
            legend.AddEntry(obsGraphLimit,"Observed","L")
            legend.AddEntry(expGraphLimit,"Expected (median)","L")
            legend.AddEntry(exp68GraphLimit,"Expected (68% quantile)","F")
            legend.AddEntry(exp95GraphLimit,"Expected (95% quantile)","F")
            legend.Draw("Same")
            
            refEntries = refLimits.getRefs(hnl_type,int(scenarioCfg['couplingWeight']))
            
            yrefOffset = cvymax-0.03
            for refEntry in refEntries:
                legendRefs = ROOT.TLegend(cvxmax-0.34,yrefOffset,cvxmax-0.05,yrefOffset-0.031*len(refEntry),"","NDC")
                rootObj.append(legendRefs)
                yrefOffset -= 0.031*len(refEntry)+0.010
                legendRefs.SetBorderSize(0)
                legendRefs.SetFillStyle(0)
                legendRefs.SetTextFont(43)
                legendRefs.SetTextSize(16)
            
                for entry in refEntry:
                    
                    if entry[0]:
                        legendRefs.AddEntry(entry[0],entry[1],"L")
                        entry[0].Draw("L")
                    else:
                        legendRefs.AddEntry("",entry[1],"")
            
                legendRefs.Draw("Same")
            
            
            ROOT.gPad.RedrawAxis()
            
            for xtick in [1,2,3,4,6,8,10,14,20]:
                tickLine = ROOT.TLine(xtick,1.5e-7,xtick,1e-7)
                rootObj.append(tickLine)
                tickLine.Draw("Same")
                
                tickLine = ROOT.TLine(xtick,1.,xtick,7e-1)
                rootObj.append(tickLine)
                tickLine.Draw("Same")
                
                tickText = ROOT.TText(xtick,4e-8,str(xtick))
                tickText.SetTextAlign(21)
                tickText.SetTextFont(43)
                tickText.SetTextSize(29)
                rootObj.append(tickText)
                tickText.Draw("Same")
                
            for xtick in [5,7,9,11,12,13,15,16,17,18,19]:
                tickLine = ROOT.TLine(xtick,1.2e-7,xtick,1e-7)
                rootObj.append(tickLine)
                tickLine.Draw("Same")
                
                tickLine = ROOT.TLine(xtick,1.,xtick,8.2e-1)
                rootObj.append(tickLine)
                tickLine.Draw("Same")
            
            #cv.Update()
            cv.Print(scenarioName+"_"+hnl_type+"_"+year+"_xsec.pdf")
            cv.Print(scenarioName+"_"+hnl_type+"_"+year+"_xsec.png")
            pPreliminary.Draw("Same")
            cv.Print(scenarioName+"_"+hnl_type+"_"+year+"_xsec_PAS.pdf")
            cv.Print(scenarioName+"_"+hnl_type+"_"+year+"_xsec_PAS.png") 
                
            '''
            cv = style.makeCanvas("cvXsec"+str(random.random()),700,670)
            cv.SetPad(0.0, 0.0, 1.0, 1.0)
            cv.SetFillStyle(4000)
            cvxmin=0.14
            #cvxmax=0.96
            cvxmax=0.81
            cvymin=0.13
            cvymax=0.92
            cv.SetBorderMode(0)
            cv.SetGridx(False)
            cv.SetGridy(False)
            cv.SetFrameBorderMode(0)
            cv.SetFrameBorderSize(1)
            cv.SetFrameFillColor(0)
            cv.SetFrameFillStyle(0)
            cv.SetFrameLineColor(1)
            cv.SetFrameLineStyle(1)
            cv.SetFrameLineWidth(1)
            cv.SetLeftMargin(cvxmin)
            cv.SetRightMargin(1-cvxmax)
            cv.SetTopMargin(1-cvymax)
            cv.SetBottomMargin(cvxmin)
            cv.SetTitle("")
            cv.SetTickx(1) 
            cv.SetTicky(1)
            
            cv.SetLogx(1)
            cv.SetLogy(1)
            cv.SetLogz(1)
            
            xsec2DHist = ROOT.TH2F(
                "xsec2DHist"+str(random.random()),";"+mHNLSym+" (GeV); "+scenarioCfg['ylabel'],
                200,np.logspace(math.log10(1.),math.log10(20.),201),
                200,np.logspace(-7,0,201),
            )
            for ibin in range(xsec2DHist.GetNbinsX()):
                massValue = xsec2DHist.GetXaxis().GetBinCenter(ibin+1)
                for jbin in range(xsec2DHist.GetNbinsY()):
                    couplingValue = xsec2DHist.GetYaxis().GetBinCenter(jbin+1)
                    if blinded:
                        logLimitValue = logLimitExpFct([massValue,couplingValue])
                    else:
                        logLimitValue = logLimitObsFct([massValue,couplingValue])
                    xsec2DHist.SetBinContent(ibin+1,jbin+1,math.exp(logLimitValue))
                    
            xsec2DHist.GetXaxis().SetMoreLogLabels()
            xsec2DHist.GetXaxis().SetTickLength(0.015/(1-cv.GetLeftMargin()-cv.GetRightMargin()))
            xsec2DHist.GetYaxis().SetTickLength(0.015/(1-cv.GetTopMargin()-cv.GetBottomMargin()))
            if blinded:
                xsec2DHist.GetZaxis().SetTitle("Expected 95% CL upper limit (pb)")
            else:
                xsec2DHist.GetZaxis().SetTitle("Expected 95% CL upper limit (pb)")
            xsec2DHist.Draw("colz")
            
            
            pCMS=ROOT.TPaveText(cvxmin+0.025,cvymax-0.025,cvxmin+0.025,cvymax-0.025,"NDC")
            pCMS.SetFillColor(ROOT.kWhite)
            pCMS.SetBorderSize(0)
            pCMS.SetTextFont(63)
            pCMS.SetTextSize(27)
            pCMS.SetTextAlign(13)
            pCMS.AddText("CMS")
            pCMS.Draw("Same")

            pLumi=ROOT.TPaveText(cvxmax,0.94,cvxmax,0.94,"NDC")
            pLumi.SetFillColor(ROOT.kWhite)
            pLumi.SetBorderSize(0)
            pLumi.SetTextFont(43)
            pLumi.SetTextSize(30)
            pLumi.SetTextAlign(31)
            pLumi.AddText(lumi[year]+"#kern[-0.5]{ }fb#lower[-0.7]{#scale[0.7]{-1}} (13 TeV)")
            pLumi.Draw("Same")
            
            cv.Print(scenarioName+"_"+hnl_type+"_"+year+"_xsec.pdf")
            '''
            
            
            
            
            
            
