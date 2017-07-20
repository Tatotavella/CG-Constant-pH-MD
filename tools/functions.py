import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def autocorr(y,tau):
    mean = np.mean(y)
    var = np.var(y)
    partsum = 0.0
    for t in range(len(y)-tau):
        partsum += (y[t] - mean)*(y[t+tau] - mean)
    corr = partsum/(len(y)-tau)
    return corr/var

def deproto(stateDirec,meanInit,resList,decorrStep):
	g = open(stateDirec,'r')


	outputDepr = {}
	for r in resList:
		outputDepr.update({r:[]}) #Deprotonated fraction

	resMeans = []
	for q in resList:
		resMeans.append(0) #Initialization for auxiliary mean calculator

	k = 0
	totalStates = 0
	next(g) #Jump Header
	for line in g:
		if(k > meanInit):
			if((k-meanInit)%decorrStep == 0): #To get decorrelated states
				inter = [x for x in line.split('\t')]
				for idx,q in enumerate(inter[1:len(inter)-1]):
					resMeans[idx]+=abs(float(q)) #Times it was charged
				totalStates += 1
		k+=1
	g.close()
	print("Data points: "+str(totalStates))
	for idx in range(len(resMeans)):
		resMeans[idx] = resMeans[idx]/totalStates
	for idx in range(len(resList)):
		outputDepr[resList[idx]].append(resMeans[idx])

	return outputDepr


def qTot(stateDirec,meanInit,decorrStep):
	g = open(stateDirec,'r')

	chargeMean = 0
	csqrtMean = 0

	k = 0
	totalStates = 0
	next(g) #Jump Header
	for line in g:
		if(k > meanInit):
			if((k-meanInit)%decorrStep == 0): #To get decorrelated states
				inter = [x for x in line.split('\t')]
				q = inter[len(inter)-1]
				chargeMean += float(q)
				csqrtMean += float(q)*float(q)
				totalStates += 1
		k+=1
	g.close()
	chargeMean = chargeMean/totalStates
	csqrtMean = csqrtMean/totalStates
	chargeVar = csqrtMean - chargeMean**2

	return [chargeMean,chargeVar]

def eTot(logDirec,meanInit,decorrStep):
	g = open(logDirec,'r')

	energyMean = 0
	esqrtMean = 0

	k = 0
	totalStates = 0
	next(g) #Jump Header
	for line in g:
		if(k > meanInit):
			if((k-meanInit)%decorrStep == 0): #To get decorrelated states
				inter = [x for x in line.split('\t')]
				e = inter[len(inter)-1]
				energyMean += float(e)
				esqrtMean += float(e)*float(e)
				totalStates += 1
		k+=1
	g.close()
	energyMean = energyMean/totalStates
	esqrtMean = esqrtMean/totalStates
	energyVar = esqrtMean - energyMean**2

	return [energyMean,energyVar]



def acidos(x, a, b):
	y = 1.0/(1.0+10.0**(-1.0*b*(x-a)))
	return y
def bases(x, a, b):
	y = 1.0/(1.0+10.0**(b*(x-a)))
	return y

def fitHenHass(phs,deprot,let):
	pka = -1
	x = np.linspace(phs[0], phs[len(phs)-1], 100)
	if let=='D' or let=='E' or let=='C' or let=='B' or let=='Y': #Acid
		try:
			popt, pcov = curve_fit(acidos, phs, deprot,p0=[5,1])
			plt.plot(phs,deprot,'*')
			plt.plot(x,acidos(x,popt[0],popt[1]))
			pka = popt[0]
		except RuntimeError:
			a=1
			plt.plot(phs,deprot,'*')
	elif let=='K' or let=='R' or let=='N' or let=='H': #Base
		try:
			popt, pcov = curve_fit(bases, phs, deprot,p0=[10,1])
			plt.plot(phs,deprot,'*')
			plt.plot(x,bases(x,popt[0],popt[1]))
			pka = popt[0]
		except RuntimeError:
			a=1
			plt.plot(phs,deprot,'*')
	return pka
