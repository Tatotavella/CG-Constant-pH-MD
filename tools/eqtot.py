import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from functions import *


totalSteps = 10000
meanInit = 500
decorrStep = 30

dataPoints = (totalSteps - meanInit)/decorrStep

print('Data points: '+str(dataPoints))


pHs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14']
phsVal = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
charges = []
chargesVar = []

for ph in pHs:
	#MC.state file
	data_dir = os.path.dirname(__file__) #<-- Absolute directory
	rel_path = "../"+ph+"/MC_result/MC.state"
	stateDirec = os.path.join(data_dir, rel_path)
	q = qTot(stateDirec,meanInit,decorrStep)
	charges.append(q[0])
	chargesVar.append(q[1])

print('q states Analyzed')


energies = []
energiesVar = []

for ph in pHs:
	#MC.log file
	data_dir = os.path.dirname(__file__) #<-- Absolute directory
	rel_path = "../"+ph+"/MC_result/MC.log"
	logDirec = os.path.join(data_dir, rel_path)
	e = eTot(logDirec,meanInit,decorrStep)
	energies.append(e[0])
	energiesVar.append(e[1])

print('e states Analyzed')


f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
#ax1.plot(phsVal,charges,'bo')
ax1.errorbar(phsVal,charges,yerr=chargesVar,fmt='bo')
ax1.plot(phsVal,charges,'b')
ax1.plot([0,15],[0,0],'g')
ax1.set_xlim([1,14])
ax1.set_xlabel('pH')
ax1.set_ylabel('Protein net charge')

#ax2.plot(phsVal,energies,'ro')
ax2.errorbar(phsVal,energies,energiesVar,fmt='ro')
ax2.plot(phsVal,energies,'r')
ax2.plot([0,15],[0,0],'g')
ax2.set_xlim([1,14])
ax2.set_xlabel('pH')
ax2.set_ylabel('Mean Electrostatic Energy')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0.2)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.show()


f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
ax1.plot(phsVal,chargesVar,'bo')
ax1.plot(phsVal,chargesVar,'b')
ax1.plot([0,15],[0,0],'g')
ax1.set_xlim([1,14])
ax1.set_xlabel('pH')
ax1.set_ylabel('Protein net charge variance')

ax2.plot(phsVal,energiesVar,'ro')
ax2.plot(phsVal,energiesVar,'r')
ax2.plot([0,15],[0,0],'g')
ax2.set_xlim([1,14])
ax2.set_xlabel('pH')
ax2.set_ylabel('Mean Electrostatic Energy variance')

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0.2)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.show()




