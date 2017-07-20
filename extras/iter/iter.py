import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import sys
from Bio.PDB import *

from methods import *

# Instantaneous pka calculation and optimization for parameters
#
#

#-----------------------------------------------------------------------------------------------------------------
# Info dictionaries

# Data in format [Bpol,Bnonpol]
wrsh = {'ASP':[-0.5,2.5],'GLU':[-0.6,2.5],'LYS':[-0.5,2.0],'ARG':[-0.6,2.2],'HIS':[-0.3,4.0]}

# Model pkas
pkamod = {'ASP': 4.0,'GLU': 4.5,'LYS': 10.6,'ARG': 12.0,'HIS': 6.3}
# Acid base classification
acidbase = {'ASP': 'a','GLU': 'a','LYS': 'b','ARG': 'b','HIS': 'b'}
# Acid list
acids = ['ASP','GLU']
# Bases list
bases = ['LYS','ARG','HIS']
# Polar list
polars = ['CYS','ASP','GLU','HIS','LYS','ASN','GLN','ARG','SER','THR','TYR']
# Non Polar list
nonpolars = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP']

#-----------------------------------------------------------------------------------------------------------------

# Delta self parameters
r_pol = 5.0
r_nonpol = 7.0
alpha_pol = 0.1
alpha_non_pol = 0.1
N_p_max = 6.0
N_np_max = 15.0
alpha_u_pol = 0.1
alpha_u_nonpol = 0.02
# Delta elec parameters
kpp = 4.15
kmm = 4.15
kpm = 4.15
ld = 10.0
eps = 1.0

#-----------------------------------------------------------------------------------------------------------------

#  Experimental pKa Res : pKa for 2lzt
pkaexp = {'7': 2.6,'15': 5.5, '18': 2.8, '35': 6.1, '48': 1.4, '52': 3.6, '66': 1.2, '87': 2.2, '101': 4.5, '119': 3.5}

#  Experimental pKa Res : pKa for 1w4h
#pkaexp = {'4':3.9,'16':4.5,'17':6.5,'20':3.7,'36':3.7,'37':3.2,'39':4.5,'41':5.4}


pdb = '2lzt'
data = coarse_beta(pdb)

betas = data[0]
res = data[1]

print(res)

bx = betas[0]
by = betas[1]
bz = betas[2]
# Relevant residues type, indexes and pkas.
chargedList=[]
typelist=[]
pka_old=[]
# Charged residues positions
x=[]
y=[]
z=[]
acid=0
base=0
total=0
for idx,r in enumerate(res):
	if r in acids:
		chargedList.append(idx)
		typelist.append(r)
		x.append(bx[idx])
		y.append(by[idx])
		z.append(bz[idx])
		pk = pkamod[r]
		pka_old.append(pk)
		acid+=1
	elif r in bases:
		chargedList.append(idx)
		typelist.append(r)
		x.append(bx[idx])
		y.append(by[idx])
		z.append(bz[idx])
		pk = pkamod[r]
		pka_old.append(pk)
		base+=1
total=acid+base

# Beta coarse grained plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(bx,by,bz)
plt.show()

# Charged residues positions plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z)
plt.show()


print(chargedList)

# Delta Self penalty
penalty=[]
chrgPos = [x,y,z]
allPos = [bx,by,bz]
wrsh_params = [alpha_pol,alpha_non_pol,r_pol,r_nonpol]
pty_params = [alpha_u_pol,alpha_u_nonpol,N_p_max,N_np_max]

for i,resi in enumerate(chargedList):
	[Npol,Nnonpol] = p_np_count(i,chargedList,res,chrgPos,allPos,wrsh_params,polars,nonpolars)
	[pol_pty,np_pty] = p_np_penalty(i,typelist,Npol,Nnonpol,pty_params,wrsh)
	penalty.append(pol_pty+np_pty)

# Instantaneous pKa calculation


elec_params = [kpp,kmm,kpm,ld,eps]

R = 0.0019872041
T = 300.0
phs = np.linspace(0.01,20,100)


new_pka = [[] for i in range(len(chargedList))]
electros = [[] for i in range(len(chargedList))]

for ph in phs:
	for i,resi in enumerate(chargedList):

		e_elec = elec_penalty(i,chargedList,typelist,ph,pka_old,acidbase,chrgPos,elec_params)
		e_elec = -1 * e_elec/(R*T*np.log(10))
		e_self = penalty[i]/(R*T*np.log(10))
		e_ph = ph - pka_old[i]
		new_pkai = pka_old[i] + e_elec + e_self

		new_pka[i].append(new_pkai)
		electros[i].append(e_elec)

plt.plot(chargedList,new_pka,'o')
plt.xlabel('Residue')
plt.ylabel('pKa')
plt.title('Instantaneous pKa for different ph')
plt.show()

for jdx,resu in enumerate(chargedList):
	let = typelist[jdx]
	aob = acidbase[let]
	compound_pka = pkamod[let]

	elec = electros[jdx]
	new = new_pka[jdx]

	plt.plot(phs,new,'o',label='Res: '+str(int(resu)+1))
	plt.xlabel('pH')
	plt.ylabel('Instantaneous pKa')
	plt.legend()
	plt.show()

	mod = []
	depr = []
	for idx,ph in enumerate(phs):
		dpka_depr = ph - new[idx]
		q_depr = qpart(aob,dpka_depr)
		dpka_mod = ph - compound_pka
		q_mod = qpart(aob,dpka_mod)
		depr.append(q_depr)
		mod.append(q_mod)


	plt.plot(phs,mod,'bo',label='Compund Model, Res: '+str(int(resu)+1))
	plt.plot(phs,depr,'ro',label='Shift Model, Res: '+str(int(resu)+1))
	plt.xlabel('pH')
	plt.ylabel('Charge')
	try:
		exp_pka = pkaexp[str(resu+1)]
		experimental = []
		for idx,ph in enumerate(phs):
			dpka_exp = ph - exp_pka
			q_exp = qpart(aob,dpka_exp)
			experimental.append(q_exp)
		plt.plot(phs,experimental,'go',label='Experimental, Res: '+str(int(resu)+1))
		# Get pka
	except KeyError:
		pass
	plt.legend()
	plt.show()



'''
#Deprotonated fraction with changing pka_new
depr = []
for idx,ph in enumerate(phs):
	dpka=ph-new[idx]
	aob='a'
	q = qpart(aob,dpka)
	depr.append(q)
plt.plot(phs,depr,'ro',label='multiple pka')
plt.plot(phs,chrg,'bo',label='simple pka')
plt.legend()
plt.show()
plt.figure(1)
plt.plot(phs,dph,'o')
plt.show()
'''

# Beta coarse grained plot
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(bx,by,bz)
plt.show()
'''
# Charged residues positions plot
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z)
plt.show()
'''
# Output
'''
npo = npnnp[i][0]
nnp = npnnp[i][1]
sep = "\t"
f.write('{:<4d}'.format(int(resi)+1)+sep+ \
		'{:<4s}'.format(typelist[i])+sep+ \
		'{:.5f}'.format(ph)+sep + \
		'{:.5f}'.format(exp)+sep + \
		'{:.5f}'.format(pka_old[i])+sep + \
		'{:.5f}'.format(e_elec)+sep + \
		'{:.5f}'.format(e_self)+sep + \
		'{:.5f}'.format(npo)+sep + \
		'{:.5f}'.format(nnp)+"\n")

print(int(resi)+1,typelist[i],ph,exp,pka_old[i],e_elec,e_self,npo ,nnp)
'''
