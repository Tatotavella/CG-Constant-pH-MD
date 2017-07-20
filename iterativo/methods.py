import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import sys
from Bio.PDB import *


def debhuck(r,q1,q2,kpp,kmm,kpm,ld,eps):
	if q1*q2<=0:
		y = (kpm*q1*q2*np.exp(-r/ld))/(eps*r)
		return y
	elif q1>0 and q2>0:
		y = (kpp*q1*q2*np.exp(-r/ld))/(eps*r)
		return y
	elif q1<0 and q2<0:
		y = (kmm*q1*q2*np.exp(-r/ld))/(eps*r)
		return y
	else:
		print("ERROR debhuck")

def qpart(aob,dpka):
	if aob=='a':
		qp = -1.0/(1.0+10**(-dpka))
		return qp
	elif aob=='b':
		qp = +1.0/(1.0+10**(+dpka))
		return qp
	else:
		print("ERROR qpart")

def coarse_beta(pdb):
    if(not isinstance(pdb, str)):
        return 0
    else:
        p = PDBParser()
        structure = p.get_structure(pdb.upper(), pdb+".pdb")
        # Get residue positions from a structure
        res = []
        bx = []
        by = []
        bz = []
        for chains in structure:
            for chain in chains:
                for residue in chain:
                    atom_list = Selection.unfold_entities(residue, 'A')
                    res_list = Selection.unfold_entities(residue, 'R')
                    name = residue.get_resname()
                    try:
                        beta = residue['CB']
                        r = beta.get_vector()
                        bx.append(r[0])
                        by.append(r[1])
                        bz.append(r[2])
                        res.append(name)
                    except KeyError:
                        p = alpha_carbon(residue)
                        if p==0:
                            pass
                        else:
                            res.append(name)
                            bx.append(p[0])
                            by.append(p[1])
                            bz.append(p[2])

		return [[bx,by,bz],res]

def alpha_carbon(residue):
    try:
        alpha = residue['CA']
        p = alpha.get_vector()
        return p
    except KeyError:
        return 0


def p_np_count(i,chargedList,res,chrgPos,allPos,wrsh_params,polars,nonpolars):
    x = chrgPos[0]
    y = chrgPos[1]
    z = chrgPos[2]
    bx = allPos[0]
    by = allPos[1]
    bz = allPos[2]
    alpha_pol = wrsh_params[0]
    alpha_non_pol = wrsh_params[1]
    r_pol = wrsh_params[2]
    r_nonpol = wrsh_params[3]
    resi = chargedList[i]
    Npol=0.0
    Nnonpol=0.0
    for j,r in enumerate(res):
    	if j!=resi:
    		#Distancia
    		dx=x[i]-bx[j]
    		dy=y[i]-by[j]
    		dz=z[i]-bz[j]
    		dis=np.sqrt(dx**2+dy**2+dz**2)
    		#Non Polar
    		if r in nonpolars:
    			if dis<=r_nonpol:
    				Nnonpol+=1
    			else:
    				Nnonpol+=np.exp(-alpha_non_pol*(dis-r_nonpol)*(dis-r_nonpol))
    		elif r in polars:
    			#Polar
    			if dis<=r_pol:
    				Npol+=1
    			else:
    				Npol+=np.exp(-alpha_pol*(dis-r_pol)*(dis-r_pol))
    return [Npol,Nnonpol]

def p_np_penalty(i,typelist,Npol,Nnonpol,pty_params,wrsh):
    alpha_u_pol = pty_params[0]
    alpha_u_nonpol = pty_params[1]
    N_p_max = pty_params[2]
    N_np_max = pty_params[3]
    pol_pty=0.0
    np_pty=0.0
    let = typelist[i]
    params = wrsh[let]
    self_pol = params[0]
    self_nonpol = params[1]
    if Npol<N_p_max:
    	pol_pty=self_pol*np.exp(-alpha_u_pol*(Npol-N_p_max)*(Npol-N_p_max))
    else:
    	pol_pty=self_pol
    if Nnonpol<N_np_max:
    	np_pty=self_nonpol*np.exp(-alpha_u_nonpol*(Nnonpol-N_np_max)*(Nnonpol-N_np_max))
    else:
    	np_pty=self_nonpol

    return [pol_pty,np_pty]

def elec_penalty(i,chargedList,typelist,ph,pka_old,acidbase,chrgPos,elec_params):
    x = chrgPos[0]
    y = chrgPos[1]
    z = chrgPos[2]
    kpp = elec_params[0]
    kmm = elec_params[1]
    kpm = elec_params[2]
    ld = elec_params[3]
    eps = elec_params[4]
    resi = chargedList[i]
    e_elec=0.0
    for j,resj in enumerate(chargedList):
        if resj!=resi:
            typei = typelist[i]
            typej = typelist[j]
            dpkai=ph-pka_old[i]
            dpkaj=ph-pka_old[j]
            aobi=acidbase[typei]
            aobj=acidbase[typej]
            # Charge j
            qj=qpart(aobj,dpkaj)
            # Distance
            dx=x[i]-x[j]
            dy=y[i]-y[j]
            dz=z[i]-z[j]
            r = np.sqrt(dx**2+dy**2+dz**2)
            e_elec+=debhuck(r,1.0,qj,kpp,kmm,kpm,ld,eps)
    return e_elec
