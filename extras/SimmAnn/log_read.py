# -⁻- coding: UTF-8 -*-
#######################################################################
# 
#  Inputs:
#  - 
#  Outputs:
#  - 
#######################################################################
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from functions import *
import sys

[step, qw, rg, qo, tc, ep] = np.genfromtxt("wolynes_thingsph7sa1", unpack = True)
[step2, qw2, rg2, qo2, tc2, ep2] = np.genfromtxt("wolynes_thingsph7sa2", unpack = True)
[step3, qw3, rg3, qo3, tc3, ep3] = np.genfromtxt("wolynes_thingsph7sa3", unpack = True)
#Temperature
TINI = 300
TFIN = 800
temp = np.zeros(len(step))
for idx in range(len(temp)):
    temp[idx] = ((TFIN - TINI)*step[idx])/step[-1] + TINI


# Graphs
plt.figure(1)
plt.plot(temp,qw,label = 'ph = 7 , SA1')
plt.plot(temp,qw2,label = 'ph = 7, SA2')
plt.plot(temp,qw3,label = 'ph = 7, SA3')
plt.xlabel('Temperature')
plt.ylabel('Q Wolynes')
plt.figure(2)
plt.plot(temp,rg,label = 'ph = 7, SA1')
plt.plot(temp,rg2,label = 'ph = 7, SA2')
plt.plot(temp,rg3,label = 'ph = 7, SA3')
plt.xlabel('Temperature')
plt.ylabel('Radius of gyration')
plt.figure(3)
plt.plot(temp,qo, label = 'ph = 7, SA1')
plt.plot(temp, qo2, label = 'ph = 7, SA2')
plt.plot(temp,qo3,label = 'ph = 7, SA3')
plt.xlabel('Temperature')
plt.ylabel('Q Onuchic')
plt.figure(4)
plt.plot(temp,tc, label = 'ph = 7, SA1')
plt.plot(temp,tc2, label = 'ph = 7, SA2')
plt.plot(temp,tc3,label = 'ph = 7, SA3')
plt.xlabel('Temperature')
plt.ylabel('Total contacts')
plt.figure(5)
plt.plot(temp,ep, label = 'ph = 7, SA1')
plt.plot(temp, ep2, label = 'ph = 7, SA2')
plt.plot(temp,ep3,label = 'ph = 7, SA3')
plt.xlabel('Temperature')
plt.ylabel('Potential energy')
plt.legend()
plt.show()
