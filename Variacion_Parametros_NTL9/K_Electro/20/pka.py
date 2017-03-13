import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

f = open('result.txt','r')
primerpaso=0
linea=[]
residuos=[] # Residuos cargados de la simulacion
results=[]
phs=[] #pHs de la simulacion
#Filtro las lineas del archivo resultado

for line in f:
	if primerpaso==0:
		linea=line.split()
		residuos=[int(i)+1 for i in linea[1:]]
		#print(residuos)
		primerpaso=1
	else:
		linea=line.split()
		phs.append(float(linea[0]))
		results.append([float(i) for i in linea[1:]])
		

#Guardo los resultados de cada linea en orden para cada residuo

deprots=[] # Fracciones deprotonadas de la simulacion
for i in range(len(residuos)):
	fracdeprot=[]
	for line in results:
		fracdeprot.append(line[i])
	deprots.append(fracdeprot)

#print(deprots[0])
#print(len(deprots[0]))
#print(phs)
	


#Henderson-Hasselbalch
def acidos(x, a, b):
	y = 1/(1+10**(b*(x-a)))
	return y
def bases(x, a, b):
	y = 1/(1+10**(-b*(x-a)))
	return y


x = np.linspace(phs[0], phs[len(phs)-1], 100)


for i in range(len(residuos)):
	deprot=deprots[i]
	if deprot[len(deprot)-1]<0.5: # Me fijo si es acido o basico para el fit
		#Es acido
		try:		
			popt, pcov = curve_fit(acidos, phs, deprot,p0=[7,1])
			plt.plot(phs,deprot,'*')
			plt.plot(x,acidos(x,popt[0],popt[1]))
		except RuntimeError:		
			plt.plot(phs,deprot,'*')
	else:
		#Es base
		try:
			popt, pcov = curve_fit(bases, phs, deprot,p0=[7,1])
			plt.plot(phs,deprot,'*')
			plt.plot(x,bases(x,popt[0],popt[1]))
		except RuntimeError:		
			plt.plot(phs,deprot,'*')
	#print(str(residuos[i])+" "+ str(popt[0]))
	print(str(popt[0]))

plt.show()





