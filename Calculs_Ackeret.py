import numpy as np


M = 2
l = 1
e = 2 * (0.5 * np.tan ((8 * np.pi) / 180))
B = (e/l)**2


Fz = 47000.7556339949
Fx = 16543.372
alpha_deg = 4

Cz = Fz / (0.5 * 1.225044 * 1 * 680.3349**2)
alpha = alpha_deg * np.pi/180

Cxp = Cz*alpha + (4*B)/np.sqrt(M**2-1)

print ("\n Coefficient de traînée de pression théorique: "+str(Cxp))

CxpS =   Fx /(0.5 * 1.225044 * 1 * 680.3349**2)

print ("\n Coefficient de traînée de pression simulé: "+str(CxpS))

print("\n Erreur relative: "+ str(abs(((Cxp-CxpS)/Cxp)*100)) + "%")
