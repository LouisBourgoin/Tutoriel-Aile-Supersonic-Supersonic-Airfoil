import numpy as np

#paramètres d'entrée
Incidence = 4
Incidence_rad = Incidence * np.pi/180
DAA = 8 #demi-angle aïgu de l'aile losangique
angle_de_deflexion_1 = DAA-Incidence
angle_de_deflexion_2 = DAA+Incidence

M_0 = 2
T_0 = 288.15
p_0 = 101325
rho_0 = 1.225044
gamma = 1.4


def Choc_Oblique(theta, M_0, T_0, p, rho_0, gamma):
  
    
    tan_theta = np.tan(theta)    
    
    for Beta in np.arange(1,500) * np.pi/1000:
        x = 2 / np.tan(Beta) * (M_0**2 * np.sin(Beta)**2 - 1) / (M_0**2 * (gamma + np.cos(2 * Beta)) + 2)
        if x > tan_theta:
            break
    
    M_1 = 1 / np.sin(Beta - theta) * np.sqrt((1 + (gamma - 1)/2 * M_0**2 * np.sin(Beta)**2) / (gamma * M_0**2 * np.sin(Beta)**2 - (gamma - 1)/2))
    h = M_0 ** 2 * np.sin(Beta) ** 2
    T_1 = T_0 * (2 * gamma * h - (gamma - 1)) * ((gamma - 1) * h + 2) / ((gamma + 1) ** 2 * h)
    p_1 = p * (2 * gamma * h - (gamma - 1)) / (gamma + 1)
    rho_1 = rho_0 * ((gamma + 1) * h) / ((gamma - 1) * h + 2)

    return Beta, M_1, T_1, p_1, rho_1

Choc1 = Choc_Oblique(angle_de_deflexion_1 * np.pi / 180, M_0, T_0, p_0, rho_0, gamma)
Choc2 = Choc_Oblique(angle_de_deflexion_2 * np.pi / 180, M_0, T_0, p_0, rho_0, gamma)


A = np.sqrt((gamma+1)/(gamma-1)) #Fonction de Prandtl-Meyer
B = (gamma-1)/(gamma+1)
PM = lambda x:  A * np.arctan(np.sqrt(B * (x**2 - 1))) - np.arctan(np.sqrt(x**2 - 1)) 

def Detente(theta, M_1, T_1, p_1, rho_1, gamma):
    
    nu_M2 = theta + PM(M_1)
    for M_2 in np.arange(1,10,0.001):
        x = PM(M_2)
        if x > nu_M2:
             break
    T_2 = T_1 * (1 + (gamma-1)/2 * M_1**2)/(1 + (gamma-1)/2 * M_2**2)
    p_2 = p_1 * (T_2/T_1)**(gamma/(gamma-1))
    rho_2 = rho_1 * (T_2/T_1)**(1/(gamma-1))
    
    
    return M_2,T_2, p_2, rho_2

Detente1 = Detente(16 * np.pi / 180, Choc1[1], Choc1[2], Choc1[3], Choc1[4], 1.4)
Detente2 = Detente(16 * np.pi / 180, Choc2[1], Choc2[2], Choc2[3], Choc2[4], 1.4)


Fn = ((Choc2[3] + Detente2[2]) - (Choc1[3] + Detente1[2]))/2
Fp = 0.5 * np.tan(DAA*np.pi/180) * ( (Choc1[3] + Choc2[3]) - (Detente1[2] + Detente2[2]) )


Fz = Fn*np.cos(Incidence_rad) - Fp*np.sin(Incidence_rad)
Fx = Fn*np.sin(Incidence_rad) + Fp*np.cos(Incidence_rad)

a = np.sqrt(gamma*287*T_0)
V = M_0*a

Cz = Fz / (0.5 * rho_0 * 1 * V**2)
Cx = Fx /(0.5 * rho_0 * 1 * V**2)

print('\nFz théorique = '+str(Fz))
print('Fx théorique = '+str(Fx)+"\n")
print('Cz théorique = '+str(Cz))
print('Cx théorique = '+str(Cx)+"\n")
