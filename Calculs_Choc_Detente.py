import numpy as np

angle_de_deflexion = 8

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

Choc = Choc_Oblique(angle_de_deflexion * np.pi / 180, 2, 288.15, 101325, 1.225044, 1.4)

print("\nChoc Oblique:\n")      
print('Angle de choc = '+str(Choc[0]*180/np.pi)+"°")
print('M_1 = '+str(Choc[1]))
print('T_1 = '+str(Choc[2])+" K")
print('p_1 = '+str(Choc[3])+" Pa")
print('rho_1 = '+str(Choc[4])+" kg/m^3")

gamma = 1.4
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

Detente = Detente(16 * np.pi / 180, Choc[1], Choc[2], Choc[3], Choc[4], 1.4)

print("\nDétente de Prandtl-Meyer:\n")  
print('M_2 = '+str(Detente[0]))
print('T_2 = '+str(Detente[1])+" K")
print('p_2 = '+str(Detente[2])+" Pa")
print('rho_2 = '+str(Detente[3])+" kg/m^3")