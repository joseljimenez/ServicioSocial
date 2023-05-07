#importando las librerias
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#Definiendo las constantes.
#Constantes fisicas.  
omega_DE=0.7
omega_m=1-omega_DE
h=.7
yr=1
Mpc=3.26e6*yr 
m=Mpc/(3.086e22)
eV=1/(0.197e-6*m)
M_p=2.4e27*eV
H_0=(100/3e5)*h/Mpc
rho_c0=3*H_0**2*M_p**2 
G=1/(8*np.pi*M_p**2)    #M_p**2=1/(8*pi*G)
t_uni=13.7e9*yr
t_eq=10e11*3e8*m
#Cosntantes para el Modelo G_3(X,phi)=mu*X**2.
k=100  #cte del lagrangiano
mu=5

#Condiciones inciales:
a_ini=1/3e3   #a_eq materia radiacion
X_ini=.001
t_ini=t_eq    #tiempo de incio

#EDO para a
def eq_for_a(y,a):
    return np.sqrt(k/(3*M_p**2))*y*a
#EDO para X
def eq_for_phi(y,a):
    return -3*np.sqrt(k/(3*M_p**2))*y**2

#Sistema de ecuaciones
def system(V,t):
    y,a=V
    eval_1=eq_for_phi(y,a)
    eval_2=eq_for_a(y,a)
    return [eval_1,eval_2]
    
#condiciones inciales
i_conditions=[X_ini,a_ini]
#print(i_conditions)

#tiempo 
#ime_array=np.linspace(t_ini,t_uni/50)
time_array=np.linspace(t_ini,t_ini,1000)

#resolvemos
sol=odeint(system,i_conditions,time_array)

#
numeric_X=sol[:,0]
numeric_a=sol[:,1]










