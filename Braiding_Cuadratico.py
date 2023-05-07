#
#Modelo del lagrangiano lineal G(X)=k*phi**2
#
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
k=1e101  #cte del lagrangiano


#Condiciones inciales:
a_ini=1/3e3   #a_eq materia radiacion
y_ini=.0001
phi_ini=.0055
t_ini=t_eq    #tiempo de incio

#EDO para a
def eq_for_a(y,a,phi):
    return np.sqrt(2*k*phi/(3*M_p**2))*y*a
#EDO para X
def eq_for_y(y,a,phi):
    return -3*np.sqrt(2*k*phi/(3*M_p**2))*y**2 + y**2/(2*phi)
#EDO para y
def eq_for_phi(y,a,phi):
    return y


#Sistema de ecuaciones
def system(V,t):
    y,a,phi=V
    eval_1=eq_for_y(y,a,phi)
    eval_2=eq_for_a(y,a,phi)
    eval_3=eq_for_phi(y,a,phi)
    return [eval_1,eval_2,eval_3]
    
#condiciones inciales
i_conditions=[y_ini,a_ini,phi_ini]
#print(i_conditions)

#tiempo 
#ime_array=np.linspace(t_ini,t_uni/50)
time_array=np.linspace(t_ini,t_uni,10000)

#resolvemos
sol=odeint(system,i_conditions,time_array)

#
numeric_y=sol[:,0]
numeric_a=sol[:,1]
numeric_phi=sol[:,2]

#print(numeric_y)
#print(numeric_a)


fig1 = plt.figure()
plt.plot(numeric_a, numeric_y, '.', label='Numeric solution phi=phi(a)')
plt.xlabel('a')
plt.ylabel('phi(a)')
plt.xscale('log')
plt.yscale('log')
plt.legend()


fig2 = plt.figure()
#Graphic for a=a(t)
plt.plot(time_array, numeric_a, '.', label='Numeric solution a=a(t)')
plt.xlabel('t')
plt.ylabel('a(t)')
#plt.xscale('log')
#plt.yscale('log')
plt.legend()

fig3 = plt.figure()
#Graphic for X=X(t)
plt.plot(time_array, numeric_y, '.', label='Numeric solution y=y(t)')
plt.xlabel('t')
plt.ylabel('phi(t)')
plt.xscale('log')
plt.yscale('log')
plt.legend()

#Evolucion de la densidad de energia
numeric_rho=k*numeric_y**2 

print("numeric_y:", numeric_y)
print(numeric_rho)

fig4 = plt.figure()
#Graphic for phi=phi(t)
plt.plot(time_array[1:], (numeric_rho)[1:], '.', label='Numeric solution rho=rho(t)')
plt.xlabel('t')
plt.ylabel('rho(t)')
plt.xscale('log')#plt.yscale('log')
plt.legend()



#Evolucion de la densidad de energia
numeric_p=k*numeric_y**2


fig5 = plt.figure()
#Graphic for phi=phi(t)
plt.plot(time_array, (numeric_p), '.', label='Numeric solution p=p(t)')
plt.xlabel('t')
plt.ylabel('p(t)')
plt.xscale('log')
plt.yscale('log')
plt.legend()


omega=numeric_p/numeric_rho
fig6 = plt.figure()
#Graphic for phi=phi(t)
plt.plot(time_array, omega, '.', label='Numeric solution w=w(t)')
plt.xlabel('t')
plt.ylabel('w(t)')
#plt.ylim(1.001,-.999)
#plt.xscale('log')
#plt.yscale('log')
plt.legend()



