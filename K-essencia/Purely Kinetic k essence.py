#!/usr/bin/env python
# coding: utf-8

# In[3]:


from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt


# In[4]:


#Constantes fisicas. Modelo de Scherer F(X)=F0+F2(X-X0)**2. F tiene un minimo en X0
omega_DE=0.7
omega_m=1-omega_DE

Mpc=1
h=.7
m=Mpc/(3.086e22)
eV=1/(0.197e-6*m)
M_p=2.4e27*eV
H_0=(100/3e5)*h/Mpc
rho_c0=3*H_0**2*M_p**2
G=1/(8*np.pi*M_p**2)    #M_p**2=1/(8*pi*G)
X_0=1
a_ini=1/3e3  #a_eq materia radiacion
e1_a1=0.0001
y_ini=np.sqrt(2*(1+ e1_a1/a_ini**3 )*X_0)  #y_ini=np.sqrt(2*X_ini)
X_ini=y_ini**2/2
#print(e1_a1/a_ini**3)
F_0=-omega_DE*rho_c0       #
F_2=omega_m*rho_c0/(4*X_0**2*e1_a1)
yr=9.467e15*m
t_uni=13.7e9*yr
t_eq=10e11*3e8*m
print("t_uni = ", t_uni)
print("F_0 = ", F_0)
print("F_2 = ", F_2)


# In[5]:


yr


# In[6]:


#Solucion analitica para X=X(a)
def X_analit(a):
    X = X_0*(1+e1_a1/a**3)
    return X


# In[7]:


#Arreglo para el factor de escala desde la epoca de igualdad materia-radiacion hasta hoy en dia.
a_arreglo=np.linspace(a_ini,1,1000) #a_ini=a_eq


# In[8]:


#ploteo para X=X(a) analitico desde la epoca de igualdad hasta hoy en dia.
plt.plot(a_arreglo,X_analit(a_arreglo))
plt.xscale("log")
plt.yscale("log")
plt.xlabel("a")
plt.ylabel("X(a)")
plt.title("Solucion analitica para X=X(a) desde la epoca de iguldad hasta hoy")


# In[9]:


#Solución analítica para t=t(a)
#A
def ta(a):
    a_eq=a_ini
    A=-8*np.pi*G*F_0/3
    B=32*np.pi*G*F_2*(X_0**2)*e1_a1/3
    
    l1=2*np.log( np.sqrt(A)*a**(3/2) + np.sqrt(A*a**3+B))
    l2=3*np.sqrt(A)
    c=-2*np.log(B)/(3*np.sqrt(A))
    #print(A*a**3,B)
    return l1/l2 + c    


# In[10]:


#ta(1)


# In[11]:


#Algo interesante es que ta(a_arreglo) da valores negativos
plt.plot(ta(a_arreglo) ,a_arreglo,'.', label='Solucion analitica para a=a(t) considerando DM & DE')   

#plt.xscale("log")
plt.yscale("log")
plt.ylabel("a(t)")
plt.xlabel("t")
plt.title("Solucion analitica para a=a(t) considerando DM & DE")
plt.legend()
plt.show()


# In[40]:


t_analitica=ta(a_arreglo)
plt.plot(t_analitica[:-1] ,(a_arreglo[1:]-a_arreglo[:-1])/(t_analitica[1:]-t_analitica[:-1]), label='Delta(a)/Delta(t) considerando DM & DE')   

#plt.xscale("log")
#plt.yscale("log")
plt.ylabel("a(t)")
plt.xlabel("t")
plt.title("Solucion analitica para a=a(t) considerando DM & DE")
plt.legend()
plt.show()


# In[13]:


#Aqui tambien ta(a_arreglo) da valores negativos
B=np.sqrt(4*F_2*X_0**2*e1_a1/(3*M_p**2))
plt.plot(2*a_arreglo**(3/2)/(3*B), a_arreglo, ".", label="a=a(t) analitico para la parte de DM")

plt.plot(np.log(a_arreglo)/(np.sqrt(-F_0/(3*M_p**2))), a_arreglo, ".", label="a=a(t) anali.para la parte de DE")

#plt.xscale("log")
#plt.yscale("log")
plt.ylabel("a(t)")
plt.xlabel("t")
#plt.xlim(-10,10000)
plt.title("Solucion analitica para a=a(t)")
plt.legend()
plt.show()


# In[14]:


plt.plot(np.log(a_arreglo)/(np.sqrt(-F_0/(3*M_p**2))), a_arreglo, ".", label="a=a(t) para la parte de DE")
plt.plot(t_analitica ,a_arreglo,'.', label='Solucion analitica para a=a(t) considerando DM & DE')
#plt.xscale("log")
#plt.yscale("log")
plt.ylabel("a(t)")
plt.xlabel("t")
#plt.xlim(-10,10000)
plt.title("Solucion analitica para a=a(t)")
plt.legend()
plt.show()


#Algo curiosos es que salen tiempos negativos


# In[46]:


#Aqui tambien ta(a_arreglo) da valores negativos
B=np.sqrt(4*F_2*X_0**2*e1_a1/(3*M_p**2))
plt.plot((2*a_arreglo**(3/2)/(3*k1))/yr,'.' , a_arreglo, label="a=a(t) para la parte de Materia")
plt.plot((ta(a_arreglo)-t_analitica[0])/yr ,a_arreglo, '.', label='Solucion analitica para a=a(t) considerando DM & DE')

plt.xscale("log")
plt.yscale("log")
plt.ylabel("a(t)")
plt.xlabel("t")
#plt.xlim(-10,10000)
plt.title("Solucion analitica para a=a(t)")
plt.legend()
plt.show()


# In[16]:


#Para el termino de materia

k1=np.sqrt(4*F_2*X_0**2*e1_a1/(3*M_p**2))

xxx=np.log(np.abs(ta(a_arreglo)))
yyy=np.log(a_arreglo)
deltaxxx=xxx[1:]-xxx[:-1]
deltayyy=yyy[1:]-yyy[:-1]
plt.plot(xxx[1:],deltayyy/deltaxxx,".", label='pendiente para el termino de materia')
plt.legend()
plt.ylim(-5,5)


# In[17]:


#Para el termino de energia
#
xxx=np.abs(ta(a_arreglo))
yyy=np.log(a_arreglo)
deltaxxx=xxx[1:]-xxx[:-1]
deltayyy=yyy[1:]-yyy[:-1]
plt.plot(xxx[1:],deltayyy/deltaxxx,".", label='pendiente para el termino de energia')
plt.ylim(-.5,.5)
plt.legend()


# In[18]:


#Solucion analítica de y=y(a), y=phi_punto desde la epoca de igualdad hasta hoy en dia.
def y_analit(a):
    return np.sqrt(2*X_analit(a))


# In[19]:


#Ploteo de y=y(a), y=phi_punto desde la epoca de igualdad hasta hoy en dia.
plt.plot(a_arreglo,y_analit(a_arreglo), '.', label='solucion y(a) analitica')
plt.xlabel("a")
plt.ylabel("\phi(a)")
plt.xscale("log")
plt.yscale("log")
plt.title("Solucion analitica para y(a)=\phi_punto(a) desde la epoca de igualdad hasta hoy en dia")
plt.legend()


# In[20]:


#Grafica \phi=phi(t), ta() no tiene valores positivos
plt.plot(ta(a_arreglo),y_analit(a_arreglo),'.', label='Solucion analitica para y=y(t)')
plt.xlabel("t")
plt.ylabel("y(t)")
#plt.xscale("log")
plt.yscale("log")
plt.title("Solucion analitica para y(t)=\phi_punto(t) desde la epoca de igualdad hasta hoy en dia")
plt.legend()


# In[21]:


#Pediente para y(t)
#
xxx=np.abs(ta(a_arreglo))
yyy=np.log(y_analit(a_arreglo))
deltaxxx=xxx[1:]-xxx[:-1]
deltayyy=yyy[1:]-yyy[:-1]
plt.plot(xxx[1:],deltayyy/deltaxxx,".", label='pendiete para y(t) analitico')
plt.ylim(-.5,.5)


# In[22]:


#La siguiente solucion numerica considera la evolucion de los parametros desde la epoca de 
#igualdad hasta hoy en dia


# In[23]:


#Solucion numerica para y=\phi_punto(t)
def ec_for_y(y,a):
    l=8*np.pi*G/3
    m=(y**2/2 -X_0)*y
    n=(3*y**2)/2-X_0
    A=m/n
    return -3*np.sqrt(l*(F_2*(3*y**4/4 + y**2*X_0 - X_0**2)-F_0))*A


# In[24]:


#Solucion numerica para a(t)
def ec_for_a(y,a):
    return np.sqrt((8*np.pi*G/3)*(F_2*(3*y**4/4 + y**2*X_0 - X_0**2)-F_0))*a


# In[25]:


#Sistema de ecuaciones
def funcion_vectorial(x,t):
    y,a=x
    eval_1=ec_for_y(y,a)
    eval_2=ec_for_a(y,a)
    return [eval_1,eval_2]


# In[26]:


condiciones_iniciales=[y_ini,a_ini]  #El factor de escala va desde la epoca de igualdad hasta hoy en dia


# In[27]:


tn=np.linspace(t_eq,t_uni,100000)    #El tiempo va desde la epoca de igualdad hasta hoy en dia


# In[28]:


#tlog=np.logspace(np.log10(t_eq),np.log10(t_uni),10000)   #tiempo logaritmico
#tlog


# In[29]:


sol=odeint(funcion_vectorial,condiciones_iniciales,tn)


# In[30]:


#sol


# In[31]:


#ploteo para \phi=\phi(t)
plt.plot(tn,sol[:,0], label='Solucion \phi=\phi(t) numerica')
plt.plot(ta(a_arreglo)-t_analitica[0],y_analit(a_arreglo), label='Solucion para y=y(t) analitica ')

plt.xlabel('time')
plt.ylabel('y(t)')
plt.xscale("log")
plt.yscale("log")
plt.title("Solución numerica para y=y(t)")
plt.legend()
#print(sol[:,0])


# In[32]:


"""
#Pediente para y(t) numerico
#
XX=tn
YY=np.log(y_analit(a_arreglo))
deltaXX=XX[1:]-XX[:-1]
deltaYY=YY[1:]-YY[:-1]
plt.plot(XX[1:],deltaYY/deltaXX,".", label='pendiete para y(t) numerico')
"""


# In[33]:


#ploteo para para a=a(t)
plt.plot(tn,sol[:,1], '.', label='Solucion numerica para a=a(t)')
plt.plot(ta(a_arreglo),a_arreglo, '.', label='Solucion analitica para a=a(t)')
plt.title("Solución numérica para a=a(t)")
#plt.ylim(0,100)
#plt.xscale("log") 
#plt.yscale("log")
plt.xlabel('t')
plt.ylabel('a(t)')
plt.legend()


# In[34]:


plt.plot(sol[:,1],sol[:,0],"-.",label="Solucion Numerica para y=y(a)")
plt.plot(a_arreglo,y_analit(a_arreglo),label="Analitica para y=y(a)")

plt.xlabel('a')
plt.ylabel('y(a)')
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.title("Solución numerica para y=y(a)")


# In[35]:


anum=sol[:,1]
ynum=sol[:,0]


# In[36]:


dloga_dt=(np.log(anum[1:])-np.log(anum[:-1]))/(tn[1:]-tn[:-1])


# In[37]:


plt.plot(anum[1:],dloga_dt,".",label='pendiente en ')
plt.xscale("log") 
plt.yscale("log")


# In[38]:


dloga_dlogt=(np.log(anum[1:])-np.log(anum[:-1]))/(np.log(tn[1:])-np.log(tn[:-1]))


# In[39]:


plt.plot(anum[1:],dloga_dlogt,".")
plt.xlim(.01,10)
plt.xscale("log") 
plt.yscale("log")


# In[ ]:




