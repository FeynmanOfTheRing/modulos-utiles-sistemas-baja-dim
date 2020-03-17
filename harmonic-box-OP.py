from Box import Box1D
from HarmonicOsc import HarmonicOsc1D
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fmin
from propiedades_opticas import *

"""  parametros y ctes  """
hbar=1.05457e-34 #J*s
qe=1.6e-19#
m=9.10938291e-31#kg
meff=0.067*m #masa efectiva GaAs

""" comprobacion con harrison y haciendo sentido fisico primero """
# omega_comp=np.sqrt(2*(0.250/10)*qe/(meff*(100e-10)**2))
omega_comp=(5.04e-21)#omega hbar que me garantiza el pico teorico sea igual
print(omega_comp*1000/qe)
### grafico el oscilador armónico y muestro que las funciones son normalizadas
Osc1=HarmonicOsc1D(omega_comp*1000/qe) ### creo objeto oscilador armónico con un hw=5meVs
Fig1 , potOsc = plt.subplots()
potOsc.set_title('funciones Oscilador $\hbar\omega = 5 meV$')

x_plot_osc=np.linspace(-400e-10,400e-10,300)
y_pot_osc=[]
for x_ in x_plot_osc:
    y_pot_osc.append(Osc1.potencial_arm_1D(x_))
potOsc.plot(x_plot_osc,y_pot_osc,'r')
potOsc.set_ylabel('V(x)[meVs]')
potOsc.set_xlabel('m')

funOsc=potOsc.twinx()
cte_escala_osc_graf=Osc1.funcion_de_onda(0,0)/2.2

for i in range(3):
    y_funcion=[]
    for x_ in x_plot_osc:
        y_funcion.append(Osc1.funcion_de_onda(x_,i)/cte_escala_osc_graf  + Osc1.Energias(i))
    funOsc.plot(x_plot_osc,y_funcion)
    print('Constante normalizacion Osc1D estado ',i,'es: ',Osc1.constante_norm(i))
funOsc.legend( ('n=0','n=1','n=2'), loc='right')
funOsc.set_ylabel('$\Psi(x)$')
plt.show()


### Grafico caja y muestro que sus autofunciones son normalizadas
Box1=Box1D(200e-10,250)#instancio el objeto caja

Fig1, potBox= plt.subplots()
potBox.set_title('funciones pozo cuadrado $L=200\AA$ & $V=100 meVs$')

x_plot_box=np.linspace(-300e-10,300e-10,300)
y_pot_box=[]
for x_ in x_plot_box:
    y_pot_box.append(Box1.potencial(x_))
potBox.plot(x_plot_box,y_pot_box)
potBox.set_ylabel('V(x)[meVs]')
potBox.set_xlabel('m')

funBox=potBox.twinx()
energias_caja1=Box1.encontrar_energias() 
cte_escala_box_graf=Box1.funcion_normalizada(0.0,1,energias_caja1[0]*qe/1000)/13

for i in range(3):
    y_funcion=[]
    for x_ in x_plot_box:
        y_funcion.append((Box1.funcion_normalizada(x_,i+1,energias_caja1[i]*qe/1000))/cte_escala_box_graf + energias_caja1[i])
    funBox.plot(x_plot_box,y_funcion)
    print('Constante de normalizacion para el estado ',i+1,'es: ',Box1.constante_ya_normalizada(i+1))
funBox.legend( ('n=1','n=2','n=3') , loc='right')
funBox.set_ylabel('$\Psi(x)$')

plt.show()

mostrar_formulas_abs()
mostrar_parametros()
""" calculo propiedades opticas Osc: """

E21_Osc=(Osc1.Energias(1)-Osc1.Energias(0))*qe/1000
M11=integrate.quad(lambda x: qe*Osc1.funcion_de_onda(x,0)*x*Osc1.funcion_de_onda(x,0),-400e-10,400e-10)[0]
print(M11)
M21=integrate.quad(lambda x: qe*Osc1.funcion_de_onda(x,1)*x*Osc1.funcion_de_onda(x,0),-400e-10,400e-10)[0]
print(M21)
M22=integrate.quad(lambda x: qe*Osc1.funcion_de_onda(x,1)*x*Osc1.funcion_de_onda(x,1),-400e-10,400e-10)[0]
print(M22)

En=np.linspace(0,0.1*qe,500,dtype=float)
absOsc=[]
for E in En:
    absOsc.append(absorption(E,2.5e8,E21_Osc,M21,M22,M11))
plt.plot(En/qe*1000,absOsc)
plt.title('absorption harmonic oscillator')
plt.ylabel('$\\alpha$ ( $\hbar$ $\omega$,I)')
plt.xlabel('photon energy in $meVs$')

""" calculo propiedades opticas Box: """
E1=energias_caja1[0]*qe/1000
E2=energias_caja1[1]*qe/1000
E21_box=E2-E1

print(E21_box)
M11=integrate.quad(lambda x: qe*Box1.funcion_normalizada(x,1,E1)*x*Box1.funcion_normalizada(x,1,E1),-300e-10,300e-10)[0]
print(M11)
M21=integrate.quad(lambda x: qe*Box1.funcion_normalizada(x,2,E2)*x*Box1.funcion_normalizada(x,1,E1),-300e-10,300e-10)[0]
print(M21)
M22=integrate.quad(lambda x: qe*Box1.funcion_normalizada(x,2,E2)*x*Box1.funcion_normalizada(x,2,E2),-300e-10,300e-10)[0]
print(M22)

En=np.linspace(0,0.1*qe,500,dtype=float)
absBox=[]
for E in En:
    absBox.append(absorption(E,2.5e8,E21_box,M21,M22,M11))
plt.plot(En/qe*1000,absBox)
plt.title('absoption infinite quantum well $200\AA$  vs harmonic oscillator with $\hbar\omega = 31.50 mevs$ ')
plt.ylabel('$\\alpha$ ( $\hbar$ $\omega$,I)')
plt.xlabel('photon energy in $meVs$')
plt.legend(('osc total','osc1','osc3','box total','box1','box3'))

plt.show()







