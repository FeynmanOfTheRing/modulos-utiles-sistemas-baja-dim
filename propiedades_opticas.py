from sympy.interactive import printing
printing.init_printing(use_latex=True)
from sympy import pprint
from sympy import pprint_try_use_unicode
from sympy import pprint_use_unicode
import sympy as sp 
import numpy as np
from IPython.display import display


""" Parametros y cts """
qe=1.6e-19
m=9.10938291e-31#kg
meff=0.067*m #masa efectiva GaAs
hb = 1.054571e-34
c = 299792458
mu = 4*np.pi*1e-7
N=3e24#m-3
Gamma12 = 1/(0.14e-12)
nr = 3.2
ep = 8.854e-12
er = (nr**2)*ep


def absorption(hw,In,E21,M21,M22,M11):
    """"esta funcion devuelve el valor la absorcion total, lineal y no lineal, 
    para un valor de hw reciviendo como parametros adicionales las integrales M21,M22,M11
    la diferencia de energias E21 y la intensidad de la lus In"""
    alpha1=sp.Function('alpha1')
    alpha3=sp.Function('alpha3')
    alpha=sp.Function('alpha')

    E,I=sp.symbols('E I')
    
    alpha1 = sp.sqrt(mu/er)*N*E*Gamma12*abs(M21)**2/((E21-E)**2+(hb*Gamma12)**2)
    alpha3 = -(sp.sqrt(mu/er)*(I/(2*nr*ep*c))*N*E*Gamma12*abs(M21)**2/((E21-E)**2+(hb*Gamma12)**2)**2)*(4*abs(M21)**2-((abs(M22-M11))**2*(3*E21**2 - 4*E21*E+E**2-(hb*Gamma12)**2))/(E21**2+(hb*Gamma12)**2))
    alpha =  alpha1 + alpha3
    #display(alpha)

    return alpha1.subs(E,hw),alpha3.subs([(E,hw),(I,In)]),alpha.subs([(E,hw),(I,In)])

def mostrar_formulas_abs():
    """esta funcion muestra por consola de manera simbolica la formula de absorcion 
    para verificar que no tenga errores (de tener errores habria que corregir este modulo)"""
    alpha1=sp.Function('alpha1')
    alpha3=sp.Function('alpha3')
    alpha=sp.Function('alpha')

    hbw,I=sp.symbols('hbw I')
    E21,M21,M22,M11,mu,er,N,Gamma12,c,hb,nr,ep=sp.symbols('E21 M21 M22 M11 mu epsilon_r N Gamma12 c hb nr epsilon')
    alpha1 = sp.sqrt(mu/er)*N*hbw*Gamma12*sp.Abs(M21)**2/((E21-hbw)**2+(hb*Gamma12)**2)
    alpha3 = -(sp.sqrt(mu/er)*(I/(2*nr*ep*c))*N*hbw*Gamma12*sp.Abs(M21)**2/((E21-hbw)**2+(hb*Gamma12)**2)**2)*(4*sp.Abs(M21)**2-((sp.Abs(M22-M11))**2*(3*E21**2 - 4*E21*hbw+hbw**2-(hb*Gamma12)**2))/(E21**2+(hb*Gamma12)**2))
    alpha =  alpha1 + alpha3
    pprint_try_use_unicode()
    pprint(alpha1)
    pprint(alpha3)
    pprint(alpha)

def mostrar_parametros():
    """ esta funcion imprime los parametros usados por consola """
    print("carga e ",qe)
    print("masa ",m)
    print("meff ",meff)
    print("hbarra ",hb)
    print("c ",c)
    print("meff ",mu)
    print("N ",N)
    print("$\Gamma_{12}$ ",Gamma12)
    print("nr ",nr)
    print("ep ",ep)
    print("er ",er)

def cambiar_parametros(m_eff,Nportadores,Gamma_12,n_r,e_p,e_r):
    """esta funcion permite modificar los parametros principales en el calculo de las
    propiedades opticas"""
    global meff,N,Gamma12,nr,ep,er
    print('cambio de meff', meff, 'a', m_eff)
    meff=m_eff
    print('cambio de N', N, 'a', Nportadores)
    N=Nportadores
    print('cambio de gamma12', Gamma12 , 'a', Gamma_12)
    Gamma12=Gamma_12
    print('cambio de nr',nr, 'a', n_r)
    nr=n_r
    print('cambio de ep',ep, 'a', e_p)
    ep=e_p
    print('cambio de er',er, 'a', e_r)
    er=e_r











