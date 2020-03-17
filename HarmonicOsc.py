class HarmonicOsc1D:  
    """ingresar @omegaHbar en meVs.Funciones principales: Energias(n)->la energia del estado n
    para el estado n  ->funcion_de_onda(x,n) da la funcion de ese estado"""
    qe=1.6e-19#carga electron
    
    def __init__(self,omegaHbar):
        """ingresar @omegaHbar en meVs"""
        self.omegaHbar=omegaHbar*self.qe/1000
    
    def Energias(self,n): 
        return 1000*((n+1/2)*self.omegaHbar)/self.qe
    def densidad_prob(self,x,n):
        return self.funcion_de_onda(x,n)**2
    def funcion_de_onda(self,x,n):
        import numpy as np 
        from numpy.polynomial.hermite import hermval 
        from math import factorial
        m=9.10938291e-31#kg
        meff=0.067*m #masa efectiva GaAs
        hbar=1.05457e-34 #J*s
        x0=np.sqrt((hbar**2)/(meff*self.omegaHbar))
        c=np.zeros(n+1)
        c[n]=np.exp(-(x**2)/(2*(x0**2)))/(np.sqrt(np.sqrt(np.pi)*x0*factorial(n)*2**n))
        return hermval(x/x0,c)
    def constante_norm(self,n):
        import scipy.integrate as integrate
        import numpy as np
        m=9.10938291e-31#kg
        meff=0.067*m #masa efectiva GaAs
        hbar=1.05457e-34 #J*s
        a=-10*np.sqrt((hbar**2)/(meff*self.omegaHbar))
        b=-a
        cons_norma=integrate.quad(self.densidad_prob,a,b,args=n)[0]
        return cons_norma
    
    def potencial_arm_1D(self,x):
        import numpy as np
        m=9.10938291e-31#kg
        meff=0.067*m #masa efectiva GaAs
        hbar=1.05457e-34 #J*s
        return 1000*(meff*((self.omegaHbar/hbar)**2)*(x**2)/2)/self.qe #mw2x2/2
        

        

    