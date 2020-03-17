class Box1D:
    """ iniciar L en metros y V en meVs. Funciones principales: encontrar_energias()->arreglo con las energias
    para el estado n con energia E ->funcion_normalizada(x,n,E) da la funcion de ese estado"""
    hbar=1.05457e-34 #J*s
    qe=1.6e-19#
    m=9.10938291e-31#kg
    meff=0.067*m #masa efectiva GaAs
    

    
    def __init__(self,L,V_0):
        self.L=L
        self.V_0=self.qe*V_0/1000
        
        energias=self.encontrar_energias()
        constantes_normalizacion=[]
        for i in range(len(energias)):
            constantes_normalizacion.append(self.consts_norm(i+1))
        self.constantes_normalizacion=constantes_normalizacion
    
    def f_par(self,E):
        import numpy as np
        k=np.sqrt(2*self.meff*E)/self.hbar
        ka=np.sqrt(2*self.meff*(self.V_0-E))/self.hbar
        return k*np.tan(k*self.L/2)-ka
    
    
    def f_impar(self,E):
        import numpy as np
        k=np.sqrt(2*self.meff*E)/self.hbar
        ka=np.sqrt(2*self.meff*(self.V_0-E))/self.hbar
        return k*( 1 / (np.tan(k*self.L/2) )) + ka
    
    
    def encontrar_energias(self):
        """encuentra las energias y devuelve un arreglo"""
        import numpy as np
        from scipy.optimize import bisect
        energias = np.linspace(0,250e-3,10000)#rango de energias a recorrer en el ciclo
        energia_ant=0.0  #varible para comparar el cambio de signo
        energias_pares=[] 
        energias_impares=[]
        i=0
        for energia_now in energias:
            energia_now=energia_now*self.qe # paso a joules
            if( self.f_par(energia_now)*self.f_par(energia_ant) < 0) and i%2==0: #cambio de signo en pares
                energias_pares.append( 1000*bisect(self.f_par,energia_ant,energia_now,xtol=1e-20)/self.qe ) #guardo en mevs
                i+=1
            if ( self.f_impar(energia_now)*self.f_impar(energia_ant) < 0 )and i%2!=0: #cambio de signo impares
                energias_impares.append( 1000*bisect(self.f_impar,energia_ant,energia_now,xtol=1e-20)/self.qe ) #guardo en mevs
                i+=1
            energia_ant=energia_now
        energias_todas=[*energias_pares,*energias_impares]
        energias_todas.sort()
        return energias_todas
    
    
    def funcion_onda(self,x,n,E):
        """ devuelve la funcion de onda n (empezando en 1) dada la energia que debe
        ser encontrada antes E con la funcion encontrar_energias() """
        import numpy as np
        k=np.sqrt(2*self.meff*E)/self.hbar
        ka=np.sqrt(2*self.meff*(self.V_0-E))/self.hbar
        cte_par=np.cos(k*self.L/2)/np.exp(-ka*self.L/2)
        cte_impar=np.sin(k*self.L/2)/np.exp(-ka*self.L/2)
        if n%2!=0: #n impar pero funcion es par (n empieza en 1)#
            if -self.L/2<=x and x<=self.L/2:
                return np.cos(k*x)
            elif -self.L/2>x:
                return cte_par*np.exp(ka*x)
            elif x>self.L/2:
                return cte_par*np.exp(-ka*x)
        elif n%2==0:
            if -self.L/2<=x and x<=self.L/2:
                return np.sin(k*x)
            elif -self.L/2>x:
                return -cte_impar*np.exp(ka*x) 
            elif x>self.L/2:
                return cte_impar*np.exp(-ka*x)

    
    def modulo_cuadrado(self,x,n,E):
        return  self.funcion_onda(x,n,E)**2


    def funcion_normalizada(self,x,n,E):
        import numpy as np
        return self.funcion_onda(x,n,E)/np.sqrt(self.constantes_normalizacion[n-1])


    def consts_norm(self,n):
        import scipy.integrate as integrate
        const_norm=[]
        energias=self.encontrar_energias()
        i=0
        for en in energias:
            i+=1
            const_norm.append( integrate.quad(self.modulo_cuadrado,-2*self.L,2*self.L,args=(i,en*self.qe/1000))[0] )
        return const_norm[n-1]
    

    def densidad_prob(self,x,n,E):
        return self.funcion_normalizada(x,n,E)**2


    def potencial(self,x):
        if -self.L/2<=x and x<=self.L/2:
            return 0.0
        elif -self.L/2>x:
            return self.V_0*(1000/self.qe)
        elif x>self.L/2:
            return self.V_0*(1000/self.qe)
    def constante_ya_normalizada(self,n):
        import scipy.integrate as integrate
        const_norm=[]
        energias=self.encontrar_energias()
        i=0
        for en in energias:
            i+=1
            const_norm.append( integrate.quad(lambda x,n,E: self.funcion_normalizada(x,n,E)**2,-2*self.L,2*self.L,args=(i,en*self.qe/1000))[0] )
        return const_norm[n-1]
