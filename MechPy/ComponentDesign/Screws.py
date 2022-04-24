# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 16:52:42 2022

@author: valte
"""
import numpy as np


alpha_table = {'sqr' : 0,
               'acme' : 14.5,
               'iso' : 15}

thread_list = [(['square','quadrado','quadrada','sqr',0],'sqr'),
               (['acme',1],'acme'),
               (['trapezoidal','trapezoid','trapezoide','trapezio','trapézio','iso',2],'iso')]

def get_thread(keyword):
    for comb in thread_list:
        err = True
        if keyword in comb[0]:
            err = False
            return comb[1]
        
    if err:
        raise ValueError(f'Perfil de rosca invalido.\nPerfis de rosca devem ser: square; acme ou iso.\nEm vez disso, foi recebido {keyword}')

class PowerScrew():
    def __init__(self,
                 thread,
                 F,
                 p,
                 nt,
                 dm,
                 f,
                 dc = 0,
                 fc = 0
                 ):
        """
        Parameters
        ----------
        thread : TYPE
            DESCRIPTION.
        F : int, float
            Carga externa.
        dm : int, float
            Diametro médio da rosca.
        p : int, float
            passo axial.
        nt: int.
            Número de filetes de rosca engajados.
        f : int, float
            coeficiente de atrito.
    
        """
        self.thread = thread
        
        thread_type = get_thread(self.thread)    
        self.alpha = np.deg2rad(alpha_table[thread_type])
        
        
        self.F = F
        
        self.p = p
        
        self.nt = nt
        
        self.l = self.p*self.nt
        
        self.dm = dm
        self.dr = dm-p/2
        self.dext = dm+p/2
        
        self.f = f
        
        self.dc = dc
        
        self.fc = fc
        
    

    def nofriction(self):
        """
        Calcula o torque sem fricção no parafuso.
    
        Returns
        -------
        float
            Torque sem atrito do parafuso.
    
        """
        return (self.F*self.l)/(2*np.pi)
    
    def raise_torque(self):
            """
        Calcula o torque para levantar uma carga externa
        
        Returns
        -------
        Tr : float
            Torque para levantar a carga.
    
        """

            Tr = ((self.F*self.dm)/2)*((self.l+self.dm*np.pi*self.f/np.cos(self.alpha))/(np.pi*self.dm-self.f*self.l/np.cos(self.alpha))) + (self.F*self.fc*self.dc)/2 #Nm, Torque para subir o parafuso (descer a prensa)
            
            return Tr
    
    def raise_efficiency(self):
    
        T = self.raise_torque()
        T0 = self.nofriction()
        return T0/T
    
    def lower_torque(self):
            """
        Calcula o torque para abaixar uma carga externa
        
        -------
        Tl : float
            Torque para descer a carga.
    
        """

            Tl = ((self.F*self.dm)/2)*((self.dm*np.pi*self.f/np.cos(self.alpha)-self.l)/(np.pi*self.dm+self.f*self.l/np.cos(self.alpha))) + (self.F*self.fc*self.dc)/2 #Nm, Torque para subir o parafuso (descer a prensa)
            
            return Tl
        
    def lower_efficiency(self):
    
        T = self.lower_torque()
        T0 = self.nofriction()
        return T0/T
    
    
    def stress_tensor(self,
                      raise_load = True,
                      wheight_threads = False):
        """
        

        Parameters
        ----------
        raise_load : bool, optional
            Quandro True, inidca que o torque analisado é de subir a carga. The default is True.
        wheight_threads : bol, optional
            Quando True, indica que o partilhamento de cargas não é igual, de modo que a força no pimeiro 
            filete é 0,38*F. The default is False.

        Returns
        -------
        S : numpy array
            Tensor de tensões na raiz do parafuso.

        """
        if wheight_threads is True:
            F = 0.38*self.F
            nt = 1
        
        else:
            F = self.F
            nt = self.nt
        
        alpha = self.alpha
        dr = self.dr
        p = self.p
        
        t = p/2 * (1+np.tan(alpha))
        
        if raise_load is True:
            T = self.raise_torque()
        else:
            T = self.lower_torque()
            
        sig_xx = (3*F)/(np.pi*dr*t*nt)
        sig_yy = (4*F)/(np.pi*dr**2)
        sig_zz = 0
        
        tau_xy = 0
        tau_xz = 0
        tau_yz = (16*T)/(np.pi*dr**3)
        
        S = np.array([[sig_xx,tau_xy,tau_xz],
                      [tau_xy,sig_yy,tau_yz],
                      [tau_xz,tau_yz,sig_zz]]
                     )

        return S

