# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 17:00:54 2022

@author: valte
"""

import numpy as np
import pandas as pd
import math as mt
from ..Units import SI,rad_velocity,m2ft,ft2m,lbf2N,N2lbf

def Tr_from_H(H,n,unit = "RPM"):
    """
    Calcula o torque em Nm baseado na rotação e potência.

    Parameters
    ----------
    H : float
        Potência de eixo, em W (Entrar em kW implica saída em Nmm).
    n : float
        Rotação do eixo. Unidades RPM ou rad/s. Deve-se indicar a unidade utili
        zando o parâmetro `unit`
     unit : str, optional
        Unidade da rotação. The default is "RPM".

    Returns
    -------
    Tr : float
        Torque no eixo.

    """
    omega = rad_velocity(n,input_unit=unit,output_unit='rad/s')
    Tr = H/omega
    return Tr

def H_from_Tr(Tr,n,unit = "RPM"):
    """
    Calcula a potência no eixo a partir do torque

    Parameters
    ----------
    Tr : Float
        Torque em Nm (Valores em Nmm resultam em H em kW).
    n : float
        Rotação do eixo. Unidades RPM ou rad/s. Deve-se indicar a unidade utili
        zando o parâmetro `unit`
     unit : str, optional
        Unidade da rotação. The default is "RPM".

    Returns
    -------
    H : FLoat
        Potência.

    """
    omega = rad_velocity(n,input_unit=unit,output_unit='rad/s')   
    H = Tr*omega
    return H

def Wt_from_Tr(Tr,d):
    """
    Calcula a carga transversal na engrenagem a partir do torque

    Parameters
    ----------
    Tr : Float
        Torque.
    d : Float
        Diametro nominal da engrenagem.

    Returns
    -------
    W_t : Float
        Carga transversal.

    """
    W_t = 2*Tr/d
    return W_t

def Wt_from_H(H,d,n,unit="RPM"):
    """
    Força tranversal na engrenagem a partir da potência de eixo

    Parameters
    ----------
    H : float
        Potência de eixo, em W (Entrar em kW implica saída em Nmm).
    d : Float
        Diametro nominal da engrenagem.
    n : Float
        Rotação do eixo. Unidades RPM ou rad/s. Deve-se indicar a unidade utili
        zando o parâmetro `unit`
     unit : str, optional
        Unidade da rotação. The default is "RPM".

    Returns
    -------
    W_t : Float
        Força transversal.

    """
    omega = rad_velocity(n,input_unit=unit,output_unit='rad/s')
    W_t = 2*H/(omega*d)
    return W_t


def adendum(P):
    return 1/P

def dedendum(P):
    return 1.25/P

def pitch_diam(m,N):
    return m*N

def base_thicknes(m):
    P = 1/m
    p = np.pi/P
    
    return p/2

def min_pinion_N(m,k=1,phi=20):
    """
    Calcula o mínimo número de dentes no pinhão para evitar interferência

    Parameters
    ----------
    m : float
        Módulo transversa.
    k : TYPE, optional
        DESCRIPTION. The default is 1.
    phi : float, optional
        Angulo de pressão. The default is 20.

    Returns
    -------
    N: int
        Mínimo número de dentes do pinhão.

    """
    phi = np.deg2rad(phi)
    return mt.ceil((2*k/((1+2*m)*np.sin(phi)**2))*(m+np.sqrt(m**2+(1+2*m)*np.sin(phi)**2)))



class SpurAGMA():
    """
    Modela engrenagens cilíndricas de dentes retos e helicoidais conforme 
    metodologia da AGMA.

    Parameters
    ----------
    F : float
        Largura de face em mm.
    P_t : float
        Passo transversal em mm.
    N_p : int
        N° de dentes do Pinhão.
    N_g : int
        N° de dentes da Coroa.
    poisson_p : float
        Coeficiente de Poisson do Material do Pinhão.
    poisson_g : float
        Coeficiente de Poisson do Material da Coroa.
    E_p : float
        Módulo de Young em GPa do Material do Pinhão.
    E_g : float
        Módulo de Young em GPa do Material da Coroa.
    Q_v : float
        Número de qualidade.
    n_p : TYPE
        DESCRIPTION.
    n_g : TYPE
        DESCRIPTION.
    N_cyc : int
        Número de ciclos para vida em fadiga.
    H_Bp : float
        Dureza Brinnel do pinhão.
    H_Bg : float
        Dureza Brinnel da Coroa.
    N_c : int, optional
        Número de dentes da cremalheira de fabricação. The default is 9999.
    r_c : float, optional
        Raio nominal da cremalheria de fabricação. The default is None.
    a_0 : float, optional
        Adendo admensional da cremalheira de fabricação. The default is 1.25.
    b_0 : float, optional
        Dedendo adimensional da cremalheira de fabricação. The default is 1.
    rho_a0 : TYPE, optional
        Raio de adoçamento da cremalheira de fabricação. The default is 0.25.
    x_0 :float, optional
        Folga adimensional de "backlash" da cremalheira de fabricação. The default is 0.
    phi_n : float, optional
        Angulo de pressão em deg. The default is 20.
    psi : float, optional
        Angulo de hélice em deg. The default is 0.
    delta_s_n : float, optional
        Parametro de modificação do dente. The default is 0.024.
    x : float, optional
        Modificador de adendo do dente. The default is 0.
    load : str, optional
        Ponto de carregamento do dente: 'tip' para ponta ou 'HPSTC' para 
        o ponto mais alto de contato em um único dente. The default is "tip".
    orientation : str, optional
        Tipo de engrenagem. 'internal' para internas e 'external' 
        para externas. The default is 'external'.
    coroamento : bol, optional
        Define se há uso de coroamento nos dentes. The default is False.
    detrimento_superficie : bol, optional
        Define se há presenseça de detrimento superfiical. The default is False.
    engrenamento : str, optional
        Define o tipo de engrenamento.
         -1 ou 'aberto': engranamento aberto;
         -2 ou 'fecahdo comercial': Unidades fechadas comerciais;
         -3 ou 'fechado preciso': Unidades fechadas de precisão;
         -4 ou 'fechado extra-preciso': Unidades fechadas extra-precisas.
        O padrão é 'aberto'.
    montagem_ajustada: bool, optional
        Define o tipo de mntagem. O padrão é Falso.
    tratamento : str, optional
        Define o tipo de tratamento superficial. The default is 'None'.
    fonte_potencia : str, optional
        Define o nível de choque esperado na fonte. The default is 'uniforme'.
    choque_maquina : str, optional
        Define o nível de choque esperado na máquina. The default is 'uniforme'.
    trust : float, optional
        Define o nível de confialbilidade adotado. The default is 99.99.
    temperatura : float, optional
        Temperatura de operação em °C. The default is 20.



        """
    def __init__(self,F,P_t,N_p,N_g,
                 poisson_p,poisson_g,E_p,E_g,
                 Q_v, n_p, n_g, N_cyc,
                 H_Bp,H_Bg,
                 N_c=9999,r_c=None,a_0 = 1.25,b_0=1,rho_a0=0.25,x_0 = 0,
                 phi_n=20,psi=0,
                 delta_s_n=0.024,x=0,
                 load = "tip", orientation = 'external',
                 coroamento = False, detrimento_superficie = False, 
                 engrenamento = 'aberto',montagem_ajustada=False,tratamento = 'None',
                 fonte_potencia = 'uniforme', choque_maquina = 'uniforme',
                 trust = 99.99, temperatura = 20):


            
        self.P = P_t
        self.F = F
        self.N_g = N_g
        self.N_p = N_p
        self.poisson_p = poisson_p
        self.poisson_g = poisson_g
        self.E_p = E_p*1e3
        self.E_g = E_g*1e3
        self.Q_v = Q_v
        self.n_p = n_p
        self.N_cyc = N_cyc
        self.H_Bp = H_Bp
        self.H_Bg = H_Bg
        self.phi_n = np.deg2rad(phi_n)
        self.psi = psi
        self.x = x
        self.delta_s_n = delta_s_n
        self.N_c = N_c
        self.rho_a0 = rho_a0
        self.loaf = load
        self.coroamento = coroamento 
        self.engrenamento = engrenamento
        self.montagem_ajustada = montagem_ajustada
        
        phi_n = np.deg2rad(phi_n)
        psi = np.deg2rad(psi)

        phi = np.arctan2(np.tan(phi_n),np.cos(psi))
        P = P_t/np.cos(psi)
        
        F = F*P #adimensional
        delta_s_n = F*delta_s_n
        
        #if psi == 0:
        m_g = N_g/N_p
        
        n_g = n_p/m_g
        
        self.n_g = n_g
        
        if r_c == None:
            r_c = N_c/(2*np.cos(psi))
            r_bc = r_c-b_0

        
        "Dimensões Pinhão"
        r_n = N_p/(2*np.cos(psi))
        a = P*adendum(P)
        b = P*dedendum(P)

        x_g = x - (delta_s_n)/(2*np.tan(phi_n))
        s_n = np.pi/2 + 2*x_g*np.tan(phi_n)
        
        r_a = (N_p/np.cos(psi) + 2*(1+x))/2
        r_b = r_n*np.cos(phi)

        self.r_n = r_n
        self.a = a
        self.b = b
        self.r_a = r_a

        self.r_b = r_b

        self.s_n = s_n
        
        
        
        "Dimensões Coroa"
        r_n2 = r_n*m_g
        a_2 = 1
        b_2 = 1.25
        c_2 = b_2-a_2
        r_r2 = r_n2-b_2
        r_b2 = r_b*m_g
        r_f2 = r_r2+c_2
        
        r_a2 = (N_g/np.cos(psi) + 2*(1-x))/2
        
        self.r_n2 = r_n2
        self.a_2 = a_2
        self.b_2 = b_2
        self.r_a2 = r_a2
        self.r_r2 = r_r2
        self.r_b2 = r_b2
        self.r_f2 = r_f2

        

        
        "Perfil de Involuta"
        
        C_r = (N_p+N_g)/(2*np.cos(psi))
        phi_r = (np.arccos((r_b2+r_b)/C_r) if orientation == 'external' 
                 else np.arccos((r_b2-r_b)/C_r))
        p_n = np.pi*np.cos(phi_n)
        p_b = 2*np.pi*r_b/N_p
        psi_b = np.arccos(p_n/p_b)
        
        
        C_6 = C_r*np.sin(phi_r)
        C_1 = (C_6 - (r_a2**2-r_b2**2)**0.5 if orientation == 'external' 
               else -(C_6 - (r_a2**2-r_b2**2)**0.5))
        C_3 = (C_6/(m_g+1) if orientation == 'external' 
               else C_6/(m_g-1))
        C_4 = C_1+p_b
        C_5 = (r_a**2-r_b**2)**0.5
        C_2 = C_5-p_b
        
        Z = C_5-C_1
        m_p = Z/p_b
        
        p_x = np.pi/np.sin(psi) if psi!=0 else 0
        
        m_F = F/p_x if psi!=0 else 0
        
        self.phi_r = phi_r
        self.psi_b = psi_b
        self.C_1 = C_1
        self.C_2 = C_2
        self.C_3 = C_3
        self.C_4 = C_4
        self.C_5 = C_5
        self.C_6 = C_6
        self.C_r = C_r
        self.Z= Z
        self.m_F = m_F
        self.m_p = m_p
        self.p_x = p_x
        
        "Comprimento Mínimo das Linhas de Contato"
        
        n_r = m_p%1
        n_a = m_F%1
        
        if psi == 0: #para engrenagens de dentes retor
            L_min = F
        
        elif n_a<=1-n_r: #engrenagens helicoidais
            L_min = (m_p*F-n_a*n_r*p_x)/np.cos(psi_b)
        else:
            L_min = (m_p*F-(1-n_a)*(1-n_r)*p_x)/np.cos(psi_b)
          
        "Partilhamento de Carga" 
        m_N = F/L_min
        
        "Angulo de Hélice de Operação"
        psi_r = np.arctan(np.tan(psi_b)/np.cos(phi_r))
        phi_nr = np.arcsin(np.cos(psi_b)*np.sin(phi_r))
        
        
        "Fator geométrico de resistência superficial I"
        d = ((2*C_r)/(m_g+1) if orientation == 'external'
             else (2*C_r)/(m_g-1))
        
        r_m1 = ((1/2)*(r_a+(C_r-r_a2)) if orientation== 'external'
                else (1/2)*(r_a-(C_r-r_a2)))
        
        "Raios de Giração"
        if psi!=0 and m_F > 1:
            rho_1 = (r_m1**2 - r_b**2)**0.5
            rho_2 = (C_6 - rho_1 if orientation=='external'
                     else C_6 + rho_1 )
        
        else:
            rho_1 = C_2
            rho_2 = (C_6 - rho_1 if orientation == 'external'
                     else C_6 + rho_1)
        self.rho_1 = rho_1
        self.rho_2 = rho_2
        
        "Fator de Sobreposição Helicoidal"
        if psi!=0 and m_F < 1:
            rho_m1 = (r_m1**2 - r_b**2)**0.5
            rho_m2 = (C_6 - rho_1 if orientation=='external'
                     else C_6 + rho_1)
            
            C_psi = (1-m_F*(1-(rho_m1*rho_m2*Z)/(rho_1*rho_2*p_n)))**0.5
            self.rho_m1 = rho_m1
            self.rho_m2 = rho_m2
            
        else:
            C_psi = 1
        
        
        self.C_psi = C_psi
        
        
        I = ((np.cos(phi_r)*C_psi**2)/
             (((1/rho_1)+(1/rho_2))*d*m_N) if orientation == 'external'
             else (np.cos(phi_r)*C_psi**2)/
             (((1/rho_1)-(1/rho_2))*d*m_N))
          
        self.I = I
        

        
        "Fator de geométrico de flexão J"
        if psi!=0:
            N_1v = N_p/(np.cos(psi)**3)
            r_1v = N_1v/2
            r_b1v = r_1v*np.cos(phi_n)
            r_a1v = r_1v+r_a-r_n
        else:
            N_1v = N_p
            r_1v = r_n
            r_b1v = r_b
   
            
        r_2v = r_1v*m_g
        r_b2v = r_b1v*m_g
        r_a2v = r_2v+r_a2-r_n2
            
            
        C_6v = (r_b2v+r_b1v)*np.tan(phi_nr)
        C_1v = (C_6v-(r_a2v**2-r_b2v**2)**0.5)
        C_4v = C_1v+p_n
        

        
        self.phi_nr =phi_nr
        self.r_1v =r_1v
        self.r_b1v =r_b1v
        self.r_2v =r_2v
        self.r_a2v =r_a2v
        self.r_b2v =r_b2v
        self.C_1v = {}
        self.C_4v = {}
        self.C_6v = {} 
        
        pinion_geometry = [N_p,r_n,r_a,r_b,n_p,'pinnion']
        gear_geometry = [N_g,r_n2,r_a2,r_b2,n_g,'gear']
        geometry_list = [pinion_geometry,
                         gear_geometry]

        self.m_p = {}
        self.tan_phi_nW = {}
        self.r_bv = {}
        self.phi_nL = {}
        self.phi_line = {}
        self.phi_n = {}
        self.phi_ni = {}
        self.phi_ns = {}
        self.lambda_nF = {}
        self.alpha = {}
        self.mu_n0 = {}
        self.r_nL = {}
        self.r_n_line = {}
        self.r_n0_line = {}
        self.lambd = {}
        self.theta_n0 = {} 
        self.theta_n = {} 
        self.beta_n = {} 
        self.eta_nF = {} 
        self.xi_nF = {} 
        self.h_f = {} 
        self.K_F = {} 
        self.K_S = {} 
        self.convergence_flag = {}       
        self.rho_min = {} 
        self.s_f = {} 
        self.H = {} 
        self.L = {} 
        self.M = {} 
        self.K_f = {} 
        self.Y = {} 
        self.J = {}
        self.K_v = {}
        self.V = {}
        self.V_max = {}
        self.Y_lewis = {}
        self.K_s = {}
        
        for geometry in geometry_list:
            N_i = geometry[0]
            r_i = geometry[1]
            r_ai = geometry[2]
            r_bi = geometry[3]
            n_i = geometry[4]
            tag = geometry[5]
            
            #Angulo de Pressão no Ponto de Aplicação
            if (psi !=0 and m_F>1) or (psi == 0 and load == "tip"):
                tan_phi_nW = ((r_a1v/r_b1v)**2-1)**0.5
            else:
                tan_phi_nW = C_4v/r_b1v  
                
            #Geometria virtual
            if psi!=0:
                N_v = N_i/(np.cos(psi)**3)
                r_v = N_v/2
                r_bv = r_v*np.cos(phi_n)

            else:
                N_v = N_i
                r_v = r_i
                r_bv = r_bi

            "Virtual cutter tool, cremalheira - appendix D"
            
            if psi != 0:
                N_0 = N_c/(np.cos(psi)**3)
                r_n0 = N_0/2
                r_nb0 = r_n0*np.cos(phi_n)

            else:
                N_0 = N_c
                r_n0 = r_c
                r_nb0 = r_bc
            

            self.x_0 = x_0
            s_n0= np.pi/2+2*x_0*np.tan(phi_n)
            
            r_s0 = r_n0 + a_0 + x_0 - rho_a0
            
            
            self.r_n0 = r_n0
            self.a_0 = a_0
            self.b_0 = b_0
            self.r_s0 = r_s0
            self.r_nb0 = r_nb0
            self.s_n0 = s_n0        
                
            
            "Angulos de Geração"
            inv_phi_n = np.tan(phi_n) - phi_n
            
            inv_phi_np = inv_phi_n + s_n/(2*r_v)
            
            phi_nL = tan_phi_nW - inv_phi_np
            
            r_nL = r_bv/np.cos(phi_nL)

            phi_ns = np.arccos(r_nb0/r_s0)
            
            inv_phi_np0 = inv_phi_n + s_n0/(2*r_n0)
            
            inv_phi_ns = np.tan(phi_ns)-phi_ns
            
            lambd = 2*(inv_phi_np0-inv_phi_ns-rho_a0/r_nb0)
            
            inv_phi_line = inv_phi_n + (2*(x_g+x_0)*np.tan(phi_n))/(N_v+N_0)
            phi_ni = (3*inv_phi_line)**0.33

            for i in range(10):
                phi_line = phi_ni + (inv_phi_line+phi_ni-np.tan(phi_ni))/(np.tan(phi_ni)**2)
                phi_ni = phi_line
            
            r_n_line = r_v*np.cos(phi_n)/np.cos(phi_line)
            r_n0_line = r_n0*np.cos(phi_n)/np.cos(phi_line)
                     
   
            "Curva de lewis"
            
            alpha_i = np.pi/2
            y = 1
            
            alpha_counter = 0
            max_it = 15
            convergence_flag = 0
            while abs(y) > 1e-3 and alpha_counter<=max_it:
                
                mu_n0 = np.arccos(r_n0_line*np.cos(alpha_i)/r_s0)-alpha_i
                
                K_S = r_n0_line*np.sin(alpha_i) - r_s0*np.sin(alpha_i+mu_n0)
                
                K_F = K_S - rho_a0
                
                theta_n0 = mu_n0 - lambd/2 + np.pi/N_0
                
                theta_n = N_0/N_v*theta_n0
                
                beta_n = alpha_i - theta_n
                
                xi_nF = r_n_line*np.sin(theta_n)+K_F*np.cos(beta_n)
                
                eta_nF = r_n_line*np.cos(theta_n)+K_F*np.sin(beta_n)
                
                #Altura da curva de Lewis
                h_f = r_nL - eta_nF
                
                y = 2*h_f*np.tan(beta_n) - xi_nF
                y_line = (2*h_f/(np.cos(beta_n)**2) - K_F*np.sin(beta_n)
                        +N_0/N_v*((r_n0_line*np.sin(alpha_i)/(r_s0*np.sin(alpha_i+mu_n0)))-1)
                        *(2*xi_nF*np.tan(beta_n)-eta_nF-(2*h_f)/np.cos(beta_n)**2)
                        -r_n0_line*(np.cos(alpha_i)-np.sin(alpha_i)/np.tan(alpha_i+mu_n0))
                        *((1+np.sin(beta_n)**2)/np.cos(beta_n)))                    
                
                alpha = alpha_i - y/y_line
                alpha_i = alpha
                alpha_counter +=1
                
                if alpha_counter >= max_it-1:
                    #print('Erro de Convergência para Alpha!!')
                    convergence_flag = 1

            
            #raio de adoçamento mínimo
            rho_min = (rho_a0 + 
                       ((r_n0_line-r_s0)**2)/
                       (r_n_line*r_n0_line/
                        (r_n_line+r_n0_line)-(r_n0_line-r_s0)))

            
            #Fator helicoidal
            if psi!=0 and m_F>1:
                omega = np.rad2deg(np.arctan(np.tan(psi)*np.sin(phi_n)))
                C_h = 1/(1-((omega/100)*(1-omega/100)))**0.5
            else:
                C_h = 1
                
            #Espessura crítica da curva de Lewis
            s_f = abs(2*xi_nF)

            
            #Fator de Forma:
            H = 0.331-0.436*phi_n
            L = 0.324-0.492*phi_n
            M = 0.261+0.545*phi_n

            
            K_f = H +((s_f/rho_min)**L)*((s_f/h_f)**M)

            #Fator de ANgulo de Hélice
            if psi!=0 and m_F>1:
                K_psi = np.cos(psi_r)*np.cos(psi) 
            else:
                K_psi = 1
                     
            #Fator Y:
            Y = K_psi/((np.cos(phi_nL)/np.cos(phi_nr))
                   *((6*h_f/(s_f**2*C_h))-np.tan(phi_nL)/s_f))
            

            #Fator de Forma:
            J = (Y*C_psi/(K_f*m_N))
      
            
            "Fator Dinâmico K_v"
            
            V = abs(np.pi*n_i*(r_n/P)/6)
            B = 0.25*(12-Q_v)**(2/3)
            A = 50+56*(1-B)
            
            K_v = ((A+np.sqrt(V))/A)**B
            
            V_max = (A + (Q_v - 3))**2
            
            
            "Fator de Tamanho"

            F = F/P
            
            u = ((s_f/P)**2)/(4*h_f/P)
            Y_lewis  = 2*u*P/3
            
            #Shigley:
            K_s = 1.192*((F*Y_lewis**0.5)/P)**0.0535
            
            #AGMA
            if K_s < 1:
                K_s = 1

            self.C_1v[tag] = C_1v
            self.C_4v[tag] = C_4v
            self.C_6v[tag] = C_6v
            self.m_p[tag] = m_p
            self.tan_phi_nW[tag] = tan_phi_nW
            self.r_bv[tag] = r_bv
            self.phi_nL[tag] = phi_nL
            self.phi_line[tag] = phi_line
            self.phi_n[tag] = phi_n
            self.phi_ni[tag] = phi_ni
            self.phi_ns[tag] = phi_ns
            self.lambda_nF[tag] = lambd
            self.alpha[tag] = alpha
            self.mu_n0[tag] = mu_n0
            self.r_nL[tag] = r_nL
            self.r_n_line[tag] = r_n_line
            self.r_n0_line[tag] = r_n0_line
            self.lambd[tag] = lambd
            self.theta_n0[tag] = theta_n0
            self.theta_n[tag] = theta_n
            self.beta_n[tag] = beta_n
            self.eta_nF[tag] = eta_nF
            self.xi_nF[tag] = xi_nF
            self.h_f[tag] = h_f
            self.K_F[tag] = K_F
            self.K_S[tag] = K_S
            self.convergence_flag[tag] = convergence_flag         
            self.rho_min[tag] = rho_min
            self.s_f[tag] = s_f
            self.H[tag] = H
            self.L[tag] = L
            self.M[tag] = M
            self.K_f[tag] = K_f
            self.Y[tag] = Y
            self.J[tag] = J
            self.K_v[tag] = K_v
            self.V[tag] = V
            self.V_max[tag] = V_max
            self.Y_lewis[tag] = Y_lewis
            self.K_s[tag] = K_s
        
        
        "Fator de Condição de Superfície"
        if detrimento_superficie == True:
            C_f = 1.5
        else:
            C_f = 1
        
        self.C_f = C_f
        
        "Fator de Distribuição de Carga K_m"
        lookup = {1:'aberto',
                  2:'fechado comercial',
                  3:'fechado preciso',
                  4:'fechado extra-preciso'
                  }
        
        if type(engrenamento) is int:
            engrenamento = lookup[engrenamento]
        
        condition = ((F/(2*r_n/P)) <= 2) and (F <= 40)

        if condition == True:
            C_mc = 1 if coroamento == False else 0.8
            C_pf = (F/(10*2*r_n/P) - 0.025 if F <= 1
                    else F/(10*2*r_n/P) - 0.0375 + 0.0125*F if F<= 17                           
                    else F/(10*2*r_n/P) - 0.1109 + 0.0207*F - 0.00338* F**2)
            C_pm = 1 #ALTERAR PARA A PARTE 2
            C_ma_dic = {'aberto' : [0.247,0.0167,-0.765e-4],
                        'fechado comercial' : [0.127,0.0158,-0.930e-4],
                        'fechado preciso' : [0.0675,0.0128,-0.926e-4],
                        'fechado extra-preciso' : [0.00360,0.0102,-0.822e-4]}
            A = C_ma_dic[engrenamento.lower()][0]
            B = C_ma_dic[engrenamento.lower()][1]
            C = C_ma_dic[engrenamento.lower()][2]              
            C_ma = A+B*F+C*F**2
            C_e = 0.8 if montagem_ajustada else 1

            K_m = 1+C_mc*(C_pf*C_pm+C_ma*C_e)
            
        else:
            K_m = 1
            
        
        self.K_m = K_m
            
        "Fator de Razaõ de Dureza"
        if m_g >= 1:
            C_H = 1
            
                             
        elif 1/(H_Bp/H_Bg) >= 1.2: 
            if 1/(H_Bp/H_Bg) > 1.7:
                A_line = 0.00698
            else:
                A_line = (9.98e-3*(H_Bp/H_Bg)-8.29e-3)
            C_H = 1 + A_line*(1/m_g - 1)
            
        else:
            C_H = 1
        
        self.C_H = C_H
        
        "Fator de Sobrecarga"
        over_dic = {'uniforme': {'uniforme': 1,
                                 'leve': 1.25,
                                 'intenso': 1.75},
                    'leve' : {'uniforme': 1.25,
                                 'leve': 1.5,
                                 'intenso': 2},
                    'medio': {'uniforme': 1.5,
                                 'leve': 1.5,
                                 'intenso': 2.25}}
        
        pot_key = fonte_potencia.lower().replace('é','e')
        machine_key = choque_maquina.lower().replace('é','e')
        
        if len(pot_key.split())>1:
            pot_key = pot_key.split()[1]
            
        if len(machine_key.split())>1:
            machine_key = machine_key.split()[1]
            
        K_0 = over_dic[pot_key][machine_key]
        
        self.K_0 = K_0
        
        "Fatores de Ciclagem de Tensão"
        
        #Y_N
        Y_N = {}
        
        if N_cyc < 3e6:
            if tratamento.lower() == 'nitretado':
                Y_N['pinnion'] = 4.9404*N_cyc**-0.1045
            elif tratamento.lower()[0:11] == 'carbonetado':
                Y_N['pinnion'] = 6.1514*N_cyc**-0.1192
            elif abs(H_Bp-160) < abs(H_Bp-250) and abs(H_Bp-160) < abs(H_Bp-400):
                Y_N['pinnion'] = 2.3194*N_cyc**-0.0538
            elif abs(H_Bp-250) < abs(H_Bp-160) and abs(H_Bp-250) < abs(H_Bp-400):
                Y_N['pinnion'] = 4.9404*N_cyc**-0.1045
            else:
                Y_N['pinnion'] = 9.4518*N_cyc**-0.148
        else:
            Y_N['pinnion'] = 1.3558*N_cyc**-0.0178 #Mais conservadora
        
        if N_cyc < 3e6:
            if tratamento.lower() == 'nitretado':
                Y_N['gear'] = 4.9404*N_cyc**-0.1045
            elif tratamento.lower()[0:11] == 'carbonetado':
                Y_N['gear'] = 6.1514*N_cyc**-0.1192
            elif abs(H_Bg-160) < abs(H_Bg-250) and abs(H_Bg-160) < abs(H_Bg-400):
                Y_N['gear'] = 2.3194*N_cyc**-0.0538
            elif abs(H_Bg-250) < abs(H_Bg-160) and abs(H_Bg-250) < abs(H_Bg-400):
                Y_N['gear'] = 4.9404*N_cyc**-0.1045
            else:
                Y_N['gear'] = 9.4518*N_cyc**-0.148
        else:
            Y_N['gear'] = 1.3558*N_cyc**-0.0178 #Mais conservadora
            
            
        self.Y_N = Y_N
            
        #Z_N
        if N_cyc < 1e7:
            if tratamento.lower() == 'nitretado':
                Z_N = 1.249*N_cyc**-0.0138
            else:
                Z_N = 2.466*N_cyc**-0.056
        else:
            Z_N = 1.4488*N_cyc**-0.023 #Mais conservadora
        
        self.Z_N = Z_N
            
        "Coeficiente Elástico C_p"
        
        C_p = np.sqrt(1/(np.pi*((1-poisson_p**2)/(self.E_p)+(1-poisson_g**2)/(self.E_g))))
        
        self.C_p = C_p 
        
        "Fator de Confiabilidade"
        
        if trust == 99.99:
            K_r = 1.5
        elif trust == 99.9:
            K_r = 1.25
        elif trust == 99:
            K_r = 1
        elif trust == 90:
            K_r = 0.85
        elif trust > 50 and trust <99:
            K_r = 0.658-0.0759*np.log(1-trust/100)
        elif trust > 99:
            K_r = 0.5-0.109*np.log(1-trust/100)
        else:
            print('Confiabilidade Muito Baixa!!!')
        
        self.K_r = K_r
        
        "Fator de Temperatura"
        if temperatura <= 120:
            K_T = 1
        else:
            K_T = 1.5
        
        self.K_T = K_T
        
        "Espessura de Aro"
        self.K_B = 1
            
    def bending_stress_g(self,W_t):
        """
        Tensão de flexão na coroa

        Parameters
        ----------
        W_t : float
            Carga transversal em N.

        Returns
        -------
        sigma : float
            Tensão em Pa.

        """
        K_0 = self.K_0
        K_v = self.K_v['gear']
        K_s =  self.K_s['gear']
        P = self.P
        F = self.F
        K_m =self.K_m
        K_B = self.K_B
        J =self.J['gear']
        sigma = W_t*K_0*K_v*K_s*(P/F)*(K_m*K_B/J)
        return sigma
    
    def bending_stress_p(self,W_t):
        """
        Tensão de flexão no pinhão

        Parameters
        ----------
        W_t : float
            Carga transversal em N.

        Returns
        -------
        sigma : float
            Tensão em Pa.

        """
        K_0 = self.K_0
        K_v = self.K_v['pinnion']
        K_s =  self.K_s['pinnion']
        P = self.P
        F = self.F
        K_m =self.K_m
        K_B = self.K_B
        J =self.J['pinnion']
        sigma = W_t*K_0*K_v*K_s*(P/F)*(K_m*K_B/J)
        return sigma
    
    
    def contact_stress_g(self,W_t):
        """
        Tensão de contato na coroa

        Parameters
        ----------
        W_t : float
            Carga transversal em N.

        Returns
        -------
        sigma : float
            Tensão em Pa.

        """
        K_0 = self.K_0
        K_v = self.K_v['gear']
        K_s =  self.K_s['gear']
        P = self.P
        F = self.F
        K_m =self.K_m
        C_p = self.C_p
        C_f = self.C_f
        r_p = self.r_n/P
        I =self.I
        sigma = C_p*(W_t*K_0*K_v*K_s*(K_m/(2*r_p*F))*(C_f/I))**0.5
        return sigma
    
    def contact_stress_p(self,W_t):
        """
        Tensão de contato no pinhão

        Parameters
        ----------
        W_t : float
            Carga transversal em N.

        Returns
        -------
        sigma : float
            Tensão em Pa.

        """
        K_0 = self.K_0
        K_v = self.K_v['pinnion']
        K_s =  self.K_s['pinnion']
        P = self.P
        F = self.F
        K_m =self.K_m
        C_p = self.C_p
        C_f = self.C_f
        r_p = self.r_n/P
        I =self.I
        sigma = C_p*(W_t*K_0*K_v*K_s*(K_m/(2*r_p*F))*(C_f/I))**0.5
        return sigma
    
    def max_bending_g(self,S_t,S_F):
        """
        Máxima tensão de flexão na coroa

        Parameters
        ----------
        S_t : float
            Limite de resistência à flexão.
        S_F :float
            Fator de segurança.

        Returns
        -------
        max_sigma : float
                Tensão em Pa.

        """
        Y_N = self.Y_N['gear']
        K_T = self.K_T
        K_r = self.K_r
        max_sigma = (S_t/S_F)*Y_N/(K_T*K_r)
        return max_sigma
    
    def max_bending_p(self,S_t,S_F):
        """
        Máxima tensão de flexão admissível no pinhão

        Parameters
        ----------
        S_t : float
            Limite de resistência à flexão.
        S_F :float
            Fator de segurança.

        Returns
        -------
        max_sigma : float
                Tensão em Pa.

        """
        Y_N = self.Y_N['pinnion']
        K_T = self.K_T
        K_r = self.K_r
        max_sigma = (S_t/S_F)*Y_N/(K_T*K_r)
        return max_sigma
    
    def max_contact_g(self,S_c,S_F):
        """
        Máxima tensão de contato admissível na coroa

        Parameters
        ----------
        S_t : float
            Limite de resistência ao contato.
        S_F :float
            Fator de segurança.

        Returns
        -------
        max_sigma : float
                Tensão em Pa.

        """
        Z_N = self.Z_N
        C_H = self.C_H
        K_T = self.K_T
        K_r = self.K_r
        max_sigma = (S_c/S_F)*(Z_N*C_H)/(K_T*K_r)
        return max_sigma
    
    def max_contact_p(self,S_c,S_F):
        """
        Máxima tensão de contato admissível no pinhão

        Parameters
        ----------
        S_t : float
            Limite de resistência ao contato.
        S_F :float
            Fator de segurança.

        Returns
        -------
        max_sigma : float
                Tensão em Pa.

        """
        Z_N = self.Z_N
        C_H = 1
        K_T = self.K_T
        K_r = self.K_r
        max_sigma = (S_c/S_F)*(Z_N*C_H)/(K_T*K_r)
        return max_sigma
    
    def safety_bending_g(self,S_t,W_t):
        """
        Fator de segurança

        Parameters
        ----------
        S_t : Float
            Tensão máxima de flexão admissível em MPa.
        W_t : float
            Carga transversal no dente em N.

        Returns
        -------
        S_F : float
            Fator de Segurança a Flexão da Coroa.

        """
        Y_N = self.Y_N['gear']
        K_T = self.K_T
        K_r = self.K_r
        sigma = self.bending_stress_g(W_t)
        S_F = ((S_t*Y_N)/(K_T*K_r))/sigma
        return S_F

    def safety_bending_p(self,S_t,W_t):
        """

        Parameters
        ----------
        S_t : Float
            Tensão máxima de flexão admissível em MPa.
        W_t : float
            Carga transversal no dente em N.

        Returns
        -------
        S_F : float
            Fator de Segurança a Flexão do Pinhão.

        """
        Y_N = self.Y_N['pinnion']
        K_T = self.K_T
        K_r = self.K_r
        sigma = self.bending_stress_p(W_t)
        S_F = ((S_t*Y_N)/(K_T*K_r))/sigma
        return S_F
    
    def safety_contact_g(self,S_c,W_t):
        """

        Parameters
        ----------
        S_c : Float
            Tensão máxima de contato admissível em MPa.
        W_t : float
            Carga transversal no dente em N.

        Returns
        -------
        S_F : float
            Fator de Segurança ao Contato da Coroa.

        """
        Z_N = self.Z_N
        C_H = self.C_H
        K_T = self.K_T
        K_r = self.K_r
        sigma_c = self.contact_stress_g(W_t)
        S_F = ((S_c*Z_N*C_H)/(K_T*K_r))/sigma_c
        return S_F
    
    def safety_contact_p(self,S_c,W_t):
        """
        

        Parameters
        ----------
        S_c : Float
            Tensão máxima de contato admissível em MPa.
        W_t : float
            Carga transversal no dente em N.

        Returns
        -------
        S_F : float
            Fator de Segurança ao Contato do Pinhão.

        """
        Z_N = self.Z_N
        C_H = 1
        K_T = self.K_T
        K_r = self.K_r
        sigma_c = self.contact_stress_p(W_t)
        S_F = ((S_c*Z_N*C_H)/(K_T*K_r))/sigma_c
        return S_F
    
    def report_agma(self,precision=2,latex=False):
        params_dict = {'Parâmetros':['J',
                                     'I',
                                     's_f',
                                     'h_f',
                                     'Y',
                                     'K_s',
                                     'K_v',
                                     'K_0',
                                     'K_m',
                                     'C_p',
                                     'C_f',
                                     'K_T',
                                     'K_r',
                                     'C_H',
                                     'Z_N',
                                     'Y_l',
                                     'K_B'],
                       'Pinhão':[f'{self.J["pinnion"]:.{precision}f}',
                                f'{self.I:.{precision}f}',
                                f'{self.s_f["pinnion"]:.{precision}f}',
                                f'{self.h_f["pinnion"]:.{precision}f}',
                                f'{self.Y["pinnion"]:.{precision}f}',
                                f'{self.K_s["pinnion"]:.{precision}f}',
                                f'{self.K_v["pinnion"]:.{precision}f}',
                                f'{self.K_0:.{precision}f}',
                                f'{self.K_m:.{precision}f}',
                                f'{self.C_p:.{precision}f}',
                                f'{self.C_f:.{precision}f}',
                                f'{self.K_T:.{precision}f}',
                                f'{self.K_r:.{precision}f}',
                                f'{self.C_H:.{precision}f}',
                                f'{self.Z_N:.{precision}f}',
                                f'{self.Y_lewis["pinnion"]:.{precision}f}',
                                f'{self.K_B:.{precision}f}'],
                       'Coroa':[f'{self.J["gear"]:.{precision}f}',
                                f'{self.I:.{precision}f}',
                                f'{self.s_f["gear"]:.{precision}f}',
                                f'{self.h_f["gear"]:.{precision}f}',
                                f'{self.Y["gear"]:.{precision}f}',
                                f'{self.K_s["gear"]:.{precision}f}',
                                f'{self.K_v["gear"]:.{precision}f}',
                                f'{self.K_0:.{precision}f}',
                                f'{self.K_m:.{precision}f}',
                                f'{self.C_p:.{precision}f}',
                                f'{self.C_f:.{precision}f}',
                                f'{self.K_T:.{precision}f}',
                                f'{self.K_r:.{precision}f}',
                                f'{self.C_H:.{precision}f}',
                                f'{self.Z_N:.{precision}f}',
                                f'{self.Y_lewis["gear"]:.{precision}f}',
                                f'{self.K_B:.{precision}f}'],
                       }
        
        params_dt = pd.DataFrame.from_dict(params_dict)
        
        if latex:
            return params_dt.to_latex(index=False,
                                      formatters={'Parâmetros':lambda x: f'${x}$',
                                                  'Pinhão': lambda x: f'{x.replace(".",",")}',
                                                  'Coroa': lambda x: f'{x.replace(".",",")}'
                                                  },
                                      escape=False)
        
        else:
            return params_dt
    
    def report_stress(self,S_t,S_c,W_t,precision=2,latex=False):
        """
        

        Parameters
        ----------
        S_t : TYPE
            DESCRIPTION.
        S_c : TYPE
            DESCRIPTION.
        W_t : TYPE
            DESCRIPTION.
        precision : TYPE, optional
            DESCRIPTION. The default is 2.
        latex : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        params_dt: DataFrame or String
            Dataframe com os dados das tensões.

        """
        
        sig_Bp = f'{self.bending_stress_p(W_t):.{precision}f}'
        sig_Bg = f'{self.bending_stress_g(W_t):.{precision}f}'
        FS_Bp = f'{self.safety_bending_p(S_t, W_t):.{precision}f}'
        FS_Bg = f'{self.safety_bending_g(S_t, W_t):.{precision}f}'
        
        sig_Cp = f'{self.contact_stress_p(W_t):.{precision}f}'
        sig_Cg = f'{self.contact_stress_g(W_t):.{precision}f}'
        FS_Cp = f'{self.safety_contact_p(S_c, W_t):.{precision}f}'
        FS_Cg = f'{self.safety_contact_g(S_c, W_t):.{precision}f}'

        params_dict = {'Parâmetros': ['sigma_B [MPa]',
                                      'sigma_C [MPa]',
                                      'FS_B',
                                      'FS_C'],
                       'Pinhão': [sig_Bp,
                                  sig_Cp,
                                  FS_Bp,
                                  FS_Cp],
                       'Coroa': [sig_Bg,
                                 sig_Cg,
                                 FS_Bg,
                                 FS_Cg]
                       }
        
        params_dt = pd.DataFrame.from_dict(params_dict)
        
        if latex:
            return params_dt.to_latex(index=False,
                                      formatters={'Parâmetros':lambda x: f'${x}$'.replace('sigma','\sigma'),
                                                  'Pinhão': lambda x: f'{x.replace(".",",")}',
                                                  'Coroa': lambda x: f'{x.replace(".",",")}'
                                                  },
                                      escape=False)
        
        else:
            return params_dt
        
    def report_dims(self,precision=2,latex=False):
        m = 1/self.P
        N_p = self.N_p
        N_g = self.N_g
        d_p = N_p*m
        d_g = N_g*m
        F = self.F
        
        dims_dict = {'Dimensões':['N',
                                  'd [mm]',
                                  'm [mm]',
                                  'F [mm]'],
                     'Pinhão': [f'{N_p:.0f}',
                                f'{d_p:.{precision}f}',
                                f'{m:.{precision}f}',
                                f'{F:.{precision}f}'],
                     'Coroa': [f'{N_g:.0f}',
                               f'{d_g:.{precision}f}',
                               f'{m:.{precision}f}',
                               f'{F:.{precision}f}'],
                     }
        
        
        dims_dt = pd.DataFrame.from_dict(dims_dict)
        
        if latex:
            return dims_dt.to_latex(index=False,
                                    formatters={'Dimensões':lambda x: f'${x}$',
                                                 'Pinhão': lambda x: f'{x.replace(".",",")}',
                                                 'Coroa': lambda x: f'{x.replace(".",",")}'
                                                  },
                                    escape=False)
        
        else:
            return dims_dt
        
        
class WormAGMA():
    def __init__(self,
                 nw:float,
                 dw:float,
                 dc:float,
                 Nw:int,
                 Nc:int,
                 px:float,
                 lambd:float,
                 casting_method:str = 'molde',
                 **threadkwargs):
        
        """
            Define os parâmetros de projeto para dimensionamento de pinhões 
            sem-fim de acorod com a metodologia da AGMA.
    
            Parameters
            ----------
            nw : float
                Rotação do pinhão sem-fim em RPM.
            dw : float
                Diametro nominal do pinhão sem-fim em mm.
            dc : float
                Diametro nominal da coroa sem-fim em mm.
            Nw : int
                N° de dentes do pinhão sem-fim em mm.
            Nc : int
                N° de dentes da coroa sem-fim em mm.
            px: float
                passo axial do pinhão em mm
            lambd : float
                Ângulo de avanço em deg.
            C : float
                Distância entre centros em mm.
            casting_method : str
                Tipo de fundição utilizado no pinhão.
                Opções: molde; resfriamento; centrifugo. O padrão é 'molde'
            **threadkwargs: kwargs
                Dados de adendo (kw: a) e dedendo (kw: b) e folga (kw: c) do 
                dente do pinhão. Somente nescessário se px < 4.06 mm. Por 
                padrão adota os valores padronizados para um angulo de pressão 
                de 20° de filente único.
    
            Returns
            -------
            None.
    
        """
        
        self.dw,self.dc,self.px,self.C = map(lambda x : m2ft(x,metric_unit='mm',imperial_unit='in'),
                                             [dw,dc,px,(dw+dc)/2]
                                             )
        
        self.nw = rad_velocity(nw,input_unit='RPM',output_unit='rad/s')
        self.lambd=np.deg2rad(lambd)
        
        self.Nw = Nw
        self.Nc = Nc
        self.mG = Nc/Nw 
        
        self.Pn = self.Nc/self.dc
        self.Pt = self.Pn*np.cos(self.lambd)
        
        a = threadkwargs.pop('a',0.3183*self.px)
        b = threadkwargs.pop('b',0.3683*self.px)
        c = threadkwargs.pop('c',b-a)

        do = self.dw+2*a
        
        "Largura de Face da Coroa"
        self.Fc = (1.125*np.sqrt((do+2*c)**2-(do-4*a)**2) if px < 0.160
                   else 2/3*self.dw)
        #Não exceder valor máximo:
        if self.Fc > 0.67*self.dc:
            self.Fc = 0.67*self.dc
        
        "Velocidade tangencial"
        self.Vs = 60*(np.pi*self.nw*self.dw/(12*np.cos(self.lambd))) #ft/min
        
        "Fator dos Materiais"
        if self.dc <= 3:
            Cs = 720+10.37*self.C**3
        elif casting_method.lower() == 'molde' or 'mold':
            Cs = (1000 if self.dc <= 2.5
                  else 1190 - 477*np.log10(self.dc)
                  )
        elif casting_method.lower() == 'resfriamento' or 'cooled':
            Cs = (1000 if self.dc <= 8
                  else 1412 - 456*np.log10(self.dc)
                  )
        elif casting_method.replace('í','i').lower() == 'centrifugo' or 'centrifugal':
            Cs = (1000 if self.dc <= 25
                  else 1251 - 180*np.log10(self.dc)
                  )
            
        self.Cs = Cs
        
        "Fator de correção de razão"
        if self.mG < 3:
            print('A razão de engrenamento Nc/Nw deve ser maior que 3. Recebido: {mG:.1f}')
        
        
        self.Cm = (0.02*np.sqrt(-self.mG**2+40*self.mG-76) + 0.46 if 3 < self.mG <= 20
                   else 0.0107*np.sqrt(-self.mG**2+56*self.mG+5145) if 20 < self.mG <= 76
                   else 1.1483-0.00658*self.mG)
        
        "Fator de velocidade"
        self.Cv = (0.659*np.exp(-0.0011*self.Vs) if self.Vs < 700
                   else 13.31*self.Vs**(-0.571) if 700 <= self.Vs < 3000
                   else 65.52*self.Vs**(-0.774)
                   )
        
    
    def Wadm(self,E=303,unit='N',FS=1.5):
        """
        Calcula a força tangencial admissível no dente da coroa (componente mais crítico)
        a partir do método AGMA, corrigido pela razão entre os módulo de elasticidade
        do material adotado e a liga de bronze adotada pela AGMA para o dimensionamento
        de ssitemas de engrenamento sem-fim.

        Parameters
        ----------
        
        E : float
            Módulo de Young do material em MPa. O padrão é 303, idêntico à liga
            de aluínio-bronze considerada para os cálculos padrões da AGMA.
        unit : str, optional
            Define a unidade desejada da força no output. The default is 'N'.
        FS : float, optional
            Fator de segurança para contabilizar parâmetros não previstos nos 
            efeitos da mudança de material. The default is 1.5.

        Returns
        -------
        float
            Força tangecial admissível na unidade específicada.

        """
        E_bronze = 303 #Mpa, modulo de elasticidade de uma liga de bronze-almínimo, utilizada no calculo da AGMA
        E_ratio = E/(E_bronze*FS)
        
        W = self.Cs*(self.dc**0.8)*self.Fc*self.Cm*self.Cv*E_ratio
        
        

        if unit == 'lbf':
            return W
        
        else:
            return lbf2N(W,metric_unit=unit)
        
        
    def get_params(self):
        """
        Gera dicionário com os parâmetros de correção do método AGMA.

        Returns
        -------
        params_dict : dict
            Dicionário com os parâmetros de correção do método AGMA.

        """
        
        params_dict = {'C_s':self.Cs,
                       'C_m':self.Cm,
                       'C_v':self.Cv,
                       'Vs':self.Vs
                       }
        
        
        return params_dict
    
    def get_dims(self):
        """
        Gera dicionário com os parâmetros dimensionais do engrenamento

        Returns
        -------
        params_dict : dict
            Dicionário com os parâmetros dimensionais do engrenamento.

        """
        params_dict = {'a':self.a,
                       'b':self.b,
                       'c':self.c,
                       'lambda':self.lambd,
                       'Pn':self.Pn,
                       'Pt':self.Pt}
        
        
        return params_dict

        

