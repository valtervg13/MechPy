# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 17:08:35 2022

@author: valte
"""
import numpy as np
from ..Dynamic.Fatigue import Fatigue as fatigue

def pitch_min_d(theta_max,L,E,F_H=[(0,0)],F_V=[(0,0)],M_H=[(0,0)],M_V=[(0,0)],FS=1,placement=1):
    """
    

    Parameters
    ----------
    theta_max : Float
        Máxima inclinação em radianos.
    L : FLoat
        Distância entre mancais em mm.
    E : Float
        Módulo de elasticidade em GPa.
    F_H : List of tuples, optional
        Em cada item, a primeira entrada é Força horizontal em N e a segunda é a 
        Distância entre Força e Mancal Esquerdo. The default is [(0,0)].
    F_V : List of tuples, optional
        Em cada item, a primeira entrada é Força verrtical  em N e a segunda é a 
        Distância entre Força e Mancal Esquerdo. The default is [(0,0)].
    M_H : List of tuples, optional
        Em cada item, a primeira entrada o Momento horizontal em Nm e a segunda é a 
        Distância entre o momento e Mancal Esquerdo. The default is [(0,0)].
    M_V : List of tuples, optional
        Em cada item, a primeira entrada é o Momento vertical em Nm e a segunda é a 
        Distância entre o momento e Mancal Esquerdo. The default is [(0,0)].
    FS : float, optional
        Fator de segurança. The default is 1.
    placement : int, optional
        Configuração das Cargas e Mancais
         -1: Cargas Centralizadas entre Mancais, inclui momentos centralizados
         -2: Cargas Sobressaltantes, a frente do mancal direito
         -0: Definido em cada carga, deve-se adiconar à F_V e F_H
          uma terceira entrada equivalente à configuração. 
          Ex: F_V = [(100,0.5,1),(200,1,2)]. The default is 1.

    Returns
    -------
    d_l : Float
        Diametro esquerdo em mm.
    d_r : Float
        Diametro direito em mm.

    """

    
    E = E*1e9 #Converte GPa em Pa
    
    def ab(array,plc=placement):
        if plc == 1:
            a_list = np.array([item[1] for item in array])
            b_list = np.array([L-a for a in a_list])
            
            a_list = a_list*1e-3
            b_list = b_list*1e-3
            return a_list,b_list
        
        elif plc == 2:
            a_list = np.array([item[1]-L for item in array])
            b_list = np.array([0 for a in a_list])
            
            a_list = a_list*1e-3
            b_list = b_list*1e-3
            return a_list,b_list
        
        elif plc == 0:
            a_list = []
            b_list = []
            
            for item in array:
                if array != [(0,0)]:

                    x = item[1]
                    plc = item[2]
                    
                    if plc == 1:
                        a_list.append(x)
                        b_list.append(L-x)
                        
                    elif plc == 2:
                        a_list.append(x-L)
                        b_list.append(0)
                else:
                    a_list.append(0)
                    b_list.append(0)
                    
            a_list = np.array(a_list)
            b_list = np.array(b_list)
            
            a_list = a_list*1e-3
            b_list = b_list*1e-3
            
            return a_list,b_list
                
    
    
    F_H_list = np.array([F[0] for F in F_H])
    a_FH, b_FH  = ab(F_H)

    
    F_V_list = np.array([F[0] for F in F_V])
    a_FV,b_FV = ab(F_V)
    
    M_H_list = np.array([M[0] for M in M_H])
    a_MH,b_MH = ab(M_H)
    
    M_V_list = np.array([M[0] for M in M_V])
    a_MV,b_MV = ab(M_V)

    H_sum_l = 0
    V_sum_l = 0
    
    H_sum_r = 0
    V_sum_r = 0
    
    #Somatorio de Forças em H
    L = L*1e-3 #Converte mm em m
    
    for i,(F,b,a) in enumerate(zip(F_H_list,b_FH,a_FH)):
        if placement == 0:
            plc = F_H[i][2]
        else:
            plc = placement
        if plc == 1:
            H_sum_l += F*b*(b**2-L**2)
            H_sum_r += F*a*(L**2-a**2)
        elif plc == 2:
            H_sum_l += F*a*L**2
            H_sum_r += F*a*3*L**2
    

    #Acréscimo dos Momentos
    H_sum_l += np.sum(M_H_list*(3*a_MH**2-6*a_MH*L+2*L**2))
    H_sum_r += np.sum(M_H_list*(3*a_MH**2-L**2))

    #Somatorio de Forças em V
    for i,(F,b,a) in enumerate(zip(F_V_list,b_FV,a_FV)):
        if placement == 0:
            plc = F_V[i][2]
        else:
            plc = placement
        if plc == 1:
            V_sum_l += F*b*(b**2-L**2)
            V_sum_r += F*a*(L**2-a**2)
        elif plc == 2:
            V_sum_l += F*a*L**2
            V_sum_r += F*a*3*L**2
            
    #Acréscimo dos Momentos
    V_sum_l += np.sum(M_V_list*(3*a_MV**2-6*a_MV*L+2*L**2))
    V_sum_r += np.sum(M_V_list*(3*a_MV**2-L**2))

    
    d_l = 1e3*(((32*FS)/
                (3*np.pi*E*L*theta_max)
                )*np.sqrt(H_sum_l**2+V_sum_l**2)
               )**0.25
    
    d_r = 1e3*(((32*FS)/
                (3*np.pi*E*L*theta_max)
                )*np.sqrt(H_sum_r**2+V_sum_r**2)
               )**0.25
        
    return d_l,d_r

def tensao_cisalhante(K_fs,T_a,T_m,d): #Shigley Ch. 7-4 eq. 7-4 
    tau_a = K_fs * (16*T_a)/(np.pi*d**3)
    tau_m = K_fs * (16*T_m)/(np.pi*d**3)
    return tau_a,tau_m

def tensao_normal(K_f,M_a,M_m,d): #Shigley Ch. 7-4 eq. 7-3
    sigma_a = K_f * (32*M_a)/(np.pi*d**3)
    sigma_m = K_f * (32*M_m)/(np.pi*d**3)
    return sigma_a,sigma_m

def tensao_vm(K_f,K_fs,M_a,M_m,T_a,T_m,d,full_output = True): #Shigley Ch. 7-4 eq. 7-5, 7-6
    tau_a,tau_m = tensao_cisalhante(K_fs,T_a,T_m,d)
    sig_a,sig_m = tensao_normal(K_f,M_a,M_m,d)
    
    sig_VMa = (sig_a**2+3*tau_a**2)**0.5
    sig_VMm = (sig_m**2+3*tau_m**2)**0.5
    
    if full_output == True:
        tension_state_a = [sig_a,tau_a]
        tension_state_m = [sig_m,tau_m]
        return sig_VMa,sig_VMm,tension_state_a,tension_state_m
    else:
        return sig_VMa,sig_VMm

def fatigue_nf(K_f,K_fs,T_a,T_m,M_a,M_m,d,S_y,S_ut,S_e,criteria = 'gerber',full_output = True): #Shigley Ch. 7-4 eq. 7-7 a 7-13 
    
    if full_output == True:
        sig_a,sig_m,tension_state_a,tension_state_m = tensao_vm(K_f,K_fs,M_a,M_m,T_a,T_m,d,full_output=True)
        
    else:
        sig_a,sig_m,tension_state_a,tension_state_m = tensao_vm(K_f,K_fs,M_a,M_m,T_a,T_m,d,full_output=True) 
    
    crit_dic = {"goodman" : fatigue.goodman_mod(S_ut = S_ut, S_e = S_e, mode ='n', sig_a = sig_a, sig_m = sig_m),
                "gerber" : fatigue.gerber(sig_a = sig_a, sig_m = sig_m, S_e = S_e, S_ut = S_ut),
                "asme" : fatigue.ASME(sig_a = sig_a, sig_m = sig_m, S_e = S_e, S_y = S_y),
                "soderberg" : fatigue.soderberg(S_y = S_y,S_e = S_e, mode='n', sig_a = sig_a, sig_m = sig_m)}
    
    sig_max = np.sqrt((tension_state_a[0]+tension_state_m[0])**2+3*(tension_state_a[1]+tension_state_m[1])**2)
    n_y = S_y/sig_max #Shigley Ch. 7-4 eq. 7-16 
    
    if full_output == True:
        return {"n_f":crit_dic.get(criteria.lower()),'n_yeld' : n_y,
                'sig_vm_a' : sig_a, 'sig_vm_m' : sig_m, 'sig_vm_max' : sig_max,
                'sig_a, tau_a' : tension_state_a, 'sig_m, tau_m' : tension_state_m}
    else:
        return crit_dic.get(criteria.lower()),n_y
    
def fatigue_d(K_f,K_fs,T_a,T_m,M_a,M_m,n,S_y,S_ut,S_e,criteria = 'gerber',full_output = True):
    
    crit_dic = {"goodman" : ((16*n/np.pi)*((1/S_e)*((4*(K_f*M_a)**2+3*(K_fs*T_a)**2))**0.5
                                    +(1/S_ut)*(4*(K_f*M_m)**2+3*(K_f*T_m)**2)**0.5))**(1/3),
                
                "gerber" : (((8*n*np.sqrt(4*(K_f*M_a)**2+3*(K_fs*T_a)**2))/(np.pi*S_e))*
                            ( 1 + (1+((2*np.sqrt(4*(K_f*M_m)**2+3*(K_fs*T_m)**2)*S_e)/
                                      (np.sqrt(4*(K_f*M_a)**2+3*(K_fs*T_a)**2)*S_ut)
                                      )**2
                                   )**0.5
                             )
                            )**(1/3),                    
                
                "asme" : ((16*n/np.pi)*(4*(K_f*M_a/S_e)**2+3*(K_fs*T_a/S_e)**2
                                        +4*(K_f*M_m/S_y)**2+3*(K_fs*T_m/S_y)**2)**0.5
                          )**(1/3),
                
                "soderberg" : ((16*n/np.pi)*((1/S_e)*(4*(K_f*M_a)**2+3*(K_fs*T_a)**2)**0.5
                                    +(1/S_y)*(4*(K_f*M_m)**2+3*(K_fs*T_m)**2)**0.5)
                               )**(1/3)
                }
    

    return crit_dic.get(criteria.lower())
 