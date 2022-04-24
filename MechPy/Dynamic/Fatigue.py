# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:57:20 2022

@author: valte
"""
import numpy as np
import pandas as pd


class Fatigue():
    def __init__(self):
        return
    def resitencia_fadiga(Sut,
                          Se_cfg,
                          Load,
                          SecType,
                          SecDim,
                          Finish,
                          Trust,
                          Temp=20,
                          Other=1,
                          print_tables = True,
                          full_output = False
                          ):
        """
        

        Parameters
        ----------
        Sut : float
            Limite de resitência à tração em MPa.
        Se_cfg : float ou string
            Configura o valor de Se'. Caso um float seja entrado, representa 
            o valor de Se' em MPa. Caso seja uma string, representa o material
            utilizado para cálculo de Se' a partir de Sut. Se a string for 'aço',
            utiliza-se os valores para o aço. Caso contrário, assume-se alumínio
        Load : str ou int
            Tipo de carregamento sofrido pelo componente. opções:
             -0/'flexao': cargas flexionais
             -1/'axial': cargas axiais
             -3/'torção': cargas torcionais
             -4/'misto': comninação de cargas
        SecType : str ou int
            Formato da seção transversal. opções:
                 -0/'retangular'
                 -1/'circular': Seção circular sem rotação
                 -2/'circular rodando': seção circualr em rotação
        SecDim : tuple
            Dimensões da seção. Para retangulos: (a,b) com "a" largura e "b' a 
            largura. Para circulos (d) com "d" o diâmetro em mm.
        Finish : str ou int
            Acabamento superficial. opções:
                 -0/'retificado'
                 -1/'usinado'/'estirado a frio'
                 -2/'laminado a quente'
                 -3/'forjado'
        Trust : float
            Confiabilidade. selecionar entre 99.9999, 99.999, 99.99, 99.9, 99, 
            95, 90'.
        Temp : float, optional
            Temperatura de operação. The default is 20.
        Other : float, optional
            Fatores diversos. The default is 1.
        print_tables : bool, optional
            Imprime as tabelas contendo os valores dos modificadores. The default is True.
        full_output : bool, optional
            Retorna um dicionário com Se e os fatores de correção obtidos.
            The default is False.

        Returns
        -------
        Se: float
            Limie de resitência à fadiga corrigido.
            
        Se_dict: dict, optional (full_output = True)
            Dicionário contendo o valor de Se, Se' e dos fatores de correção

        """
        
        Sut = Sut*1e6
        
        Load = str(Load)
        
        SecDim = tuple(map(lambda x: 1e-3*x,SecDim))
        
        if Se_cfg == "aço":
            if Sut < 1400*10**6:
                Se_line = 0.504*Sut
            else:
                Se_line = 700*1e6
        
        elif type(Se_cfg) is int or float:
            Se_line = Se_cfg*1e6 
            
        else: #Assume-se alumínio
            if Sut < 400*1e6:
                Se_line = 0.4*Sut
            else:
                Se_line = 160*1e6
                
        
        #C_Load (k_c)
        
        if Load.lower() in [0,'flexão','flexao',4,'misto']:
            Cload = 1
        elif Load.lower() in [1,'axial']:
            Cload = 0.85
        elif Load.lower() in (2,'torção','torçao','torcao','torcão'):
            Cload = 0.59

        
        
        #C_size (k_b)
        if Load.lower() in [1,'axial']:
                Csize = 1
                
        else:
            if SecType in [0,'retangular']:
                A95 = 0.05*SecDim[0]*SecDim[1]
                d = np.sqrt(A95/0.0766)
    
            elif SecType in [1,'circular']:
                A95 = 0.01046*(SecDim[0])**2
                d = np.sqrt(A95/0.0766)
            
            elif SecType in [2,'circular rodando','circular girando','eixo']:
                d = SecDim[0]
                
            if d < 2.79e-3 or d>254e-3:
                Csize = 1
            elif d<51e-3:
                Csize = 1.243*(d*10**3)**(-0.107)
            else:
                Csize = 1.51*(d*10**3)**(-0.157)
        
        #C_surface (k_a)
        if Finish  in [0,'retificado']:
            a=1.58
            b= -0.085
        elif Finish in [1,'usinado','estirado a frio']:
            a=4.51
            b= -0.265
        elif Finish in [2,'laminado a quente']:
            a=57.7
            b= -0.718
        elif Finish == [3,'forjado']:
            a=272
            b= -0.995
            

        Csurface = a*(Sut*1e-6)**(b)

        
        #C_trust (k_e)
        if Trust == 99.9999:
            Ctrust = 0.620
        elif Trust == 99.999:
            Ctrust = 0.659
        elif Trust == 99.99:
            Ctrust = 0.702
        elif Trust == 99.9:
            Ctrust = 0.753
        elif Trust == 99:
            Ctrust = 0.814
        elif Trust == 95:
            Ctrust =0.868
        elif Trust == 90:
            Ctrust =0.897
        else:
            raise ValueError('Enter a valid trust value from 99.9999, 99.999, 99.99, 99.9, 99, 95, 90')
        
        #C_temp(k_d)
        T_F = 9*Temp/5+32 #temperatura em farenheit
        if Temp == 20:
            Ctemp=1
        else: 
            Ctemp = 0.975+0.432e-3*T_F-0.115e-5*T_F**2+0.104e-8*T_F**3-0.595e-12*T_F**4
            
        
        #C_other (k_f)
        
        Cother = Other
        
        
        #Cálculo do Limite de Resistência à Fadiga
        S_e = Se_line*Cload*Csize*Csurface*Ctrust*Cother
        
        if print_tables == True:
            Cols = ['Fator de Correção','Valor','Variáveis','Valores']
            InfoNames = ['Tipo de Carregamento',
                          'Diametro Efetivo',
                          'Acabamento',
                          'a',
                          'b',
                          'Nível de Confiabilidade',
                          'Temperatura [°F]',
                          '-']
            InfoValues = [Load,
                          str('%.1f mm'%(d*10**3)),
                          Finish,
                          f"{a}",
                          f"{b}",
                          str('%s\%%'%(Trust)),
                          f'{T_F:.1f}',
                          "-"]
            Cnames= ['$C_{carregamento} (k_c)$',
                    '$C_{tamanho} (k_b)$',
                    '$C_{superficie} (k_a)$',
                    '{}',
                    '{}',
                    '$C_{confiabilidade} (k_e)$',
                    '$C_{Temperatura} (k_d)$',
                    '$C_{outros} (k_f)$}']
            
            Cvalues = [f'{Cload:.3f}',
                       f'{Csize:.3f}',
                       f'{Csurface:.3f}',
                       '{}',
                       '{}',
                       f'{Ctrust:.3f}',
                       f'{Ctemp:.3f}',
                       f'{Cother:.3f}']
            
            print('Tabela de Fatores de Correção: \n\n')
            pd.options.display.float_format = '{:.3f}'.format
            
            
            Cdata = pd.DataFrame(data=np.transpose([Cnames,
                                  Cvalues,
                                  InfoNames,
                                  InfoValues]),
                                 columns=[Cols])
            CdataLTX = Cdata.to_latex(index=False,decimal=',',escape=False,multicolumn=True)
            #print(CdataLTX)
            print(Cdata)
            
        if full_output == True:
            return  {'S_e' : S_e*1e-6, 'Se_line' : Se_line*1e-6, 'd': d,
                    'C_load' : Cload, 'C_size' : Csize, 'C_surface' : Csurface, 'C_trust' : Ctrust, 'C_temp': Ctemp,'C_other' : Cother}
        else:
            return S_e*1e-6
        
    def S_N(sigma,Sut,Type,Load,SecType,SecDim,Finish,Trust,
           print_tables = True,full_output = False,*args):

        if Load=="flexão":
            Sm = 0.9*Sut
        else:
            Sm = 0.75*Sut

        Sf_dic = Fatigue.resitencia_fadiga(Sut, Type, Load, SecType, SecDim, Finish, Trust, print_tables = print_tables, full_output=True)
        Se = Sf_dic['Se_line']
        Sf = Sf_dic['S_e']
        d = Sf_dic['d']
        Cload = Sf_dic['C_load']
        Csize = Sf_dic['C_size']
        Csurface = Sf_dic['C_surface']
        Ctrust = Sf_dic['C_trust']
        
        b = (1/(-3))*np.log10(Sm/Sf)
        a = 10**(np.log10(Sm)-3*b)
        
        if print_tables == True:
            StressDic ={}
            
            Snames = ["$S_e'$",'$S_m$','$S_e$','$a$','$b$']
            Svalues = ['%.2f'%(Se*10**(-6)),'%.2f'%(Sm*10**(-6)),'%.2f'%(Sf*10**(-6)),'%.2f'%((a*10**(-6))),'%.4f'%(b)]
            Sunit = ['Mpa','Mpa','Mpa','Mpa','-']
            
            print('Tabela de Tensões: \n\n')
            StressDic['Grandeza']=Snames
            StressDic['Valor']=Svalues
            StressDic['Unidade']=Sunit
            pd.options.display.float_format = '{:,.2f}'.format
            Sdata = pd.DataFrame(StressDic)
            SdataLTX = Sdata.to_latex(index=False,decimal=',',escape=False)
            print(SdataLTX)
            
        N = (sigma/a)**(1/b)
        
        S_e = Sf
        if full_output == True:
            return { 'N' : N, 'a' :a, 'b' : b, 'S_e' : S_e, 
                    'C_load' : Cload, 'C_size' : Csize, 'C_surface' : Csurface, 'C_trust' : Ctrust}
        
        else:
            return N
    
    
    def goodman_mod(S_ut,S_e,mode,**kwargs):
        if mode == 'n':
            sig_a = kwargs['sig_a']
            sig_m = kwargs['sig_m']
            
            n = 1/(sig_a/S_e+sig_m/S_ut)
            
            return n
        
        if mode == 'sig_a':
            n = kwargs['n']
            sig_m = kwargs['sig_m']
            
            sig_a = (1/n-sig_m/S_ut)*S_e
            
            return sig_a
            
        if mode == 'sig_m':
            n = kwargs['n']
            sig_a = kwargs['sig_a']
            
            sig_m = (1/n-sig_a/S_e)*S_ut
            
            return sig_m
    
    def soderberg(S_y,S_e,mode,**kwargs,):           
            return Fatigue.goodman_mod(S_y,S_e,mode = mode,**kwargs)
    
    def gerber(sig_a,sig_m,S_e,S_ut):            
            return (1/2)*((S_ut/sig_m)**2)*(sig_a/S_e)*(
                    (-1 + np.sqrt(1+((2*sig_m*S_e)/(S_ut*sig_a))**2))) #sig_m>0
        
    def ASME (sig_a,sig_m,S_e,S_y):
        return np.sqrt(1/((sig_a/S_e)**2+(sig_m/S_y)**2))
    
    def first_cycle(sig_a,sig_m,S_y, mode = 'n'):
        if mode == 'n':
            n = S_y/(sig_a + sig_m)
            
            return n
        
        elif mode == 'sig_a':
            sig_a = n/S_y -sig_m
            
            return sig_a
        
        elif mode == 'sig_m':
            sig_m = n/S_y-sig_a
            
            return sig_m
