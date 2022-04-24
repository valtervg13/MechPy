# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 07:05:14 2022

@author: valter
"""
import numpy as np

prefix = {'m':1e-3,
          'c':1e-2,
          'd':1e-1,
          '':1,
          'da':1e1,
          'h':1e2,
          'k':1e3,
          'M':1e6,
          'G':1e9
          }

imperial = {'in':1/12,
            'ft':'',
            'yd':'3'}

def SI(value:float,
       input_unit:str,
       output_unit:str):
    """
    Converte unidades do SI

    Parameters
    ----------
    value : float
        Valor a ser convertido.
    input_unit : str
        Unidade da entrada. Ex: mm2, GPa, cm3.
    output_unit : str
        Unidade da saída.

    Returns
    -------
    float
        Valor convertido da unidade de entrada para a unidade de saída.

    """
    if input_unit=='mm' or output_unit=='mm':
        si_unit = 'm'
        
    else:
        si_unit=''
        for char in input_unit:
            if char in output_unit:
                si_unit+=char
    
    if input_unit != si_unit:
        
        if 'mm' in input_unit:
            input_split = input_unit.split('mm')
            input_prefix = 'm'
            input_power = (1 if input_split[1] == ''
                           else float(input_split[1])
                           )
        else:
            input_split = input_unit.split(si_unit)
           
            if len(input_split) < 2:
                input_prefix = input_split[0]
                input_power = 1
            else:
                print(input_split)
                input_prefix = input_split[0]
                input_power = (1 if input_split[1] == ''
                               else float(input_split[1])
                               )
       
        input_rate = prefix[input_prefix]**input_power
    
    else:
        input_rate = 1
    
    
    if output_unit != si_unit:
        
        if 'mm' in output_unit:
            output_split = output_unit.split('mm')
            output_prefix = 'm'
            output_power = (1 if output_split[1] == ''
                           else float(output_split[1])
                           )
         
        else:
            output_split = output_unit.split(si_unit)
           
            if len(output_split) < 2:
                output_prefix = output_split[0]
                output_power = 1
            else:
                output_prefix = output_split[0]
                output_power = (1 if output_split[1] == ''
                               else float(output_split[1])
                               )
           
        output_rate = prefix[output_prefix]**output_power
           

    
    else:
        output_rate = 1
        
    
    rate = input_rate/output_rate
    
    return rate*value


def rad_velocity(n:float,input_unit:str,output_unit:str):
    """
    Converte velocidades radiais de RPM para rad/s e vice versa

    Parameters
    ----------
    n : float
        velocidade radial.
    input_unit : str
        unidade de entrada, pode ser 'RPM' ou 'rad/s'.
    output_unit : str
        unidade de saída, pode ser 'RPM' ou 'rad/s'.

    Returns
    -------
    float
        velocidade radial convertida.

    """
    
    if input_unit == output_unit:
        return n
    
    elif input_unit.lower() == "rpm" and output_unit=='rad/s':
        omega = n*2*np.pi/60
    

    elif input_unit == "rad/s" and output_unit.lower() == "rpm":
        omega = 60*n/(2*np.pi)
    
    return omega

def m2ft(x:float,metric_unit:str='m', imperial_unit:str='ft'):
    """
    Converte comprimento métrico para pés.

    Parameters
    ----------
    x : float
        valor a ser convertido.
    metric_unit : str, optional
        Unidade métrica adotada. permite utilizar valores diretamente em 
        unidades como mm, dm, km e outras. The default is 'm'.

    Returns
    -------
    float
        Valor convertido para unidade imperial especificada.

    """
    x = SI(x,metric_unit,'m')
    imp_rate = 3.28084/imperial[imperial_unit]
    
    return x * imp_rate

def ft2m(x:float, metric_unit:str='m', imperial_unit:str='ft'):
    """
    Converte pés para comprimento métrico.

    Parameters
    ----------
    x : float
        valor a ser convertido.
    metric_unit : str, optional
        Unidade métrica adotada. permite utilizar valores diretamente em 
        unidades como mm, dm, km e outras. The default is 'm'.

    Returns
    -------
    float
        Valor convertido para unidade métrica específicada.

    """
    imp_rate = 3.28084*imperial[imperial_unit] 
    
    m = x/imp_rate
    
    return SI(m,'m',metric_unit)

def N2lbf(x:float, metric_unit:str='N'):
    
    x = SI(x,metric_unit,'N')
    imp_rate = 1/4.448
    
    return x*imp_rate

def lbf2N(x:float, metric_unit:str='N'):
    
    imp_rate = 1/4.448
    
    F = x/imp_rate
    
    return SI(F,'N',metric_unit)


    