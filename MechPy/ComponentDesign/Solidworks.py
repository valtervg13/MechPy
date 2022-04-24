# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 08:11:26 2022

@author: valte
"""

from numpy import *

def degsin(x):
    return sin(rad2deg(x))


def degcos(x):
    return cos(rad2deg(x))


def degtg(x):
    return tan(rad2deg(x))

def ARCTAN(x):
    return rad2deg(arctan(x))

def ARCSIN(x):
    return rad2deg(arcsin(x))


def ARCCOS(x):
    return rad2deg(arccos(x))

#Dicionário de equivalências entre funções do SW e da Numpy
funcs = {'sqr':'sqrt',
         'atn':'ARCTAN',
         'arcsin':'ARCSIN',
         'arccos':'ARCCOS',
         '^':'**',
         'cos':'degcos',
         'sin':'degsin',
         'tan':'degtan'}
    
def importDims(path,
               keep_coments=False,
               eval_eqs=True,
               **kwargs):
    """
    Importa dados de variáveis armazenados em arquivos externos de modelos 
    do Solidworks

    Parameters
    ----------
    path : str.
        Caminho do arquivo de equações do modelo do Solidworks.
    keep_coments : bool, optional
        Mantém os textos de descrição do arquivo. O padrão é False.
    eval_eqs : bool, optional
        Avalia o valor numérico de equações. Se falso, equações são retornadas 
        como strings. O padrão é True.
    **kwargs : keyword arguments.
        Argumentos adicionais.

    Returns
    -------
    dim_dict : dict.
        Dicionário contendo os valores dimensionais do arquivo do Solidworks.
        As chaves são definidas como os nomes das variáveis presentes no 
        arquivo

    """
    #Unidades presentes no arquivo de dimensões do SW
    #Nescessário para que sejam removidas ao realizar os cálculos
    
    default_units = ['mm', #comprimento
                    'deg','rad' #angulação
                    ]
    
    units = kwargs.pop('units',default_units)

    
    #Abre o arquivo de texto com as equações do SW
    with open(path,'r',encoding='utf-8-sig') as dim_file:
        
        
        dim_lines = dim_file.read().replace('\ufeff','').split('\n') #Divide por linhas, os três primeiros, caracteres possuem erros e são removidos
        dim_args = [line.split('=') 
                    for line in dim_lines if line != ''] #Separa nome da variável e equação, através do sinal de igual '=' 

        dim_args = [[key.replace(' ','').replace('"',''),
                     eq.split('\'')
                     ] 
                    if (keep_coments and len(eq.split('\''))==2)
                    else [key.replace(' ','').replace('"',''),
                          eq.split('\'')[0]
                          ] 
                    for key,eq in dim_args
                    ] #Formata como uma lista de listas no formato: [Nome da Variávle, Equação da Variavel]. Remove os espaços de tudo e aspas do nome da variável    
        

        
        dim_args_dict = dict(dim_args) #transforma em dicionário a lista de listas
        
        dim_dict = {} #Inicializa dicionário que armazenara os valores calculados

        
        #Loop que formata as equações em fomato compreensível pela Numpy
        if eval_eqs:
            for key,eq in dim_args_dict.items():
                
                comments_on = False
                comment = ''
                if type(eq) != str:
                    comments_on = True
                    
                    comment = eq[1]
                    eq = eq[0]
                
                eq = eq.replace(' ','')
                for unit in units:
                    eq=eq.replace(unit,'')
                for swfunc,func in funcs.items():
                    eq=eq.replace(swfunc,func)
                
                dim_args_dict[key] = ([eq,comment] if comments_on
                                      else eq)
        
        #Loop que substitui os nomes de variáveis por seus valores e calcula o valor das equações
            for key,eq in dim_args_dict.items(): 
                
                comments_on = False
                comment = ''
                if type(eq) != str:
                    comments_on = True
                    
                    comment = eq[1]
                    eq = eq[0]
                    
                
                #Substitui nomes de variáveis pelas respectivas equações até haverem só números
                while ('"' in eq) == True:
                    for var in dim_args_dict.keys():
                        if type(dim_args_dict[var]) != str:
                            var_value = str(dim_args_dict[var][0])
                        else: 
                            var_value = str(dim_args_dict[var])
                            
                        eq=eq.replace(f'"{var}"',var_value)
        
                
                val = eval(eq) #Cálculo do valor
                
                dim_dict[key] = ([val,comment] if comments_on 
                                 else val)#Salva resultado no dicionário
                
                
        else:
            dim_dict = dim_args_dict
    
            
    return dim_dict

def updateDims(path,vars_dict: dict,**kwargs):
    """
    

    Parameters
    ----------
    path : path : str.
        Caminho do arquivo de equações do modelo do Solidworks.
    vars_dict : dict
        Dicionário contendo os nomes das variáveis a serem exportadas como 
        chaves e tuplas (<Valor>,<Descrição>) ou apenas o <Valor> como os 
        valores.
    **kwargs : keyword arguments.
        Argumentos adicionais.

    Returns
    -------
    None.

    """
    
    #Unidades presentes no arquivo de dimensões do SW
    #Nescessário para que sejam removidas ao realizar os cálculos
    
    default_units = ['mm',
                     'deg',
                     'rad'
                     ]
    
    units = kwargs.pop('units',default_units)
    
    #Dicionário de equivalências entre funções do SW e da Numpy
    
    old_vars = importDims(path,
                          units=units,
                          keep_coments=True,
                          eval_eqs=False)
    
    for key,val in old_vars.items():
        
        if key in vars_dict.keys():
            new_val = vars_dict[key]
            update_val = ([str(new_val)] 
                          if type(new_val) not in (list,tuple,set)
                          else list(map(str,(new_val)))
                          )
        else:
            update_val = ([str(val)] 
                          if type(val) not in (list,tuple,set)
                          else list(map(str,(val)))
                          )
        
        old_vars[key] = "'".join(update_val)

    new_vars = "\n".join([f'"{key}"={value}' 
                          for key,value in old_vars.items()
                          ]
                         )
    
    
    with open(path,'w',encoding='utf-8-sig') as file:
        file.write(new_vars)
        
        


