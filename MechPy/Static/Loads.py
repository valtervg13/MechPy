# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 16:06:42 2022

@author: valte
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle as pc
import copy as cp
from scipy.optimize import fsolve

def bmatrix(a,precision):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = np.array2string(a, 
                            max_line_width=np.infty, 
                            precision=precision).replace('[', '').replace(']', '').splitlines()
    
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)


class Loads():
    def __init__(self):
        pass
    
    class Load():
        """

        Parameters
        ----------
        tag : str
            Nome da variável da carga.
        position : array_like.
            Vetro posição em coordenadas cartesianas (x,y,z).
        value : array-like ou str, optional
            Vetor de intensidade da carga. Se `value` contiver uma string,
            considera-se que a componente é uma variável desconhecida. se
            for `None`, todas as componentes são consideradas variáveis. 
            The default is `None`.

        Returns
        -------
        None.

        """
        def __init__(self,tag, position, value=None, load_type = 'F'):

            self.tag = tag
            self.position = np.array(position)
            self.load_type = load_type
            self.isknown = ([True if type(f_i) != str
                             else False
                             for f_i in value
                             ] if value != None
                            else (False,False,False)
                            )
            
            if value == None:
                self.value = np.ones(3)
            else:
                self.value = np.array([f_i if type(f_i) != str
                                       else -1 if '-' in f_i
                                       else 1
                                       for f_i in value
                                       ]
                                      )
            

            self.value_x = self.value[0]
            self.value_y = self.value[1]
            self.value_z = self.value[2]
            
        def set_a(self,a_list):
            """
            Parameters
            ----------
            a_list : list or tuple
                Lista contendo os coeficientes escalares a serem 
                multiplicados pela magnitude de forma a gerar os momentos 
                causados pela carga. 
                M_x = a_list[0]*F, M_y = a_list[1]*F, M_z = a_list[2]*F.

            Returns
            -------
            None.

            """
            self.a_x = a_list[0]
            self.a_y = a_list[1]
            self.a_z = a_list[2]
            
        def set_value(self,val_list):
            """
            

            Parameters
            ----------
            val_list : list or tuple
                Lista contendo os valores de magnitude na forma (). 
                M_x = a_list[0]*F, M_y = a_list[1]*F, M_z = a_list[2]*F.

            Returns
            -------
            None.

            """
            self.isknown = ([True if type(f_i) != str
                             else False
                             for f_i in val_list
                             ])
            
            self.value = np.array([f_i if type(f_i) != str
                                   else -1 if '-' in f_i
                                   else 1
                                   for f_i in val_list
                                   ]
                                  )
            
            self.value_x = val_list[0]
            self.value_y = val_list[1]
            self.value_z = val_list[2]
            
            
        def get_values(self,axis=''):
            axis_dict = {'x':0,
                         'y':1,
                         'z':2
                         }
            
            idx = []
            for ax in axis:
                if ax not in idx:
                    idx.append(axis_dict[ax])
                    
            values = [cp.copy(self.value[ax]) for ax in idx]
            
            if len(values) == 1:
                return values[0]
            
            else:
                return values
        
        def get_pos(self,axis=''):
            axis_dict = {'x':0,
                         'y':1,
                         'z':2
                         }
            
            idx = []
            for ax in axis:
                if ax not in idx:
                    idx.append(axis_dict[ax])
                    
            values = [cp.copy(self.position[ax]) for ax in idx]
            
            if len(values) == 1:
                return values[0]
            
            else:
                return values
            
            
            
            
                
    def solve_reactions(*loads:Load,
                        load_init = None,
                        coords='ijk',
                        add_equations=[],
                        print_matrix=False,
                        update_loads=False,
                        **kwargs):
        """
        Parameters
        ----------
        loads: Load
            Carga, conforme definida na classe Load
        coords : TYPE, optional
            Indica a orientação das coordenadas de origem através da notação
            de Einstein, de modo que Se `coords` = ijk, temos, para a1 = 1e_x e a2 = 1e_y
            a1 x a2 = a1*a2 e_z. Já para ikj têm-se a1 x a2 = -a1*a2 e_z 
            O padrão é 'ijk'.

        Returns
        -------
        Dicionário contendo as reações. Chaves são definidas pela tag de cada 
        objeto "Load" acrescidos de um sufixo "_x","_y","_z" a depender da
        componente.
        """
        
        F_loads = [load for load in loads if load.load_type.lower()=='f']
        M_loads = [load for load in loads if load.load_type.lower()=='m']
        
        #========================================================================
        # MATRIZES DO SISTEMA
        #========================================================================
        
        F = [[],
             [],
             []]
        
        Vars = []
        
        R = [[],
             [],
             []]
        
        M = [[],
             [],
             []]
        
        M_r = [[],
               [],
               []]
        
        #========================================================================
        # OBTENÇÃO DAS DIREÇÕES E SENTIDOS DOS MOMENTOS
        #========================================================================
        for load in loads:
            subs = (f'{coords[0]}{coords[1]}{coords[2]},{coords[1]},{coords[2]}->{coords[0]}'
                      )
            
            levi_civita = np.array([[[int((i - j) * (j - k) * (k - i) / 2)
                                      for k in range(3)
                                      ]
                                     for j in range(3)
                                     ]
                                    for i in range(3)
                                    ]
                                   )
            


            x = load.position
            
            origin = np.ones(3)
            i = np.array([1,0,0])
            j = np.array([0,1,0])
            k = np.array([0,0,1])
            
            
            a_i = np.einsum(subs,levi_civita,x,i)
            a_j = np.einsum(subs,levi_civita,x,j)
            a_k = np.einsum(subs,levi_civita,x,k)
            
            a_list = [a_i,
                      a_j,
                      a_k]
            load.set_a(a_list)
            
            tag_sufix = ['_x','_y','_z']
            for n in range(3):
                if load.load_type.lower() == 'f':
                    if load.isknown[n]:
                        for m in range(3):
                            if m==n:
                                F[m].append(load.value[n])
                            else:
                                F[m].append(0)
                    
                            M[m].append(a_list[n][m]*load.value[n])
                    else:
                        for m in range(3):
                            if m==n:
                                Vars.append(load.tag+tag_sufix[n])
                                R[m].append(load.value[n])

                            else:
                                R[m].append(0)  
                            
                            M_r[m].append(a_list[n][m]*load.value[n])
                            
                if load.load_type.lower() == 'm':
                    if load.isknown[n]:
                        for m in range(3):
                            if m==n:
                                F[m].append(0)
                                M[m].append(load.value[n])
                            else:
                                F[m].append(0)
                                M[m].append(0)
                    else:
                        for m in range(3):
                            if m==n:
                                Vars.append(load.tag+tag_sufix[n])
                                R[m].append(0)
                                M_r[m].append(load.value[n])
                            else:
                                M_r[m].append(0)
                                R[m].append(0)
        
        var_len = len(Vars)
        
        #======================================================================
        # MATRIZ DE FUNÇÕES
        #======================================================================
        
        #Checagem da dependencia linear
        
        exp_matrix = [[] for n in range(6)]
        
        for n in range(3):
            exp_matrix[n] = [*R[n],sum(F[n])]
            exp_matrix[n+3] = [*M_r[n],sum(M[n])]
        
        exp_matrix = np.array(exp_matrix)
        
        
        lin_dep = []
        for n in range(6):
            if all([exp_matrix[n,k]==0 
                    for k in range(exp_matrix.shape[1])
                    ]
                   ):
                lin_dep.append(n)
                
            else:
                for m in range(6):
                    if n != m:
                        inner_product = np.inner(exp_matrix[n,:],
                                                 exp_matrix[m,:]
                                                 )
                        norm_n = np.linalg.norm(exp_matrix[n,:])
                        norm_m = np.linalg.norm(exp_matrix[m,:])
            
                        if np.abs(np.abs(inner_product) - np.abs(norm_n * norm_m)) < 1E-4: #inequação de Cauchy
                            lin_dep.append(m)
        
        
        R = [R[n] 
             for n in range(3)
             if n not in lin_dep 
             ]
        
        F = [F[n] 
             for n in range(3)
             if n not in lin_dep 
             ]
        
        
        M_r = [M_r[n] 
               for n in range(3)
               if n+3 not in lin_dep 
               ]
        
        M = [M[n] 
             for n in range(3)
             if n+3 not in lin_dep 
             ]
        
        
        def equations(R_0):
            
            R_eval = np.matmul(np.array(R),R_0)
            M_reval = np.matmul(np.array(M_r),R_0)
            
            
            F_eq = [R_eval[n]+sum(F[n]) for n in range(len(R_eval))]
            M_eq = [M_reval[n]+sum(M[n]) for n in range(len(M_reval))]
            
            
            add_eval = [func(R_0) for func in add_equations]
            
            eq_list = [*F_eq,
                       *M_eq]
                
            if add_eval != []:
                for eq_eval in add_eval: 
                    eq_list.append(eq_eval)
                
            return eq_list
        
        R_0 = (np.ones(var_len) 
               if load_init == None
               else load_init
               )
        if print_matrix:
            
            precision = kwargs.pop('precision',2)
                        
            print(f"""
                  
Matriz do sistema:
    
{bmatrix(np.array([*R,*M_r]),
         precision=precision)
 }{bmatrix(np.transpose([Vars]),precision=precision)
   } = {bmatrix(np.transpose([[*np.sum(F,axis=1),*np.sum(M,axis=1)]]),
                precision=precision)}
                  
""")
           
        Roots = fsolve(equations,R_0)

        
        root_dict = dict([(Vars[n],Roots[n]) 
                          for n in range(var_len)
                          ]
                         )
        
        if update_loads == True:
            for root in list(root_dict.keys()):
                axis = {'x':0,
                        'y':1,
                        'z':2}
                        
                tag,ax = root.split('_')
                idx = axis[ax]
    
                update_val = [0,0,0]
                
                for load in loads:
                    if tag == load.tag: 
                        for i in range(3):
                            if i==idx:
                                update_val[i] = root_dict[f'{tag}_{ax}']
                            else:
                                update_val[i] = load.value[i]

                        load.set_value(update_val)

                
                

        
        return root_dict
    
    def solve_internal(*loads:Load,
                       axis='x',
                       pos_shift = (0,0),
                       l_divs=500,
                       **graph_kwargs
                       ):
        """
        

        Parameters
        ----------
        *loads : Load
            Cargas externas.
        axis : TYPE, optional
            DESCRIPTION. The default is 'x'.
        pos_shift : TYPE, optional
            DESCRIPTION. The default is (0,0).
        #axis : TYPE, optional
            axis=x -> pos=(y,z), axis=y -> pos=(z,x), axis=z -> pos=(x,y). The default is (0,0).
        axis : TYPE, optional
            DESCRIPTION. The default is y -> pos=(z,x).
        axis : TYPE, optional
            DESCRIPTION. The default is z -> pos=(x,y)                       l_divs=500.
        **graph_kwargs : keyword arguments
            Argumentos adicionais para criação de gráficos e impressão de resultados.
                -F_graphs: str
                 string contendo "x", "y" e/ou "z". Gráficos dos esforços internos serão plotados.
                -M_graphs: str
                 string contendo "x", "y" e/ou "z". Gráficos dos momentos internos serão plotados.
                -print_reactions: bool
                Imprime no console as reações em cada ponto de aplicação das forças
                -pckl: str
                Salva os gráficos como arquivos editáveis. pckl deve conter o caminho com nome e
                extensão do arquivo

        Returns
        -------
        dict
            Dicionário contendo as forças internas em cada seção.
        dict
            Dicionário contendo os momentos internos em cada seção
        list
            Lista de indices das repostas internas nos pontos de aplicação das forças

        """
        
        axis_dict = {'x':0,
                     'y':1,
                     'z':2}
        
        ax = axis_dict[axis]
        
        n_load = len(loads)
        
        
        sort_loads = sorted(loads,key = lambda load: load.position[ax])
        l_section = [load.position[ax] for load in sort_loads]
        
        l_range = []
        for k in range(n_load):

            if l_section[k] != l_section[-1]:
                l_inner =  np.linspace(l_section[k],
                                       l_section[k+1],
                                       l_divs
                                       )
                for l in l_inner:
                    l_range.append(l)
            else:
                l_range.append(l_section[k])
                

        def remove_dups(List):
            
            return list(map(float,
                            list(dict.fromkeys(map(str,
                                                   List)
                                               )
                                     )
                            )
                        )
        
        l_range = np.array(remove_dups(l_range))
        
        
        F_roots = {axis: [],
                   'R_x': [],
                   'R_y': [],
                   'R_z': []
                   }
        M_roots = {axis: [],
                   'Mr_x': [],
                   'Mr_y': [],
                   'Mr_z': []
                   }
        
        for k in range(len(l_range)):
            
            pos = [0,0,0]
            pos[ax] = l_range[k]
            if ax == 0:
                pos[1] = pos_shift[0]
                pos[2] = pos_shift[1]
            elif ax == 1:
                pos[2] = pos_shift[0]
                pos[0] = pos_shift[1]
            else:
                pos[0] = pos_shift[0]
                pos[1] = pos_shift[1]

            
            R_int = Loads.Load(f'R{k}',
                               position=pos)
          
            M_int = Loads.Load(f'Mr{k}',
                               position=pos,
                               load_type='M')
            
            load_list = [*[load 
                           for load in loads 
                           if load.position[ax] <= l_range[k]
                           ],
                         R_int,
                         M_int
                         ]

            node_solution = Loads.solve_reactions(*load_list)
            
            F_roots[axis].append(l_range[k])
            F_roots['R_x'].append(node_solution[f'R{k}_x'])
            F_roots['R_y'].append(node_solution[f'R{k}_y'])
            F_roots['R_z'].append(node_solution[f'R{k}_z'])

            
            M_roots[axis].append(l_range[k])
            M_roots['Mr_x'].append(node_solution[f'Mr{k}_x'])
            M_roots['Mr_y'].append(node_solution[f'Mr{k}_y'])
            M_roots['Mr_z'].append(node_solution[f'Mr{k}_z'])
            
        #==================================================================
        # PLOTAGEM DE DIAGRAMAS
        #==================================================================
        
        F_graphs = graph_kwargs.pop('F_graphs','')
        units = graph_kwargs.pop('units',{'distance':'m',
                                          'force':'N',
                                          'torque':'Nm'})
        
        pckl_dict = graph_kwargs.pop('pckl',None)
        
        if pckl_dict is not None:
            pckl_title = pckl_dict['title']
            pckl_path = pckl_dict['path']
            pckl_tag = pckl_dict.pop('tag','')
            
            

        F_idx = []

        for graph_axis in ['x','y','z']:
            graph_label = ('$N$' 
                           if graph_axis == axis
                           else f'$V_{graph_axis}$'
                           )
            
            x_label = f'{graph_axis} [{units["distance"]}]'
            y_label = (f'N: Esforço normal [{units["force"]}]' 
                       if graph_axis == axis
                       else f'$V_{graph_axis}$: Esforço Cortante em {graph_axis}'
                       )
            
            l = np.array(F_roots[axis])
            F = np.array(F_roots[f'R_{graph_axis}'])
            
            l_ticks = l_section
            f_indexes = np.where([l_value in l_ticks for l_value in l])[0]
            f_ticks = [F[i] for i in f_indexes]
            
            F_idx.append(f_indexes)

            if graph_axis in F_graphs:                

                fig = plt.figure(figsize=(9,6))
                plt.rc('font',**{'family':'serif','serif':['Times New Roman']})
                ax = fig.gca()
                ax.plot(l,F,'b-.',label=graph_label)
                ax.legend([graph_label])
                plt.fill_between(l,F,color='blue',**{'alpha':0.2})
                ax.set_xlim(0,np.max(l)+0.005*np.max(l))
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                plt.xticks(l_ticks,[f'{tick:.2f}' for tick in l_ticks])
                plt.yticks(f_ticks,[f'{tick:.2f}' for tick in f_ticks])
                plt.grid()
                
                if pckl_dict is None:
                    plt.show()
                    
                else:
                    pckl_file = f'{pckl_path}/{pckl_title}_{graph_label}_{pckl_tag}.txt'
                    
                    
                    with open(pckl_file,'wb') as file:
                        pc.dump(fig,file)
                
                
        M_graphs = graph_kwargs.pop('M_graphs','')
        
        M_idx = []

        for graph_axis in ['x','y','z']:
            graph_label = ('$T$' 
                           if graph_axis == axis
                           else f'$Mf_{graph_axis}$'
                           )
            
            x_label = f'{graph_axis} [{units["distance"]}]'
            y_label = (f'T: Torque [{units["torque"]}]' 
                       if graph_axis == axis
                       else f'$Mf_{graph_axis}$: Momento fletor em {graph_axis}'
                       )
            
            l = np.array(M_roots[axis])
            M = np.array(M_roots[f'Mr_{graph_axis}'])
            
            l_ticks = l_section
            m_indexes = np.where([l_value in l_ticks for l_value in l])[0]
            m_ticks = [M[i] for i in m_indexes]
            
            M_idx.append(m_indexes)
            
            if graph_axis in  M_graphs:
                    
                    fig = plt.figure(figsize=(9,6))
                    plt.rc('font',**{'family':'serif','serif':['Times New Roman']})
                    ax = fig.gca()
                    ax.plot(l,M,'b-.',label=graph_label)
                    ax.legend([graph_label])
                    plt.fill_between(l,M,color='blue',**{'alpha':0.2})
                    ax.set_xlim(0,np.max(l)+0.005*np.max(l))
                    ax.set_xlabel(x_label)
                    ax.set_ylabel(y_label)
                    plt.xticks(l_ticks,[f'{tick:.2f}' for tick in l_ticks])
                    plt.yticks(m_ticks,[f'{tick:.2f}' for tick in m_ticks])
                    plt.grid()
                    
                    if pckl_dict is None:
                        plt.show()
                        
                    else:
                        pckl_file = f'{pckl_path}/{pckl_title}_{graph_label}_{pckl_tag}.txt'
                        
                        with open(pckl_file,'wb') as file:
                            pc.dump(fig,file)
                
        
        indexes = {'F_idx': list(map(int,F_idx[0])),
                   'M_idx': list(map(int,M_idx[0]))
                   }
        
        
        print_reactions = graph_kwargs.pop('print_reactions',False)
        if print_reactions == True:
            print('\nEsforços Internos')
            for R in ['R_x','R_y','R_z']:
                if R.split('_')[1] in F_graphs:
                    print('\n')
                    for i in indexes['F_idx']:
                        x = F_roots[axis][i]
                        
                        if i!=0:
                            print(f'{R}({axis}={x} "->") = {F_roots[R][i-1]:.3f} {units["force"]}')
                            
                        print(f'{R}({axis}={x} "<-") = {F_roots[R][i]:.3f} {units["force"]}')

            for Mr in ['Mr_x','Mr_y','Mr_z']:
                if Mr.split('_')[1] in M_graphs:
                    print('\n')
                    for i in indexes['M_idx']:
                        x = M_roots[axis][i]
                        
                        if i!=0:
                            print(f'{Mr}({axis}={x} "->") = {M_roots[Mr][i-1]:.3f} {units["torque"]}')
                            
                        print(f'{Mr}({axis}={x} "<-") = {M_roots[Mr][i]:.3f} {units["torque"]}')

        
        return F_roots,M_roots,indexes
            
