# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:02:17 2022

@author: valter
"""

from mechpy import Loads

M_0 = Loads.Load(tag='M0',
                 position=(0,0,0),
                 value=(0,'?',0),
                 load_type='m'
                 )

F_1 = Loads.Load(tag='F1',
                 position=(0,0,0),
                 value=('','',0)
                 )

F_2 = Loads.Load(tag='F2',
                 position=(3,0,0),
                 value=(0,'',0)
                 )



F_3 = Loads.Load(tag='F3',
                 position=(2,0,0),
                 value=(5,2,0)
                 )

T1 = Loads.Load(tag='T1',
                 position=(0,0,0),
                 value=(5,0,0),
                 load_type='m'
                 )

T2 = Loads.Load(tag='T2',
                 position=(2,0,0),
                 value=(-5,0,0),
                 load_type='m'
                 )

M1 = Loads.Load(tag='M1',
                 position=(2,0,0),
                 value=(0,1,0),
                 load_type='m'
                 )


R = Loads.solve_reactions(M_0,F_1,F_2,F_3,T1,T2,M1,
                          coords='ijk')
for key,value in R.items():
    print(f'{key} = {value}')

F_1.set_value([R['F1_x'],
               R['F1_y'],
               F_1.value_z])

F_2.set_value([F_2.value_x,
               R['F2_y'],
               F_2.value_z])

R_in,M_in = Loads.solve_internal(F_1,F_2,F_3,T1,T2,M1,
                                 F_graphs='xyz',
                                 M_graphs='xyz')

