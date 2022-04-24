# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:23:38 2022

@author: valter

"""

from mechpy import Stress as st


circle = st.Section('circle',30e-3)
print(circle.properties())


load_dict = {'N': 5000,
             'Mf_z': 10,
             'T': 50}

S = st.get_stress(circle, (0,15e-3),normal='y',**load_dict)

print(S)

sig_vm = st.VonMisses(S)
print(sig_vm*1e-6)

