# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 15:48:01 2021

@author: valte
"""
import mechpy as mc

s = mc.gears.GearTrain.Component(1,drives = 3, first_gear = True, n = -100, N = 20)
p = mc.gears.GearTrain.Component(2,driven_by=1,drives=3,connects=5,N = 30)
a = mc.gears.GearTrain.Component(3,orient="internal",driven_by = 2,N = 80)
b = mc.gears.GearTrain.Component(4,Type="arm",connects = 2)

p2 = mc.gears.GearTrain.Component(5,drives=6,N = 30)
s2 = mc.gears.GearTrain.Component(6,driven_by=5,N = 60)
train = mc.gears.GearTrain([s,p,a,b,p2,s2])

print(train.e_total)
print(train.e)
data = train.n_data_matrix
print(data)
