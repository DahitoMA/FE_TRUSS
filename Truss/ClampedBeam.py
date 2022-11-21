#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from InitTest import *
from Truss import FEfunctions as FE
import numpy as np
import copy
import math
import timeit
import random
from mpl_toolkits import mplot3d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib.pyplot import cm
# import cma # From https://github.com/CMA-ES/pycma
# import platypus # From https://github.com/Project-Platypus/Platypus
# from platypus import *
# from platypus.algorithms import NSGAII, NSGAIII
# from platypus.core import Problem
# from platypus.types import Real, Integer
# from platypus import Hypervolume, calculate


xMin = 0
xMax = 1
yMin = 0
yMax = 0.25 * (xMax - xMin) #1
nEls = 13 # Number of elements
nN = 8 # Number of nodes
Fa = 10000 # Force applied
### Young's moduli
Esteel = 200e6 # steel
Ealu = 69e6 # aluminum
Emag = 45e6 # magnesium
Etitan = 105e6 # titanium
### Young's moduli *10^-6
Esteel2 = 200
Ealu2 = 69
Emag2 = 45
Etitan2 = 105
# Categorical variables, 1: titan   2: mag   3: steel   4: alu
Youngs = [Etitan, Emag, Esteel, Ealu]
Categ = [1,2,3,4] 
Youngs2 = [Etitan2, Emag2, Esteel2, Ealu2]

### Mass costs from: https://www.mineralinfo.fr
Csteel = 0.5 # 0,60€ / kg
Calu = 1.46 # 1,66€ / kg
Cmag = 1.88 # 1,88€ / kg
Ctitan = 3.6 # 3,60€ / kg
### Densities from https://fr.wikipedia.org/wiki/Masse_volumique
rhosteel = 7.85e3 # kg / m^3
rhoalu = 2.7e3 # kg / m^3
rhomag = 1.75e3 # kg / m^3
rhotitan = 4.5e3 # kg / m^3

# Constraints 
Uy3max =  6.99772047e-02 #case 1
Ux5max = -0.05797101 #case 2
    
" Build connectivity and global coordinates matrices "
connect = np.zeros((2, nEls))  # Connectivity matrix
for k in range(0,7):
    connect[0,:] = range(1,nEls+1) 
    connect[1,:] = connect[0,:] + 1
connect[0,7] = connect[0,8] = connect[0,9] = 8
connect[0,10] = 7
connect[0,11] = connect[0,12] = 6
connect[1,7] = 1
connect[1,8] = 2
connect[1,9] = connect[1,10] = connect[1,11] = 3
connect[1,12] = 4

gcor = np.zeros((2,nN))
gcor[1,0:5] = yMin
gcor[1,5:8] = yMax
for k in range(0,5):
    gcor[0,k] = xMin + k * (xMax-xMin)/4
    
gcor[0,7] = gcor[0,1]
gcor[0,6] = gcor[0,2]
gcor[0,5] = gcor[0,3]

alphas, vectors = FE.elementvectors(connect, gcor, nEls)
        
def objectives(X):
    """ compute the global cost, weight and compliance of the structure and round integer variables to the closest feasible """
    # X = [Young's modulus, Thicknesses]
    n = len(X) # 2*nEls
    if n%nEls != 0:
        return 'Error: Dimension must be divisible by nEls'
    
    Sections = [0] * nEls
    for k in range(0,nEls):
        Sections[k] = X[k+nEls]**2
    Costs = np.zeros(nEls)
    Weights = np.zeros(nEls)
    
    # Compute the volume of each elements
    for el in range(nEls):
        # volume of an element
        if vectors.any == None:
            n1 = int(connect[0,el]) # Node number 1
            n2 = int(connect[1,el]) # Node number 2
            V = Sections[el] * np.linalg.norm(gcor[:,n2-1] - gcor[:,n1-1])
        else: 
            V = Sections[el] * np.linalg.norm(vectors[:,el])
    
        #Rounding costs to nearest feasible ones
        X[el] = min(Youngs, key=lambda v:abs(v - X[el]))
        if X[el] == Emag:
            Weights[el] = V * rhomag
            Costs[el] = Cmag * Weights[el]
        elif X[el] == Etitan:
            Weights[el] = V * rhotitan
            Costs[el] = Ctitan * Weights[el]
        elif X[el] == Ealu:
            Weights[el] = V * rhoalu
            Costs[el] = Calu * Weights[el]
        elif X[el] == Esteel:
            Weights[el] = V * rhosteel
            Costs[el] = Csteel * Weights[el]
    U, F = FE.solveq(case,Fa,X[0:nEls],alphas,X[nEls:n],gcor,nEls,nN,connect)
            
    if case == 1:
        constr = abs(U[5]) - abs(Uy3max)
    elif case == 2:
        constr = abs(U[8]) - abs(Ux5max)
    if solver == 'nsga2' or solver == 'nsga3':
        if opt_pb == 'min cost':
            return [sum(Costs)], [constr]
        elif opt_pb == 'min weight':
            return [sum(Weights)], [constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs), sum(Weights)], [constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs), sum(Weights), np.dot(F, U)], [constr]
    elif solver == 'mads' or solver == 'miso' or solver == 'boa':
        if opt_pb == 'min cost':
            return [sum(Costs),constr]
        elif opt_pb == 'min weight':
            return [sum(Weights),constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs),sum(Weights),constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights),np.dot(F, U),constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs),np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs),sum(Weights),np.dot(F, U),constr]
    elif solver == 'cmaes':
        if opt_pb == 'min cost':
            return sum(Costs)
        elif opt_pb == 'min weight':
            return sum(Weights)
        elif opt_pb == 'min compliance':
            return np.dot(F, U)
        
        
def objectives_categ(X):
    """ Compute the global cost, weight and compliance of the structure with categorical variables
    If cycle == True, categorical variables can take any integer
    """
    # X = [Material numbers, Thicknesses]
    n = len(X) # 2*nEls
    if n%nEls != 0:
        return 'Error: Dimension must be divisible by nEls'
    # Material numbers (from 1 to 4)
    mat_nb = copy.copy(X[0:nEls])
#    if cycle:
#        for i in range(nEls):
#            mat_nb[i] = mat_nb[i]%4
#            if mat_nb[i] == 0:
#                mat_nb[i] = 4
        
    Sections = [0] * nEls
    for k in range(0,nEls):
        Sections[k] = X[k+nEls]**2
    Costs = np.zeros(nEls)
    Weights = np.zeros(nEls)
    Es = np.zeros(nEls)
    
    # Compute the volume of each elements
    for el in range(nEls):
        # volume of an element
        if vectors.any == None:
            n1 = int(connect[0,el]) # Node number 1
            n2 = int(connect[1,el]) # Node number 2
            V = Sections[el] * np.linalg.norm(gcor[:,n2-1] - gcor[:,n1-1])
        else: 
            V = Sections[el] * np.linalg.norm(vectors[:,el])
    
        #Rounding costs to nearest feasible ones
#        X[el] = min(Categ, key=lambda v: abs(v - X[el]))
        if mat_nb[el]%4 == Categ[0]%4: #titan
            Es[el] = Etitan
            Weights[el] = V * rhotitan
            Costs[el] = Ctitan * Weights[el]
        elif mat_nb[el]%4 == Categ[1]%4: #mag
            Es[el] = Emag
            Weights[el] = V * rhomag
            Costs[el] = Cmag * Weights[el]
        elif mat_nb[el]%4 == Categ[2]%4: #steel
            Es[el] = Esteel
            Weights[el] = V * rhosteel
            Costs[el] = Csteel * Weights[el]
        elif mat_nb[el]%4 == Categ[3]%4: #alu
            Es[el] = Ealu
            Weights[el] = V * rhoalu
            Costs[el] = Calu * Weights[el]

    U, F = FE.solveq(case,Fa,Es,alphas,X[nEls:n],gcor,nEls,nN,connect)
        
    if case == 1:
        constr = abs(U[5]) - abs(Uy3max)
    elif case == 2:
        constr = abs(U[8]) - abs(Ux5max)
    if solver == 'nsga2' or solver == 'nsga3':
        if opt_pb == 'min cost':
            return [sum(Costs)], [constr]
        elif opt_pb == 'min weight':
            return [sum(Weights)], [constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs), sum(Weights)], [constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs), sum(Weights), np.dot(F, U)], [constr]
    elif solver == 'mads' or solver == 'miso' or solver == 'boa':
        if opt_pb == 'min cost':
            return [sum(Costs),constr]
        elif opt_pb == 'min weight':
            return [sum(Weights),constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs),sum(Weights),constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights),np.dot(F, U),constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs),np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs),sum(Weights),np.dot(F, U),constr]
    elif solver == 'cmaes':
        if opt_pb == 'min cost':
            return sum(Costs)
        elif opt_pb == 'min weight':
            return sum(Weights)
        elif opt_pb == 'min compliance':
            return np.dot(F, U)
        
def objectives2(X):
    """ compute the global cost, weight and compliance of the structure and round integer variables to the closest feasible """
    """ To use for scaled Young's modulus """
    # X = [Young's modulus, Thicknesses]
    n = len(X) # 2*nEls
    if n%nEls != 0:
        return 'Error: Dimension must be divisible by nEls'
    
    Sections = [0] * nEls
    for k in range(0,nEls):
        Sections[k] = X[k+nEls]**2
    Costs = np.zeros(nEls)
    Weights = np.zeros(nEls)
    RealEs = np.zeros(nEls)
    
    # Compute the volume of each elements
    for el in range(nEls):
        # volume of an element
        if vectors.any == None:
            n1 = int(connect[0,el]) # Node number 1
            n2 = int(connect[1,el]) # Node number 2
            V = Sections[el] * np.linalg.norm(gcor[:,n2-1] - gcor[:,n1-1])
        else: 
            V = Sections[el] * np.linalg.norm(vectors[:,el])
    
        #Rounding costs to nearest feasible ones
        X[el] = min(Youngs2, key=lambda v:abs(v - X[el]))
        if X[el] == Emag2:
            Weights[el] = V * rhomag
            Costs[el] = Cmag * Weights[el]
        elif X[el] == Etitan2:
            Weights[el] = V * rhotitan
            Costs[el] = Ctitan * Weights[el]
        elif X[el] == Ealu2:
            Weights[el] = V * rhoalu
            Costs[el] = Calu * Weights[el]
        elif X[el] == Esteel2:
            Weights[el] = V * rhosteel
            Costs[el] = Csteel * Weights[el]
        for idx in range(nEls):
            RealEs[idx] = X[idx]*1e6
    U, F = FE.solveq(case,Fa,RealEs,alphas,X[nEls:n],gcor,nEls,nN,connect)
            
    if case == 1:
        constr = abs(U[5]) - abs(Uy3max)
    elif case == 2:
        constr = abs(U[8]) - abs(Ux5max)
    if solver == 'nsga2' or solver == 'nsga3':
        if opt_pb == 'min cost':
            return [sum(Costs)], [constr]
        elif opt_pb == 'min weight':
            return [sum(Weights)], [constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs), sum(Weights)], [constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs), np.dot(F, U)], [constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs), sum(Weights), np.dot(F, U)], [constr]
    elif solver == 'mads' or solver == 'miso':
        if opt_pb == 'min cost':
            return [sum(Costs),constr]
        elif opt_pb == 'min weight':
            return [sum(Weights),constr]
        elif opt_pb == 'min compliance':
            return [np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight':
            return [sum(Costs),sum(Weights),constr]
        elif opt_pb == 'min weight&compliance':
            return [sum(Weights),np.dot(F, U),constr]
        elif opt_pb == 'min cost&compliance':
            return [sum(Costs),np.dot(F, U),constr]
        elif opt_pb == 'min cost&weight&compliance':
            return [sum(Costs),sum(Weights),np.dot(F, U),constr]
    elif solver == 'cmaes':
        if opt_pb == 'min cost':
            return sum(Costs)
        elif opt_pb == 'min weight':
            return sum(Weights)
        elif opt_pb == 'min compliance':
            return np.dot(F, U)


#### Define objective functions according to variable type
if type_mat in ('Cycle_type', 'Categ_type'):
    cycle_bool = True
else:
    cycle_bool = False
" Problem definition "
if type_mat == 'Int_type':
    object_func = objectives
    Titan = Etitan
    Mag = Emag
    Steel = Esteel
    Alu = Ealu
elif type_mat == 'Int_type2':
    object_func = objectives2
    Titan = Etitan2
    Mag = Emag2
    Steel = Esteel2
    Alu = Ealu2
elif type_mat in ('Categ_type', 'Cycle_type'):
    object_func = objectives_categ
    Titan = 1
    Mag = 2
    Steel = 3
    Alu = 4
elif type_mat == 'Real_type':
    object_func = objectives
    Titan = Etitan
    Mag = Emag
    Steel = Esteel
    Alu = Ealu
