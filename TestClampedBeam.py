# -*- coding: utf-8 -*-

######### Simple example to test ClampedBeam.py

import os
import time

folder = os.getcwd() + "/"

case = 1 # 1 or 2
func = 'weight' # 'cost', 'weight', 'compliance', 'cost&weight', 'weight&compliance',
            # 'cost&compliance' or 'cost&weight&compliance'
exec(open(folder+'Truss/ClampedBeam.py').read())

" Plot nodes "
FE.displaynodes(gcor)
  
" Plot mesh "
FE.displaymesh(connect, gcor)

Thicknesses = 0.027 * np.ones(nEls)
" Initial design "
Es = np.zeros(nEls) # Vector of Young modules for each element
Es[0:4] = Etitan
Es[4] = Es[10:13] = Ealu
Es[7:10] = Emag
Es[5:7] = Esteel

X = np.concatenate((Es,Thicknesses))
" Display materials on mesh "
FE.displaymaterials(X, Etitan, Emag, Esteel, Ealu, connect, gcor)
U, F = FE.solveq(case,Fa,Es,alphas,Thicknesses,gcor,nEls,nN,connect)

" New coordinates "
gcor = FE.changegcor(gcor, U, nN)

" Print the modification of the mesh "
FE.displaymesh(connect, gcor)
FE.displaymaterials(X, Etitan, Emag, Esteel, Ealu, connect, gcor)

