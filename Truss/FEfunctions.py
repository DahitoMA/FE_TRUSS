#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import scipy.sparse.linalg as sl
import numpy as np
import copy
import matplotlib.pyplot as plt

def changegcor(gcor, U, nN=None):
    if nN == None:
        nN = int(len(U)/2)
    for j in range(nN):
        gcor[0,j] += U[2*j]
        gcor[1,j] += U[2*j+1]
    return gcor
        
def elementvectors(connect, gcor, nEls=None):
    if nEls == None:
        nEls = int(np.shape(connect)[1])
    vectors = np.zeros((2,nEls)) # Corresponding vector for each element
    alphas = np.zeros(nEls) # Vector of angles for each element
    for ne in range(nEls):
        n1 = int(connect[0,ne]) # Node number 1
        n2 = int(connect[1,ne]) # Node number 2
        vectors[:,ne] = gcor[:,n2-1] - gcor[:,n1-1]
        v = 1/np.linalg.norm(vectors[:,ne]) * vectors[:,ne] # Normalized vector
        alphas[ne] = np.sign(v[1]) * np.arccos(np.dot(v, [1,0]))
    return alphas, vectors

def stiffness(Es,alphas,Thicknesses,gcor,nEls,nN,elDof,connect,dFreedom):
    # Build stiffness matrix K
    K = np.zeros((2*nN,2*nN))
    
    for ne in range(0,nEls):
        Sections = [0] * nEls
        for k in range(0,nEls):
            Sections[k] = Thicknesses[k]**2
        S = Sections[ne]
        n1 = int(connect[0,ne]) # Node number 1
        n2 = int(connect[1,ne]) # Node number 2
        le = np.linalg.norm(gcor[:,n2-1] - gcor[:,n1-1])
        E = Es[ne] # Young's modulus
        R = np.zeros((2,2))
        alpha = alphas[ne]
        
        R[0,0] = np.cos(alpha)
        R[1,1] = R[0,0]
        R[0,1] = np.sin(alpha)
        R[1,0] = - R[0,1]
                            
        Kiiloc = np.zeros((2,2))
        Kijloc = np.zeros((2,2))
        Kjjloc = np.zeros((2,2))
        Kiiloc[0,0] = (E*S)/le
        Kiiglo = np.matmul(np.transpose(R) , np.matmul(Kiiloc,  R))
        Kijloc = - Kiiloc
        Kijglo = np.matmul(np.transpose(R) , np.matmul(Kijloc,  R)) 
        Kjjloc = Kiiloc
        Kjjglo = np.matmul(np.transpose(R) , np.matmul(Kjjloc,  R))
                    
        idx1 = 2*(n1-1)
        idx2 = 2*(n2-1)
        K[idx1:idx1+2,idx1:idx1+2] += Kiiglo
        K[idx2:idx2+2,idx2:idx2+2] += Kjjglo
        K[idx1:idx1+2,idx2:idx2+2] += Kijglo
        K[idx2:idx2+2,idx1:idx1+2] = K[idx1:idx1+2,idx2:idx2+2]
        
    return K

def solveq(case,Fa,Es,alphas,Thicknesses,gcor,nEls,nN,connect):
    if case != 1 and case != 2:
        print('error: case must be int 1 or 2')
        return
                
    " Set boundary conditions "    
    if case == 1:
        nB = 4 # Number of boundary elements
        bEls = np.zeros((nB,1)) # Boundary elements
        bEls[0] = 1
        bEls[1] = 8
        bEls[2] = 4
        bEls[3] = 5
        bNod = np.zeros((2,1))  # Boundary nodes
        bNod[0] = 1
        bNod[1] = 5
        
    elif case == 2:
        nB = 2 # Number of boundary elements
        bEls = np.zeros((nB,1)) # Boundary elements
        bEls[0] = 1
        bEls[1] = 8
        bNod = [1]
    
    " Set degrees of freedom for each element "
    pDeg = np.ones((nEls, 1)) # Polynomial degree
    elDof = 2 * (pDeg + 1)
    dFreedom = connect

    " Build K "
    Fi = np.zeros(2*nN)
    K = stiffness(Es,alphas,Thicknesses,gcor,nEls,nN,elDof,connect,dFreedom)

    " Eliminations in K and building of F"    
    if case == 1:
        # Delete columns and lines relative to nodes 1 and 5
        K = np.delete(K, 8, 0)
        K = np.delete(K, 8, 1)
        K = np.delete(K, 8, 0)
        K = np.delete(K, 8, 1)
        K = np.delete(K, 0, 0)
        K = np.delete(K, 0, 1)
        K = np.delete(K, 0, 0)
        K = np.delete(K, 0, 1)
        # Build F
        Fi[5] = Fa
        F = Fi
        F = np.delete(F, 8, 0)
        F = np.delete(F, 8, 0)
        F = np.delete(F, 0, 0)
        F = np.delete(F, 0, 0)
        
    elif case ==2:
        # Delete columns and lines relative to node 1
        K = np.delete(K, 0, 0)
        K = np.delete(K, 0, 1)
        K = np.delete(K, 0, 0)
        K = np.delete(K, 0, 1)
        # Build F
        Fi[8] = -Fa
        F = Fi
        F = np.delete(F, 0, 0)
        F = np.delete(F, 0, 0)
        
    " Solution without boundary conditions "
    U = sl.cg(K,F)[0]

    " Insertion of boundary conditions "
    F = Fi
    
    if case == 1:
        Ui = np.zeros(2*nN)
        Ui[0] = Ui [1] = 0
        Ui[8] = Ui[9] = 0
        Ui[2:8] = U[0:6]
        Ui[10:len(Ui)] = U[6:len(U)]
        U = Ui
    elif case ==2:
        Ui = np.zeros(2*nN)
        Ui[0] = Ui[1] = 0
        Ui[2:len(Ui)] = U
        U = Ui
        
    return U, F

def dispviolation(case, U, Uy3max, Ux5max):
    if case == 1:
        if abs(U[5]) > abs(Uy3max):
            return abs(U[5]) - Uy3max
        else:
            return 0
    elif case == 2:
        if abs(U[8]) > abs(Ux5max):
            return abs(U[8]) - abs(Ux5max)
        else:
            return 0    
    
def displaymesh(connect, gcor):
    plt.figure()
    for k in range(0, np.shape(connect)[1]):
        a = int(connect[0,k] - 1)
        b = int(connect[1,k] - 1)
        plt.plot([gcor[0,a], gcor[0,b]], [gcor[1,a], gcor[1,b]])  
    

def displaynodes(gcor):
    plt.figure()
    for k in range(0,np.shape(gcor)[1]):
        plt.plot(gcor[0, k], gcor[1, k], 'o')
        
def displaymaterials(sol, Titan, Mag, Steel, Alu, connect, gcor, cycle=False):
    coltitan = 'green'
    colmag = 'red'
    colsteel = 'blue'
    colalu = 'purple'
    nEls = np.shape(connect)[1]
    mat = copy.copy(sol[0:nEls])
    thick = sol[nEls:2*nEls]
    colormat = []
    thickcolor = []
    if cycle:
        for k in range(nEls):
            mat[k] = mat[k]%4
            if mat[k] == 0:
                mat[k] = 4
    for l in range(nEls):
        if thick[l] < 0.02:
            thickcolor.append(0.5)
        elif thick[l] < 0.03:
            thickcolor.append(1.5) 
        elif thick[l] < 0.04:
            thickcolor.append(2.5)
        else:
            thickcolor.append(3.5)
    for j in range(nEls):
        if mat[j] == Titan:
            colormat.append(coltitan)
        elif mat[j] == Mag:
            colormat.append(colmag)
        elif mat[j] == Steel:
            colormat.append(colsteel)
        elif mat[j] == Alu:
            colormat.append(colalu)                  
    plt.figure()
    for k in range(nEls):
        a = int(connect[0,k] - 1)
        b = int(connect[1,k] - 1)
        plt.plot([gcor[0,a], gcor[0,b]], [gcor[1,a], gcor[1,b]], colormat[k], linewidth=thickcolor[k])
    plt.title('Colored mesh according to material type and thickness range') 
    legend_elements = [plt.Line2D([0], [0], color='g', label='titan'),
                   plt.Line2D([0], [0], color='r', label='mag'),
                   plt.Line2D([0], [0], color='b', label='steel'),
                   plt.Line2D([0], [0], color='purple', label='alu')]
    plt.legend(handles=legend_elements, loc='upper right', prop={'size': 12})
    
def strtofloatmatrix(filename, nb_runs):
    " Convert a data file into a float matrix "
    with open(filename, 'r+') as f:
        lines = f.readlines()
        length = len(lines) # convert file into matrix
        mtx = np.zeros((length, nb_runs))
        mtxmean = np.zeros(length) # mean of all runs
        mtxmed = np.zeros(length) # median of all runs
        mtx25 = np.zeros(length) # fisrt quartile
        mtx75 = np.zeros(length) # third quartile
        for k in range(length):
            for j in range(nb_runs):
                mtx[k,j] = float(lines[k].split()[j])
        for i in range(length):
            vect = [mtx[i,j] for j in (np.nonzero(mtx[i,:])[0])] # nonzero components
            mtxmean[i] = np.mean(vect)
            mtxmed[i] = np.median(vect)
            mtx25[i] = np.percentile(vect,25)
            mtx75[i] = np.percentile(vect,75)
    return mtx, mtxmean, mtxmed, mtx25, mtx75

