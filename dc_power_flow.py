#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 17:47:18 2020

"""

from gurobipy import *
import numpy as np
import scipy.io as sio


#%% Get data from .mat file

baseMVA = 100 # MVA base for per-unit conversion
caseData = sio.loadmat('case145.mat')
pgl_max = np.array(caseData['pglmax'])
pg_min = np.array(caseData['pgmin'])
pl = np.array(caseData['PLT'])
cg = np.array(caseData['cg'])
phi_max = np.array(caseData['phi_max'])
DF = np.array(caseData['DF'])
V = np.array(caseData['V'])
SU = np.array(caseData['SU'])
SD = np.array(caseData['SD'])
DFs = np.array(caseData['DFs'])

#%% Initialize data arrays

ng = np.size(pg_min,0) # Number of generators
nseg = np.size(pgl_max,1) # Number of segments in linear approximation cost fcn
nsc = np.size(DFs,2) # Number of scenarios
nt = np.size(pl,1) # Number of periods in planning horizon
nb = np.size(pl,0) # Number of busses in the system

G = range(ng)
L = range(nseg)
S = range(nsc)
T = range(nt)
N = range(nb)

#%% Gurobi model

m = Model('DCOPF')

#============
# Decision Vars
#============
u = m.addVars(G,T, vtype=GRB.BINARY,name='u',
              lb=0,ub=1.0) # u_{g}
pg = m.addVars(G,T, vtype=GRB.CONTINUOUS,name='pg',
              lb=-GRB.INFINITY,ub=GRB.INFINITY) # p_{g,t}
pgl = m.addVars(G,T,L, vtype=GRB.CONTINUOUS,name='pgl',
              lb=0,ub=GRB.INFINITY) # p_{g,t,l}
sug = m.addVars(G,T, vtype=GRB.CONTINUOUS,name='sug',
              lb=0,ub=GRB.INFINITY) # SU_{g,t}
sdg = m.addVars(G,T, vtype=GRB.CONTINUOUS,name='sdg',
              lb=0,ub=GRB.INFINITY) # SD_{g,t}


#============
# Constraints
#============
# Coupling constraints
m.addConstrs((pgl[g,t,l]<=u[g,t]*pgl_max[g,l] for g in G for t in T for l in L))

# ED constraints
m.addConstrs((quicksum(pgl[g,t,l] for l in L)-baseMVA*pg[g,t]==-pg_min[g]
                                 for t in T for g in G),name='pcwLin') #pcwLin
m.addConstrs((quicksum(pg[g,t] for g in G)==quicksum(pl[i,t] 
                           for i in N) for t in T),'power_bal') #power balance
for s in S:
    for t in T:
        DF = DFs[:,:,s]
        bV = np.reshape(DF@pl[:,t],[np.size(phi_max,0),1])
        m.addMConstrs(DF@V,pg.select('*',t),'<=',phi_max + bV,'flowLim+') # Transmission flow
        m.addMConstrs(-DF@V,pg.select('*',t),'<=',phi_max - bV,'flowLim-') # Transmission flow

# UC constraints
m.addConstrs((sug[g,0]>=u[g,0]*SU[g] for g in G),'SU0') # startup at t=0
m.addConstrs((sug[g,t]>=(u[g,t]-u[g,t-1])*SU[g] for g in G for t in T[1:]),'SUT') # startup at t>0
m.addConstrs((sdg[g,t]>=(u[g,t-1]-u[g,t])*SD[g] for g in G for t in T[1:]),'SDT') # shutdown


#============
# Objective
#============
m.modelSense = GRB.MINIMIZE
m.setObjective(quicksum(quicksum(sdg[g,t] + sug[g,t]+ u[g,t]*cg[g,0] +
                                 quicksum(pgl[g,t,l]*cg[g,1+l] for l in L) for g in G)for t in T))


# Run optimization
m.update()
m.optimize()

#%% Print solution
print('Operating cost: %g' % m.getObjective().getValue())
#print('Unit commitment solution')
#for v in u.select():
#    print('%s %g' % (v.varName,v.x))
#print('Economic dispatch solution')
#for v in pg.select():
#    print('%s %g' % (v.varName,baseMVA*v.x))
