from gurobipy import *
import numpy as np
import scipy.io as sio

#dynamic class
class Expando(object):
    pass


# ### Master Problem

# In[22]:


class Master:
    '''
    Parameters are 
    Pmax_gl : array GxL       pgl_max = np.array(caseData['pglmax'])
    Pmin_g : array G          pg_min = np.array(caseData['pgmin'])
    Cmin_g : array G
    C_gl : array GxL          cg = np.array(caseData['cg'])
    SU_g : array G            SU = np.array(caseData['SU'])
    SD_g : array G            SD = np.array(caseData['SD'])
    PL_it : array NxT         pl = np.array(caseData['PLT'])
    DF_s_bi : array SxBxN (eg S arrays of size BxN)   DF = np.array(caseData['DF'])
    V_ig : array NxG          V = np.array(caseData['V'])
    Phimax_b : array B        phi_max = np.array(caseData['phi_max'])
    along with the sets G,N,B,T,L,S 
    '''
    def __init__(self,Pmax_gl,Pmin_g,Cmin_g,C_gl,SU_g,SD_g,PL_it,DF_s_bi,V_ig,Phimax_b):#,M,N,P,b,A,c,G,h):
        #sets
        self.G=range(Pmin_g.shape[0])
        self.N=range(PL_it.shape[0])
        self.B=range(Phimax_b.shape[0])
        self.T=range(PL_it.shape[1])
        self.L=range(C_gl.shape[1])
        self.S=range(DF_s_bi.shape[2])
        
        #variables and constraints
        self.variables = Expando()
        self.constraints = Expando()
        #Data
        self.data = Expando()
        self._load_data()
        self._init_algo_data()
        #Model
        self._build_model()
        
    #Load the Data and Parameters
    def _load_data(self):
        self._load_coupling_params()
        self._load_continuous_params()
        
        self.data.SU_g=SU_g
        self.data.SD_g=SD_g
        
    def _load_coupling_params(self):
        self.data.Pmax_gl=Pmax_gl
        self.data.Cmin_g=Cmin_g
    def _load_continuous_params(self):
        self.data.Pmin_g=Pmin_g
        self.data.PL_it=PL_it
        self.data.DF_s_bi=DF_s_bi
        self.data.V_ig=V_ig
        self.data.Phimax_b=Phimax_b
        self.data.C_gl=C_gl
        

    #Initialize Algorithm successive variables
    def _init_algo_data(self):
        #stores the successive cuts added to the master pb
        self.data.cutlist = []
        
        #current bounds for the optimal solution to the original problem
        self.data.ub = GRB.INFINITY #upperbound of the optimal solution to the original pb
        self.data.lb = -GRB.INFINITY #lowerbound of the optimal solution to the original pb
        
        #stores the sequence of bounds for the optimal solution to the original problem
        self.data.ubs = []
        self.data.lbs = []
        #stores the sequence of variables x,y,and SPu of the successive (RMP) and (SP)
        self.data.xs = []
        self.data.ys = []
        self.data.SPus = []
        

    """
    Build the model
    """
    #Variables
    def _set_variables(self):
        m = self.model
        self.variables.u_gt = m.addVars(self.G,self.T,vtype=GRB.BINARY, name='u_gt')
        self.variables.SD_gt = m.addVars(self.G,self.T,vtype=GRB.CONTINUOUS, name='SD_gt')
        self.variables.SU_gt = m.addVars(self.G,self.T,vtype=GRB.CONTINUOUS, name='SU_gt')
        self.variables.mu = m.addVars(self.T,name='mu') #this will correspond to cost vector ^T the (Pgtl)s for each t in T
        m.update()
        
        
    #Objective function
    def _set_objective(self):
        self.model.setObjective(
            quicksum(self.data.Cmin_g[g]*self.variables.u_gt[g,t]+self.variables.SU_gt[g,t]+self.variables.SD_gt[g,t] for g in self.G for t in self.T)
            +quicksum(self.variables.mu[t] for t in self.T) ,
            GRB.MINIMIZE)
        

    #Constraints
    def _set_constraints(self):

        u_gt=self.variables.u_gt 
        SD_gt=self.variables.SD_gt
        SU_gt=self.variables.SU_gt
        SU_g = self.data.SU_g
        SD_g = self.data.SD_g

        # Pure UC constraints (independent of sub-problems)
        self.constraints.uc ={}
        self.constraints.uc['SU0'] = self.model.addConstrs((SU_gt[g,0]>=u_gt[g,0]*SU_g[g] for g in self.G),'SU0') # startup at t=0
        self.constraints.uc['SUT'] = self.model.addConstrs((SU_gt[g,t]>=(u_gt[g,t]-u_gt[g,t-1])*SU_g[g] for g in self.G for t in self.T[1:]),'SUT') # startup at t>0
        self.constraints.uc['SDT'] = self.model.addConstrs((SD_gt[g,t]>=(u_gt[g,t-1]-u_gt[g,t])*SD_g[g] for g in self.G for t in self.T[1:]),'SDT') # shutdown

        # Cuts that will be added by sub-problems
        self.constraints.cuts = {} #empty set of constraints for the initial RPM
    
    #Model
    def _build_model(self):
        self.model = Model()
        self._set_variables()
        self._set_objective()
        self._set_constraints()
        self.model.setParam(GRB.Param.OutputFlag,0)
        self.model.update()
    
            
    #Updates bounds on the optimal solution to the original problem
    def _update_bounds(self):
        #z_sub = self.submodel.model.ObjVal
        z_subs = [self.submodel[t].model.ObjVal for t in self.T]
        z_master = self.model.ObjVal
        #self.data.ub = z_master - self.variables.mu.x + z_sub
        self.data.ub = z_master + quicksum(- self.variables.mu[t].x + z_sub[t] for t in self.T)
        self.data.lb = self.model.ObjBound

        self.data.ubs.append(self.data.ub)
        self.data.lbs.append(self.data.lb)

    
    def optimize(self, simple_results=False):
        
        m = self.model
        
        # Initial solution
        cont = 1
            
        dontStop = True
        while dontStop:
            dontStop = False
            
            #==========
            # Solve master problem
            #==========
            m.update() # Update since we added cuts
            m.optimize() #from this we get an optimal solution mu bar, x bar
            if (m.status==GRB.INFEASIBLE) or (m.status==GRB.INF_OR_UNBD) :
                print('Problem infeasible!!!!!')
                return
            
            
            if not hasattr(self,'submodel'):
                # Initialize submodels
                self.submodel = {}
                for t in self.T:
                    self.submodel[t] = Subproblem(self,t) # Build an instance of the subproblem (SP) from the initial solution

            else:
                # Update values of variables u_gt in subproblem
                for t in self.T:
                    for g in self.G:
                        for l in self.L:
                            self.submodel[t].constraints.coupl[g,l].rhs = self.variables.u_gt[g,t].x*self.data.Pmax_gl[g,l]
                    self.submodel[t].model.update()
                m.update()
            m.optimize()
                    
            for t in self.T:
                subp = self.submodel[t] # For brevitiy below
                # Solve subproblem
                subp.optimize()
                
                # Case 1: subproblem infeasible
                if subp.model.status==GRB.INFEASIBLE or subp.model.status==GRB.INF_OR_UNBD :
                    dontStop = True # We'll have to iterate once more
                    rBar = np.array(subp.model.FarkasDual) # Unbounded ray of the dual
                    rCoupling = np.reshape(rBar[0:len(self.G)*len(self.L)],[len(self.G),len(self.L)]) # Get components corresponding to coupling constraints
                    constrs = subp.model.getConstrs()
                    # Add cut r(b-Ax)<=0, 
                    # We are skipping the first |G|*|L| values of rBar, since these are the coupling constraints
                    # , which have rhs b equal to zero but constrs.rhs nonzero
                    # because of the Ax part of (b-Ax)
                    m.addConstr(-quicksum(rBar[mi]*constrs[mi].rhs for mi in range(len(self.G)*len(self.L)+1,len(constrs)))
                                -quicksum(self.variables.u_gt[g,t]*self.data.Pmax_gl[g,l]*rCoupling[g,l] for g in self.G for l in self.L)<=0)
#                    print('%d\t%d\t%d' % (cont,t,subp.model.status))
                elif subp.model.status==GRB.OPTIMAL:
                # Case 2: subproblem solved to optimality
                    if self.variables.mu[t].x-subp.model.ObjVal<-1E-3: # Using arbitrary tolerance of 1E-3
                        # Optimality constraint violated, add cut
                        dontStop = True # We'll have to iterate once more
                        uBar = np.array(subp.model.Pi) # Shadow prices (optimal solution of the dual)
                        uCoupling = np.reshape(uBar[0:len(self.G)*len(self.L)],[len(self.G),len(self.L)]) # Get components corresponding to coupling constraints
                        constrs = subp.model.getConstrs()
                        # Add cut u(b-Ax)<= mu
                        # We are skipping the first |G|*|L| values of rBar, since these are the coupling constraints
                        # , which have rhs b equal to zero but constrs.rhs nonzero
                        # because of the Ax part of (b-Ax)
                        m.addConstr(quicksum(uBar[mi]*constrs[mi].rhs for mi in range(len(self.G)*len(self.L)+1,len(constrs)))
                                +quicksum(self.variables.u_gt[g,t]*self.data.Pmax_gl[g,l]*uCoupling[g,l] for g in self.G for l in self.L)<=
                                self.variables.mu[t])
#                    print('%d\t%d\t%d\t%4.2f' % (cont,t,subp.model.status,self.variables.mu[t].x))
                cont = cont + 1
                    



# ### Subproblem

# In[6]:


# Subproblem
class Subproblem:
    def __init__(self, RMP,t):
        #sets
        self.G=RMP.G
        self.N=RMP.N
        self.B=RMP.B
        self.T=RMP.T
        self.L=RMP.L
        self.S=RMP.S
        
        # List of contingency scenarios that will be added to the set of 
        # constraints. Initially, only the base case is considered (s=0)
        self.contScenarios = []
        self.contScenarios.append(0)
        
        # Time period
        self.t = t
        
        # Base MVA
        self.baseMVA = 100
        
        #variables and constraints
        self.variables = Expando()
        self.constraints = Expando()
        #Data
        self.data = Expando()
        #RMP
        self.RMP = RMP
        #Model
        self._build_model()

    def optimize(self):
        m = self.model
        P_gtl = self.variables.P_gtl
        Pmin_g = self.RMP.data.Pmin_g
        PL_it = self.RMP.data.PL_it
        DF_s_bi = self.RMP.data.DF_s_bi
        V_ig = self.RMP.data.V_ig
        Phimax_b = self.RMP.data.Phimax_b
        P_gt = self.variables.P_gt
        t = self.t
        
        start = True
        overload = False
        while (start or overload):
            start = False
            if overload:
                m.update()
                overload = False
        
            m.optimize()
            if (m.status==GRB.INFEASIBLE) or (m.status==GRB.INF_OR_UNBD):
                # If infeasible, get back to master problem
                # Otherwise, go on
                return

            for s in self.S[1:]:
                # For each contingency, check overloads of branches
                if s in self.contScenarios:
                    # No need to check scenarios that are already in the subproblem
                    # formulation, those are guaranteed to be fine
                    continue
                
                DF = DF_s_bi[:,:,s]
                bV = np.reshape(DF@PL_it[:,t],[np.size(Phimax_b,0),1])
                P_gt_val = np.array([P_gt[g].x for g in self.G])# value of P_gt from solution
                flowVector = DF@V_ig@P_gt_val-bV
                if np.any(abs(flowVector)>Phimax_b):
                    # If there is overload, add scenario to problem and re-solve
                        overload = True
                        self.contScenarios.append(s)
                        m.addMConstrs(DF@V_ig,P_gt.select('*'),'<=',Phimax_b + bV,'flowLim+') # Transmission flow
                        m.addMConstrs(-DF@V_ig,P_gt.select('*'),'<=',Phimax_b - bV,'flowLim-') # Transmission flow
                    
        
    """
    Build Subproblem (SP)
    """
    def _set_variables(self):
        m = self.model

        self.variables.P_gt = m.addVars(self.G, vtype=GRB.CONTINUOUS,name='pg',
              lb=-GRB.INFINITY,ub=GRB.INFINITY) # p_{g,t}
        self.variables.P_gtl = m.addVars(self.G,self.L, vtype=GRB.CONTINUOUS,name='pgl',
              lb=0,ub=GRB.INFINITY) # p_{g,t,l}
        m.update()

    def _set_objective(self):
        m = self.model

        C_gl = self.RMP.data.C_gl
        P_gtl = self.variables.P_gtl
        
        m.setObjective(quicksum(C_gl[g,l]*P_gtl[g,l] for g in self.G for l in self.L))

    def _set_constraints(self):

        m = self.model
        P_gtl = self.variables.P_gtl
        Pmin_g = self.RMP.data.Pmin_g
        PL_it = self.RMP.data.PL_it
        DF_s_bi = self.RMP.data.DF_s_bi
        V_ig = self.RMP.data.V_ig
        Phimax_b = self.RMP.data.Phimax_b
        P_gt = self.variables.P_gt
        Pmax_gl = self.RMP.data.Pmax_gl
        t = self.t
        
        # Coupling constraints
        # ====================
        # Important: the uc_var constraints must be added first
        # We are counting on this to retrieve the unbounded rays when adding
        # cuts
        # ====================
        self.constraints.coupl = m.addConstrs((P_gtl[g,l]<=self.RMP.variables.u_gt[g,t].x*Pmax_gl[g,l] for g in self.G for l in self.L),'uc_var')

        # ED constraints
        # ==============
        # Pg = sum of piecewise linear bits
        self.constraints.pcwLin = m.addConstrs((quicksum(P_gtl[g,l] for l in self.L)-self.baseMVA*P_gt[g]==-Pmin_g[g] for g in self.G),name='pcwLin') #pcwLin
        # Power balance
        self.constraints.powBal = m.addConstr(quicksum(P_gt[g] for g in self.G)==quicksum(PL_it[i,t] 
                                   for i in self.N),'power_bal') #power balance
        
        # Flow constraints in all scenarios considered
        for s in self.contScenarios:
            DF = DF_s_bi[:,:,s]
            bV = np.reshape(DF@PL_it[:,t],[np.size(Phimax_b,0),1])
            m.addMConstrs(DF@V_ig,P_gt.select('*'),'<=',Phimax_b + bV,'flowLim+') # Transmission flow
            m.addMConstrs(-DF@V_ig,P_gt.select('*'),'<=',Phimax_b - bV,'flowLim-') # Transmission flow
            

    def _build_model(self):
        self.model = Model()
        self._set_variables()
        self._set_objective()
        self._set_constraints()
        self.model.setParam(GRB.Param.OutputFlag,0)
        self.model.setParam(GRB.Param.InfUnbdInfo,1)
        self.model.update()

    def update_fixed_vars(self, RMP=None):
        pass


# ### To run the algorithm




#m = Master(M,N,P,b,A,c,G,h)
 
    
baseMVA = 100 # MVA base for per-unit conversion
caseData = sio.loadmat('case145.mat')
Pmax_gl = np.array(caseData['pglmax'])
Pmin_g = np.array(caseData['pgmin'])
PL_it = np.array(caseData['PLT'])
cg = np.array(caseData['cg'])
Phimax_b = np.array(caseData['phi_max'])
#DF = np.array(caseData['DF'])
V_ig = np.array(caseData['V'])
SU_g = np.array(caseData['SU'])
SD_g = np.array(caseData['SD'])
DF_s_bi = np.array(caseData['DFs'])    
C_gl = cg[:,1:]
Cmin_g = cg[:,0]

m=Master(Pmax_gl,Pmin_g,Cmin_g,C_gl,SU_g,SD_g,PL_it,DF_s_bi,V_ig,Phimax_b)
m.optimize()
#%%
if m.model.status ==GRB.OPTIMAL:
    totCost = m.model.ObjVal
#    for t in m.T:
#        totCost += m.submodel[t].model.ObjVal
    print('Cost Benders Solution %4.3f' % totCost)





