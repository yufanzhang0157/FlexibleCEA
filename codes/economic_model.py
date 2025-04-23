#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyomo.environ import *

import numpy as np

def PFALopt(Vout,Cout,Hout,XT0,XH0,XC0,XN0,XS0,pe=0.2051,pc=0.2,T=24,M=100000):
    
   
    #Vout: outdoor temperature
    #Cout: outdoor CO2
    #Hout: outdoor humidity
    
    ## photosynthesis parameter
    cepsilon = 17*10**(-9)
    cGamma = 7.32*10**(-5)
    cQ10Gamma = 2
    cpar = 0.5
    cradrf = 0.4
    cbnd = 0.004
    cstm = 0.007
    ccar1 = -1.32*10**(-5)
    ccar2 = 5.94*10**(-4)
    ccar3 = -2.64*10**(-3)
    ck = 0.9
    clars = 75
    ctau = 0.07
    
    # respiration parameter
    cresps = 3.47*10**(-7)
    crespr = 1.16*10**(-7)
    cQ10resp = 2
    
    # growth rate coefficient
    crgrmax = 5*10**(-6)
    cQ10gr = 1.6
    
    # uptaken carbon coefficient
    calpha = 0.68
    cbeta = 0.8
    
    ## transpiration parameter
    capl = 62.8
    cvplal = 3.6*10**(-3)
    cv0 = 0.85
    cv1 = 611
    cH20 = 18
    cr = 8314
    ctabs = 273
    cv2 = 17.4
    cv3 = 239
    
    ## heat parameter
    crad = 0.2
    
    
    c_lat_water = 2256.4 # [kJ/kg]
    c_U = 0.3 # [W/m^2 C]
    c_dehum_eev = 3
    fan_eff = 1000.*0.00047194745*15
    
    ## state dynamic parameter
    ccapq = 30000

    c_Length = 12.2 # [m]
    c_Width = 2.5 # [m]
    c_Height = 3.0 # [m]
    c_surface_area = (c_Length*c_Width +2*(c_Width*c_Height + c_Length*c_Height))
    c_volume = c_Length*c_Width*c_Height
    c_grow_area = 80.0 # m^2
    
    model = ConcreteModel()
    
    
    model.Uqs = Var(range(T),bounds=(0,212*2)) # heat supply
    model.Uqr = Var(range(T),bounds=(0,212*2)) # heat removal
    #model.Uhs = Var(range(T),bounds=(0,0)) # humidity supply
    model.Uhr = Var(range(T),bounds=(0,2.1*10**(-5)*2)) # humidity removal
    model.Uc = Var(range(T),bounds=(0,4.6*10**(-6)*2)) # CO2 supply
    model.Ul = Var(range(T),bounds=(0,100*2)) # electricity input for lighting
    model.Uv = Var(range(T),bounds=(0,10**(-2)*2)) # ventilation rate
    

    
    model.Xn = Var(range(T),bounds=(10**(-8),None))
    model.Xs = Var(range(T),bounds=(10**(-8),None))
    model.Xt = Var(range(T),bounds=(18,25))
    model.Xc = Var(range(T),bounds=(800*0.04401/24.45/10**3,1200*0.04401/24.45/10**3))
    model.Xh = Var(range(T),bounds=(0,None))
    model.Rh = Var(range(T),bounds=(70,80))
    
    model.Epsilont = Var(range(1,T),bounds=(0,0))
    model.Epsilonc = Var(range(1,T),bounds=(0,0))
    model.EpsilonRh = Var(range(1,T),bounds=(0,0))
    

    def TempepsDark_UB_constraint_rule(model, t):
        return model.Xt[t] <= 20
    model.TempepsDark_UB_constraint = Constraint([1,2,3,20,21,22,23],                                              rule=TempepsDark_UB_constraint_rule)
    
    def TempepsDark_LB_constraint_rule(model, t):
        return 18 <= model.Xt[t] 
    model.TempepsDark_LB_constraint = Constraint([1,2,3,20,21,22,23],                                              rule=TempepsDark_LB_constraint_rule)
    
    def TempepsDay_UB_constraint_rule(model, t):
        return model.Xt[t] <= 25
    model.TempepsDay_UB_constraint = Constraint(range(4,20),                                              rule=TempepsDay_UB_constraint_rule)
    
    def TempepsDay_LB_constraint_rule(model, t):
        return 22 <= model.Xt[t] 
    model.TempepsDay_LB_constraint = Constraint(range(4,20),                                              rule=TempepsDay_LB_constraint_rule)
    
    def Lightdark_constraint_rule(model, t):
        return model.Ul[t] == 0
    model.Lightdark_constraint = Constraint([0,1,2,3,20,21,22,23],                                              rule=Lightdark_constraint_rule)  
    
    def Carboneps_UB_constraint_rule(model, t):
        return model.Xc[t] <= 1200*0.04401/24.45/10**3
    model.Carboneps_UB_constraint = Constraint(range(1,T), rule=Carboneps_UB_constraint_rule)
    
    def Carboneps_LB_constraint_rule(model, t):
        return 800*0.04401/24.45/10**3 <= model.Xc[t] 
    model.Carboneps_LB_constraint = Constraint(range(1,T), rule=Carboneps_LB_constraint_rule)
    
    def Rheps_UB_constraint_rule(model, t):
        return model.Rh[t] <= 80
    model.Rheps_UB_constraint = Constraint(range(1,T), rule=Rheps_UB_constraint_rule)
    
    def Rheps_LB_constraint_rule(model, t):
        return 70 <= model.Rh[t] 
    model.Rheps_LB_constraint = Constraint(range(1,T), rule=Rheps_LB_constraint_rule)
    

    
    
    model.Gamma = Var(range(T),bounds=(10**(-18),None))
    model.epsilon = Var(range(T),bounds=(10**(-18),None))
    model.sigmaco2 = Var(range(T),bounds=(10**(-18),None))
    model.phi_photomax = Var(range(T),bounds=(10**(-18),None))
    model.phi_photo = Var(range(T),bounds=(10**(-18),None))
    model.phi_resp = Var(range(T),bounds=(10**(-18),None))
    model.rgr = Var(range(T),bounds=(10**(-18),None))
    model.phi_calpl = Var(range(T),bounds=(10**(-18),None))
    
    model.phi_hplal = Var(range(T),bounds=(10**(-18),None))
    
    #model.EpsilonRh = Var(range(1,T),bounds=(10**(-18),5))
    
    ### photosynthesis
    def Gamma_constraint_rule(model, t):
        return cGamma*cQ10Gamma**((model.Xt[t]-20)/10) == model.Gamma[t] 
    
    model.Gamma_constraint = Constraint(range(T), rule=Gamma_constraint_rule)
    
    def epsilon_constraint_rule(model, t):
        return cepsilon*(model.Xc[t]-model.Gamma[t])/(model.Xc[t]+2*model.Gamma[t]) ==     model.epsilon[t] 
    
    model.epsilon_constraint = Constraint(range(T), rule=epsilon_constraint_rule)
    
    def sigmaco2_constraint_rule(model, t):
        return 1/model.sigmaco2[t] == 1/cbnd+1/cstm+1/(ccar1*model.Xt[t]*model.Xt[t]+                                                      ccar2*model.Xt[t]+ccar3)
    model.sigmaco2_constraint = Constraint(range(T), rule=sigmaco2_constraint_rule)
    
    def phi_photomax_constraint_rule(model, t):
        return model.phi_photomax[t] ==     model.epsilon[t]*cpar*(0.52*model.Ul[t])*    model.sigmaco2[t]*(model.Xc[t]-model.Gamma[t])/(model.epsilon[t]*cpar*                                                    (0.52*model.Ul[t])+                                                    model.sigmaco2[t]*(model.Xc[t]-model.Gamma[t]))
    
    model.phi_photomax_constraint = Constraint(range(T), rule=phi_photomax_constraint_rule)
    
    def phi_photo_constraint_rule(model, t):
        return model.phi_photo[t] == model.phi_photomax[t]*(1-exp(-ck*clars*(1-ctau)*model.Xs[t]))
    
    model.phi_photo_constraint = Constraint(range(T), rule=phi_photo_constraint_rule)
    
    def phi_resp_constraint_rule(model, t):
        return model.phi_resp[t] ==     (cresps*(1-ctau)+crespr*ctau)*model.Xs[t]*cQ10resp**((model.Xt[t]-25)/10)
    model.phi_resp_constraint = Constraint(range(T), rule=phi_resp_constraint_rule)
    
    def rgr_constraint_rule(model, t):
        return model.rgr[t] ==     crgrmax*model.Xn[t]/(model.Xs[t]+model.Xn[t])*cQ10gr**((model.Xt[t]-20)/10)
    model.rgr_constraint = Constraint(range(T), rule=rgr_constraint_rule)
    
    def phi_calpl_constraint_rule(model, t):
        return model.phi_calpl[t] == model.phi_photo[t] - model.phi_resp[t]/calpha -     (1-cbeta)/(calpha*cbeta)*model.rgr[t]*model.Xs[t]
    model.phi_calpl_constraint = Constraint(range(T), rule=phi_calpl_constraint_rule)
    
    ### canopy transpiration
    def phi_hplal_constraint_rule(model, t):
        return model.phi_hplal[t] ==     (1-exp(-capl*model.Xs[t]))*cvplal*(cv0*cv1*cH20/(cr*(model.Xt[t]+ctabs))*exp(cv2*model.Xt[t]/(model.Xt[t]+cv3))-model.Xh[t])
    model.phi_hplal_constraint = Constraint(range(T), rule=phi_hplal_constraint_rule)
    
    
    ## from absolute humidity to relative humidity
    def Rh_constraint_rule(model, t):
        return model.Rh[t]*(cv1*cH20/(cr*(model.Xt[t]+ctabs))*exp(cv2*model.Xt[t]/(cv3+model.Xt[t]))) ==     100*model.Xh[t]
    model.Rh_constraint = Constraint(range(T), rule=Rh_constraint_rule)
    
    def temp_dynamics_constraint_rule(model, t):
        return ccapq*(model.Xt[t+1]-model.Xt[t]) ==     3600*(model.Uqs[t]-model.Uqr[t]+0.48*model.Ul[t]-    c_U*(c_surface_area/c_grow_area)*(model.Xt[t]-Vout[t])-    model.Uv[t]*1290*(model.Xt[t]-Vout[t])-c_lat_water*1000*model.phi_hplal[t]) 
    model.temp_dynamics_constraint = Constraint(range(T-1), rule=temp_dynamics_constraint_rule)
    
    def co2_dynamics_constraint_rule(model, t):
        return model.Xc[t+1]-model.Xc[t] ==     3600*(model.Uc[t]-model.phi_calpl[t]-model.Uv[t]*(model.Xc[t]-Cout[t]))*c_grow_area/c_volume
    model.co2_dynamics_constraint = Constraint(range(T-1), rule=co2_dynamics_constraint_rule)
    
    def humidity_dynamics_constraint_rule(model, t):
        return model.Xh[t+1]-model.Xh[t] ==     3600*(model.phi_hplal[t] - model.Uv[t]*(model.Xh[t]-Hout[t])-     model.Uhr[t])*c_grow_area/c_volume
    model.humidity_dynamics_constraint = Constraint(range(T-1), rule=humidity_dynamics_constraint_rule)
    
    def Xn_dynamics_constraint_rule(model, t):
        return model.Xn[t+1]-model.Xn[t] == 3600*(calpha*model.phi_photo[t] -                                                   model.rgr[t]*model.Xs[t]-    model.phi_resp[t]-(1-cbeta)/cbeta*model.rgr[t]*model.Xs[t])
    model.Xn_dynamics_constraint = Constraint(range(T-1), rule=Xn_dynamics_constraint_rule)
    
    def Xs_dynamics_constraint_rule(model, t):
        return model.Xs[t+1]-model.Xs[t] == 3600*model.rgr[t]*model.Xs[t]
    model.Xs_dynamics_constraint = Constraint(range(T-1), rule=Xs_dynamics_constraint_rule)
    
    # constraint initialize Xt
    model.Xt_initial = Constraint(expr=model.Xt[0] == XT0)
    # constraint initialize Xh
    model.Xh_initial = Constraint(expr=model.Xh[0] == XH0)
    #constraint initialize Xc
    model.Xc_initial = Constraint(expr=model.Xc[0] == XC0*0.04401/24.45/10**3)
    # constraint initialize Xn
    model.Xn_initial = Constraint(expr=model.Xn[0] == XN0)
    # constraint initialize Xs
    model.Xs_initial = Constraint(expr=model.Xs[0] == XS0)
    
    def Objfunc(model):
        return c_grow_area*(3600*pc*sum(model.Uc[t] for t in range(T))+    pe*sum(model.Ul[t]/1000+3600*(model.Uhr[t])/c_dehum_eev+           model.Uv[t]/fan_eff+(model.Uqs[t]+model.Uqr[t])/(1000*3) for t in range(T)))

    
    model.obj = Objective(rule=Objfunc, sense=minimize)
    #model.pprint()
    solver = SolverFactory("ipopt")
    solver.options['acceptable_tol'] = 1e-2  # Relax feasibility requirement
    solver.options['acceptable_iter'] = 50  # Allow some infeasibility
    solver.solve(model)
    
    

    print("Cost:", model.obj()) 

        
    print("Xt:", min([model.Xt[t].value for t in range(T)]),             max([model.Xt[t].value for t in range(T)]))
    print("Xc in ppm:", min([model.Xc[t].value*(24.45/0.04401)*1000 for t in range(T)]),             max([model.Xc[t].value*(24.45/0.04401)*1000 for t in range(T)]))
    print("Rh:", min([model.Rh[t].value for t in range(T)]),              max([model.Rh[t].value for t in range(T)]))
        
    print("Optimal solution:")
    print("Uqs:", min([model.Uqs[t].value for t in range(T)]),              max([model.Uqs[t].value for t in range(T)]))
    print("Uqr:", min([model.Uqr[t].value for t in range(T)]),              max([model.Uqr[t].value for t in range(T)]))
    print("Uhr:", min([model.Uhr[t].value for t in range(T)]),              max([model.Uhr[t].value for t in range(T)]))
    print("Uc:", min([model.Uc[t].value for t in range(T)]),              max([model.Uc[t].value for t in range(T)]))
    print("Ul:", min([model.Ul[t].value for t in range(T)]),              max([model.Ul[t].value for t in range(T)]))
    print("Uv:", min([model.Uv[t].value for t in range(T)]),              max([model.Uv[t].value for t in range(T)]))
 
        
    return model.obj(),np.array([c_grow_area/3000*(model.Uqs[t].value) for t in range(T)]),    np.array([c_grow_area/3000*(model.Uqr[t].value) for t in range(T)]),    np.array([c_grow_area*3600/c_dehum_eev*(model.Uhr[t].value) for t in range(T)]),    np.array([c_grow_area*3600*model.Uc[t].value for t in range(T)]),    np.array([c_grow_area/1000*model.Ul[t].value for t in range(T)]),    np.array([c_grow_area/fan_eff*model.Uv[t].value for t in range(T)]),    np.array([model.Xn[t].value for t in range(T)]),    np.array([model.Xs[t].value for t in range(T)]),    np.array([model.Xt[t].value for t in range(T)]),    np.array([model.Xc[t].value*(24.45/0.04401)*1000 for t in range(T)]),    np.array([model.Rh[t].value for t in range(T)]),    np.array([model.Xh[t].value for t in range(T)])
    

    

    
    
    
    


# In[ ]:




