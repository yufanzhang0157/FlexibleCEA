#!/usr/bin/env python
# coding: utf-8

# In[1]:


def cost_baseline_cal(baseline,weather,c_grow_area=80,pe_up=0.15,pe_down=0.1,pe=0.2051,pc=0.2,                 D=28,T=24):
    
    costEqs={}
    costEqr={}
    costEhs={}
    costEhr={}
    costC={}
    costEl={}
    costEv={}
    costE={}
    
    
    
    cost={}

    
    FW = {}
    
    Elec_capacity = {}
    
    for key in weather.keys():
        costEqs[key] = np.sum(baseline[key]['Eqs']*pe)/c_grow_area
        costEqr[key] = np.sum(baseline[key]['Eqr']*pe)/c_grow_area
        costEhr[key] = np.sum(baseline[key]['Ehr']*pe)/c_grow_area
        costC[key] = np.sum(baseline[key]['Ec']*pc)/c_grow_area
        costEl[key] = np.sum(baseline[key]['El']*pe)/c_grow_area
        costEv[key] = np.sum(baseline[key]['Ev']*pe)/c_grow_area
        costE[key] = costEqs[key]+costEqr[key]+costEhr[key]+costEl[key]+costEv[key]
        
        Elec_capacity[key] = np.max(baseline[key]['Eqs']+baseline[key]['Eqr']+baseline[key]['Ehr']+                                   baseline[key]['El']+baseline[key]['Ev'])/c_grow_area
        
        FW[key] =         (baseline[key]['Xn'].reshape(-1)[-1]+baseline[key]['Xs'].reshape(-1)[-1])*22.5*(1-0.07)
    
    
    print('Annual Heating cost per m2:',sum(costEqs.values()))
    print('Annual Cooling cost per m2:',sum(costEqr.values()))
    print('Annual Dehumidification cost per m2:',sum(costEhr.values()))
    print('Annual Lighting cost per m2:',sum(costEl.values()))
    print('Annual Ventilation cost per m2:',sum(costEv.values()))
    print('Annual CO2 cost per m2:',sum(costC.values()))
    print('Annual Electricity cost per m2:',sum(costE.values()))
    
    
    print()
    print('Annual net cost per m2:',sum(costE.values())+sum(costC.values()))
    print('FW kg',sum(FW.values()))
    print('Annual net cost/per kg FW:',(sum(costE.values())+sum(costC.values()))/sum(FW.values()))
    
    print()
    print('Mamimum electricity consumption per m2',max(Elec_capacity.values()))
    
    


# In[ ]:




