#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
def from_relative_humidity(xT, xH):
    cv1=611
    cv2=17.4
    cv3=239
    cH20=18
    ctabs=273
    cr=8314
    
    return xH/100*(cv1*cH20/(cr*(xT+ctabs))*np.exp(cv2*xT/(xT+cv3)))

def continuous_opt(DVout,DCout,DHout,D=28,Xt0=20,Xh0=from_relative_humidity(20,0.63),                   Xc0=356,Xn0=0.003,Xs0=0.003,T=24):
    
    Eqsx=np.zeros((D,T))
    Eqrx=np.zeros((D,T))
    Ehsx=np.zeros((D,T))
    Ehrx=np.zeros((D,T))
    Ecx=np.zeros((D,T))
    Elx=np.zeros((D,T))
    Evx=np.zeros((D,T))
    Xnx=np.zeros((D,T))
    Xsx=np.zeros((D,T))
    Xtx=np.zeros((D,T))
    Xcx=np.zeros((D,T))
    Rhx=np.zeros((D,T))
    Xhx=np.zeros((D,T))
    cost=np.zeros((D,))
    
    for d in range(D):
        print(d)
        if d == 0:
            cost[d,],Eqsx[d,],Eqrx[d,],Ehrx[d,],Ecx[d,],Elx[d,],Evx[d,],            Xnx[d,],Xsx[d,],Xtx[d,],Xcx[d,],Rhx[d,],Xhx[d,] =             PFALopt(Vout=DVout[d,:],Cout=DCout[d,:],Hout=DHout[d,:],                  XT0=Xt0,XH0=Xh0,XC0=Xc0,XN0=Xn0,XS0=Xs0)
        else:
            cost[d,],Eqsx[d,],Eqrx[d,],Ehrx[d,],Ecx[d,],Elx[d,],Evx[d,],            Xnx[d,],Xsx[d,],Xtx[d,],Xcx[d,],Rhx[d,],Xhx[d,] =             PFALopt(Vout=DVout[d,:],Cout=DCout[d,:],Hout=DHout[d,:],                  XT0=Xtx[d-1,T-1],XH0=Xhx[d-1,T-1],XC0=Xcx[d-1,T-1],XN0=Xnx[d-1,T-1],                   XS0=Xsx[d-1,T-1])
    return cost,Eqsx,Eqrx,Ehrx,Ecx,Elx,Evx,Xnx,Xsx,Xtx,Xcx,Rhx,Xhx        
      
    


# In[ ]:




