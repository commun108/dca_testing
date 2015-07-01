# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 15:32:36 2014

@author: commun

Python module with usefull functions, eg minimizer
"""

import numpy as np

def minimize(m_mean, ipw, fit):
    """
    fit is a at the ip windows ipw evaluated function, hence an array
    m_mean is the mean decay curve of the data set, the y values have to 
    be on the same position as ipw
    """
     
    #mean values of decay -> global chargeability 
    #used for determination of up- or downshift
    mg = np.mean(m_mean)
    ac = np.mean(fit)
    shift=0.001 #shift value
    
    #downshift
    if ac>mg:
        delta_i = m_mean-fit #compute difference between mean and fit
        rmse_i = sum(delta_i**2)/len(delta_i) #rmse 
        cfs = fit-shift #shift the fit downwards
        
        for kk in range(10**9): #loop for minimization
            delta = m_mean-cfs #new difference
            rmse = sum(delta**2)/len(delta) #new rmse
            if rmse>rmse_i:
               rmse_f = rmse_i 
               break
            else:
               shift += 0.001 #shift again and iterate
               cfs = cfs-shift
               rmse_i = rmse
    #upshift    
    else:  
        delta_i = m_mean-fit
        rmse_i = sum(delta_i**2)/len(delta_i)
        cfs = fit+shift #upshift
        
        for kk in range(10**9):
           delta = m_mean-cfs
           rmse = sum(delta**2)/len(delta)
           if rmse>rmse_i:
               rmse_f = rmse_i
               break
           else:
               shift += 0.001
               cfs = cfs+shift
               rmse_i = rmse
    
    return rmse_f, delta         
        
        
        
