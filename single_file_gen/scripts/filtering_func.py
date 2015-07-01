# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 15:32:36 2014

@author: commun

Script to filter tomographic TDIP data on the analysis of the decay curve.
3 models are fitted to the sampled decay. A fit of these which is lower than
a specified treshold are dropped, the global chargeability is set to nan.
Preliminary filtering is done computing linear regression models and 
rejecting measurements at certain parameters.

CHANGELOG:
23/01/15:
	- np.loadtxt -> genfromtxt, more versatile, stable
	- move the filtering process in own loop, prior to that compute statistical paramters
	- add linear regression slope to stat loop
03/02/15:    
    - add condition to check if fit is not a line, if so use pow2m
    - catch expection in case alternative pow2m throws runtimeerror
    - move filtering to seperate file
    - implement cross correlation between fit on data and mean decay
23/02/15:
    - changed folder paths
11/03/15:
    - implemented deviation of resistance (misfit between measured resistance
    and modeled one!)
18/03/15:
    - implemented misfit of measured global chargeability and global 
    chargeability of fit
    - minor improvements of fitting routine
25/03/15
    - output file with curve fit parameters (for usability)
    - change computation of mean decay -> preliminary filtering is changed
    - make fitting on measurements for misfit mean decay/measured decay
    more robust (exception!)
    - file with window-wise misfit of measured decay and fitted decay
2/04/15
    - catch a exception when fitting of first model doesn't work -> set all
    error values to 9999 or else
30/04/15
    - implementation of masters approach -> derivation of 3 master curves valid 
    for data set, rms between the nearest master curve and measured decay is minimized
"""

#import of needed modules
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.optimize import curve_fit
import scipy
import minimize as mz
from scipy import stats

#defining the used fitting models
#linear
def linear(x, a, b):
    return a*x+b
#pow2
def pow2(x, a, b, c):
    return a*x**-b+c
#pow2m    
def pow2m(x, a, b, c):
    return a*x**b+x**-c 


def filtering(path, folder_out, path_figures):

########################################################################
#   FILENAME
########################################################################
    
    
    
    
    path_filt = folder_out #relative path to folder with slash at end
    
    #loading the data set
    data = np.genfromtxt(path, delimiter='\t', skip_header=1)
    
    #setting the ip window gating
    mdelay = data[0,29]
    gates = data[0,30:50]
    ipw = (mdelay+np.cumsum(data[0,30:50]))-data[0,30:50]/2
    

    #starting values for coefficients of models
    p0_linear = [1, 1]
    p0_pow2 = [1, 0, 1]
    #vectors for rms values
    rms_1, rms_2 = np.ones((len(data),1)), np.ones((len(data),1))
    #vectors for xcorr values
    xc1, xc2 = np.ones((len(data),1)), np.ones((len(data),1))
    #vector for dev(fit/mean) and dev(measured/mean)
    rmsfm, rmsmm = np.zeros((len(data),1)), np.zeros((len(data),1))
    #vector for rms misfit between data (mi) and fit
    rms_misfit = np.ones((len(data),1))
    #array for linear regression parameter (slope)
    linrg = np.ones((len(data),1))
    #array for rms value between measured decay and master decay
    rmsmmaster = np.ones((len(data),1))
    #array for rms values between fit on measured decay and master decay
    rmsfmaster = np.ones((len(data),1))
    #array for resistance misfit
    dev_res = np.zeros((len(data), 1))
    #array for chargeability misfit (single fit)
    dev_pha_sf = np.zeros((len(data), 1))
    #array for chargeability misfit (all fit)
    dev_pha_af = np.zeros((len(data), 1))
    
    #array for curve fit parameters: 3 parameters per fit and 3 fits + 
    #1 identifier wich model was fitted (0 -> pow2 / 1 -> pow2m)
    fit_param = np.zeros((len(data), 10))
    
    #array for window-wise misfit of measured and fitted decay
    #20 windows -> 20 columns
    ipw_misfit = np.zeros((len(data), 20))
    
    #array for window-wise misfit of mean decay and measured decay (minimized)
    ipw_misfit_mm = np.zeros((len(data), 20))
    
   

    #loop over data set to derive error parameters and do the fitting
    for line in range(len(data)):
        mi = data[line,9:29]
        
        try:            
            #fit using all sample points
            p, c = curve_fit(pow2, ipw, mi, p0_pow2, maxfev=100000)
            #fit using only even points
            pe, ce = curve_fit(pow2, ipw[0::2], mi[0::2], p0_pow2, maxfev=100000)
            #fit using only odd points
            po, co = curve_fit(pow2, ipw[1::2], mi[1::2], p0_pow2, maxfev=100000)
            
            if ((p[0]>0 and p[1]>0 and p[2]<0) and
                (pe[0]>0 and pe[1]>0 and pe[2]<0) and
                (po[0]>0 and po[1]>0 and po[2]<0)):
                    
                #fitting with pow2 is ok -> eval fit on window mids
                f = pow2(ipw, p[0], p[1], p[2])
                fe = pow2(ipw, pe[0], pe[1], pe[2])
                fo = pow2(ipw, po[0], po[1], po[2])
                #fitting for resistance misfit (model voltage)
                vm = pow2(7, p[0], p[1], p[2])*100
                #write the fit parameters to file
                fit_param[line, 0:3] = p
                fit_param[line, 3:6] = pe
                fit_param[line, 6:9] = po
            else:
                try:
                        
                    #fitting pow2m model
                    #use whole data for fitting
                    p, c = curve_fit(pow2m, ipw, mi, p0_pow2, maxfev=100000)
                    #only use even points
                    pe, ce = curve_fit(pow2m, ipw[0::2], mi[0::2], p0_pow2, maxfev=100000)
                    #only use odd points
                    po, co = curve_fit(pow2m, ipw[1::2], mi[1::2], p0_pow2, maxfev=100000)
                                    
                except RuntimeError:
                    p = p0_pow2
                    pe = p0_pow2
                    po = p0_pow2
                    
                f = pow2m(ipw, p[0], p[1], p[2])
                fe = pow2m(ipw, pe[0], pe[1], pe[2])
                fo = pow2m(ipw, po[0], po[1], po[2])
                #fitting for resistance misfit (model voltage)
                vm = pow2m(7, p[0], p[1], p[2])*100
                #write the fit parameters to file
                fit_param[line, 0:3] = p
                fit_param[line, 3:6] = pe
                fit_param[line, 6:9] = po
                fit_param[line, 9] = 1
           
    
            #compute rmse and cross correlation between fit1 and fit2, fit1 and fit3
            
            #misfit between fit and measured data (goodness of fit)
            rms_misfit[line] = np.sqrt(np.mean((mi-f)**2))
                
            #normalize parameters to get xcorr in range -1, 1
            fn = (f - np.mean(f)) / (np.std(f) * len(f))
            fen = (fe - np.mean(fe)) /  np.std(fe)
            fon = (fo - np.mean(fo)) /  np.std(fo)
                    
            #assign xcor values to vector  
            xc1[line] = max(np.correlate(fn, fen, 'full'))
            xc2[line] = max(np.correlate(fn, fon, 'full'))
            
            #assign rmse values to vector
            rms_1[line] = np.sqrt(np.mean((f - fe)**2))
            rms_2[line] = np.sqrt(np.mean((f - fo)**2))
            
            #compute linear regression parameters, only use pl[0] (slope)
            pl, cl = curve_fit(linear, ipw, mi, p0_linear, maxfev=100000) 
            linrg[line] = pl[0] 
            
########################################################################
# COMPUTE DEVIATION OF RESISTANCE FOR ~T=0 OF MODELED VOLTAGE (FIT) 
######################################################################## 
    
            #get the voltage in vicinity to t=0 (*100 is needed)
            #evaluate the fit function for the given parameters
            # -> is done in fitting section
            
            #get the current of the measurement
            current = data[line, 8]
            
            #compute the resistance for the fit
            resm = vm/current
            
            #calculate the misfit between measured resistance and the modeled one
            res_real = data[line, 4]
            dev_res[line] = res_real - resm
        
########################################################################
# COMPUTE MISFIT OF GLOBAL CHARGEABILITY 
######################################################################## 
    
            #measured global chargeability
            gc_meas = data[line, 5]
            
            #compute global chargeability of fit
            gc_mod_sf = np.mean(f)
            gc_mod_af = (np.mean(f) + np.mean(fo) + np.mean(fe))/3
            
            #get the misfit by subtracting measured and modelled values
            dev_pha_sf[line] = gc_meas - gc_mod_sf
            dev_pha_af[line] = gc_meas - gc_mod_af

########################################################################
# GET WINDOW-WISE MISFIT (MEASURED/FITTED) 
######################################################################## 
    
            ipw_misfit[line] = mi - f
    
######################################################################## 
    #if even the fitting of a robust model isn't possible set all the errors 
    #to nonsense values
        except RuntimeError:
            fit_param[line, 0:3] = 9999
            fit_param[line, 3:6] = 9999
            fit_param[line, 6:9] = 9999
            
            rms_misfit[line] = 9999
            
            xc1[line] = 0
            xc2[line] = 0
        
            rms_1[line] = 9999
            rms_2[line] = 9999
            
            linrg[line] = 9999
            
            dev_res[line] = 9999
            
            dev_pha_sf[line] = 9999
            dev_pha_af[line] = 9999
            
            ipw_misfit[line] = 9999
        
########################################################################
########################################################################

########################################################################
# COMPUTE THE MEAN DECAY OF DATA SET
########################################################################

    #in order to compute a mean decay of the data set the set is filtered 
    #before the calculation
    
    #array of filter indices (boolean)
    ind = np.ones((len(data),1), dtype=bool)
    
    #filter for slope
    ind[linrg>0] = 0
    
    #filter for negative chargeabilities
    ind[data[:,5]<=0] = 0
    
    #filter for rms misfit
    ind[rms_misfit>0.04] = 0
    
    #apply filter
    data_c = data[ind[:,0]]
    
    a, b = np.shape(data_c)
    
    if a <= 20:
        
        ind = np.ones((len(data),1), dtype=bool)
    
        #filter for slope
        ind[linrg>0] = 0
        
        #filter for negative chargeabilities
        ind[data[:,5]<=0] = 0
        
        #filter for rms misfit
        ind[rms_misfit>np.percentile(rms_misfit, 25)] = 0
        
        #apply filter
        data_c = data[ind[:,0]]
        
    
    #compute 3 master curves of data set
    masters = np.zeros((3, 20))
    medge = np.zeros((4, 20))
    n_b = 3
    
    #bin the gate values of all curves at every gate
    for gate in range(20):
        bmd, bed, c = stats.binned_statistic(np.sort(data_c[:,9+gate], axis=0),
        np.sort(data_c[:,9+gate], axis=0), statistic=np.median, bins=n_b)
        masters[:,gate] = bmd
        medge[:,gate] = bed
    
    #compute the integral chargeability of the master curves -> used for next steps   
    masters_ints = np.zeros((3,1))
    for ints in range(3):
        masters_ints[ints] = np.mean(masters[ints,:])
    
    #with filtered data, compute mean decay curve of data set
    m_mean = np.zeros(20)
    ms = data_c[:,9:29]
    for ll in range(len(m_mean)):
        m_mean[ll] = np.median(ms[:,ll])
            
    #compute deviation of single decay to mean decay, before calculate fit    
    for line in range(len(data)):
        mi = data[line,9:29]    
        if fit_param[line, -1] == 0:
            f = pow2(ipw, fit_param[line, 0],  fit_param[line, 1], fit_param[line, 2])
        else:
            f = pow2m(ipw, fit_param[line, 0],  fit_param[line, 1], fit_param[line, 2])
        
            
        #compute rms between mean decay and fit on data    
        rmsfm[line], nn = mz.minimize(m_mean, ipw, f)
    
        #compute rms between mean decay and measured data  
        rmsmm[line], ipw_misfit_mm[line] = mz.minimize(m_mean, ipw, mi)
        
        #compute distances of measured intregal chargeability to the master curves
        #in order to find the nearest master curve
        dists = np.zeros((3,1))
        for dist in range(3):
            dists[dist] = abs(np.mean(mi)-masters_ints[dist])
            
        #get index of shortest distance    
        idx = np.argmin(dists) 
        
        #compute rms between measured decay and nearest master curve
        rmsmmaster[line], x = mz.minimize(masters[idx], ipw, mi)
        
        #compute rms between fit on measured decay and nearest master curve
        rmsfmaster[line], x = mz.minimize(masters[idx], ipw, f)
            
    #storing rms/deviation values
    mean_xc = (xc1 + xc2)/2
    mean_rms = (rms_1 + rms_2)/2
    error = np.concatenate( \
    (rms_1, rms_2, mean_rms, xc1, xc2, mean_xc, rmsfm, rms_misfit,
    linrg, rmsmmaster, dev_res, rmsfmaster, dev_pha_af, rmsmm, ipw_misfit, ipw_misfit_mm),
    axis=1)


    frags = path.split('/')
    lid = frags[-1][:-4]
    
    #write error parameters to file
    np.savetxt(path_filt + lid + '_error_t.dat',
    error, fmt='%3.6f',
    delimiter='\t',
    header='rms1 rms2 mean_rms xc1 xc2 mean_xc rmsfm rms_misfit slope_lrg rmsmmaster dev_res rmsfmaster dev_pha_af rmsmm ipw_misfit ipw_misfit_mm',
    comments='')    
    
    #write fit parameters to file
    np.savetxt(path_filt + lid + '_fit_param.dat',
    fit_param, fmt='%3.6f',
    delimiter='\t',
    header='p_1    p_2    p_3    pe_1    pe_2    pe_3    po_1    po_2    po_3    id',
    comments='') 
    
    #plotting
    plt.ioff()
    fig = plt.figure(figsize=(8, 5))
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    
    for line in range(len(data)):
        plt.plot(ipw, data[line, 9:29], 'b', linewidth=0.5)
        plt.hold('True')
    for line in range(len(data_c)):
        plt.plot(ipw, data_c[line, 9:29], 'k', linewidth=1)
    for iter in range(3):
        plt.plot(ipw, masters[iter, :], 'r', linewidth=2)
    
    plt.ylim([0, np.max(masters)+5])
    plt.xlim([ipw[0]-50, ipw[-1]+50])
    plt.xlabel('Time [ms]', fontsize=14)
    plt.ylabel('Apparent Chargability [mV/V]', fontsize=14)
    plt.title('Mastercurves (red), strictly filtered data (black)', fontsize=16)
    
    fig.savefig(path_figures + lid + '_masters.png', bbox_inches='tight', dpi=200)
    plt.close()
    
 
    
