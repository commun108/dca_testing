# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 15:32:36 2014

@author: commun

Script to filter tomographic TDIP data on the analysis of the decay curve.

CHANGELOG:
    
    
TODO:
    
"""
#import of needed modules
import numpy as np
import matplotlib.pyplot as plt
import pylab
import matplotlib.gridspec as gridspec
import scipy
import matplotlib



#close open figures
plt.close('all')


########################################################################
#    FUNCTION DEFINITION                                               #
########################################################################

def slope(data, error, log=None):
    """
    Filtering for positive slope. 
    
    Measurements with positive slope of the linear fit model are filtered
    based on the assumption that curves with positive values represent
    extensive noise.
    """
    ind = np.ones(len(data), dtype=bool)
    ind[error[:,8]>0] = 0
    data_f, error_f = data[ind], error[ind]
    
    #plot filtered percentage
    print(' SLOPE '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))
    
    if log is not None:
        f = log
        f.write(' SLOPE '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))

    return data_f, error_f
    
def negcha(data, error, log=None):
    """
    Filtering for negative or equal zero chargeabilities.
    
    Measurements are removed where the integral chargeability is below
    or equal zeros.
    """
    ind = np.ones(len(data), dtype=bool)
    ind[data[:,5]<=0] = 0
    data_f, error_f = data[ind], error[ind]
    
    #plot filtered percentage
    print(' NEGCHA '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))
    
    if log is not None:
        f = log
        f.write(' NEGCHA '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))
    
    return data_f, error_f
    
def meanrms(data, error, log=None):
    """
    Filtering for mean RMS of fit to measured window chargeabilities and
    fit only on the even and odd gate numbers. 
    
    By looking at the histogram of the misfit values gaps are detected
    and values above the gap are rejected. A higher mean RMS between the
    fits can be seen as a measure of noise in the decay curve.
    """
    dev = error[:,2]
    dev = dev[dev<200] 
    nhist = 20
        
    v, ed = np.histogram(dev, nhist)
    for idx, call in enumerate(v):
        if call == 0:
            rule =  error[:,2]<ed[idx]
            data_f = data[rule]
            error_f = error[rule]
            break
            
    #plot filtered percentage
    print(' MEANRMS '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))
    
    if log is not None:
        f = log
        f.write(' MEANRMS '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))
    
    return data_f, error_f
    
def zerogate(data, error, log=None):
    """
    Filtering for the number of gates with a chargeability value below 
    zeros.
    
    The number of gates with a chargeability below zero are determined. 
    If the number exceeds or is equal to 6 the measurement is filtered.
    """
    ind = np.ones(len(data), dtype=bool)

    for line in range(len(data)):
        s = np.sum(data[line,9:29]<0)
        if s>=6:
            ind[line] = 0          
            
    data_f, error_f = data[ind], error[ind]

    #plot filtered percentage
    print(' ZEROGATE '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))
    
    if log is not None:
        f = log
        f.write(' ZEROGATE '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))
    
    return data_f, error_f
    
def mastercurves(data, error, log=None):
    """
    Filtering for the deviation of a mastercurve to the measured curve.
    
    Iterative methodology. The histogramm is calculated and is searched for
    gaps. Values above the gap are rejected. Then the next iteration is 
    run. Should the first gap has a values below 1.5 the process is stopped.
    """
    dev = error[:,11]
    dev = dev[dev<25]
    nhist = 30
    tresh = 1.5
    data_f = data.copy()
    error_f = error.copy()
    
    for it in range(5):
        
        v, ed = np.histogram(dev, nhist)
        for idx, call in enumerate(v):
            #~ if call == 0:
            if call <= 1:
                if ed[idx]<tresh and it>0:
                    continue
                else:    
                    rule =  error_f[:,11]<ed[idx]
                    rule2 = dev<ed[idx]
                    dev = dev[rule2]
                    data_f = data_f[rule]
                    error_f = error_f[rule]
                    
    #plot filtered percentage
    print(' MASTERCURVES '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))               

    if log is not None:
        f = log
        f.write(' MASTERCURVES '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))  

    return data_f, error_f
      
def rmsmisfit(data, error, log=None):
    """
    Filtering for RMS between fit and measured decay of every gate.
    
    As the RMS can be seen as a value of noisyness of the decay curve, 
    curves which exceed a treshold are filtered. Eg only curves with low
    noise should be kept. 
    """    
    
    data_f = data.copy()
    error_f = error.copy()
    
    for it in range(5):
        ind = np.ones(len(data_f), dtype=bool)
        
        if np.percentile(error_f[:,7], 90)>1:
            ind[error_f[:,7]>np.percentile(error_f[:,7], 90)] = 0
            
        else:
            if max(error_f[:,7])<1.5:
                break
            else:    
                ind[error_f[:,7]>np.percentile(error_f[:,7], 99.5)] = 0
                
        data_f, error_f = data_f[ind], error_f[ind]  
        
    #plot filtered percentage
    print(' RMSMISFIT '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))               

    if log is not None:
        f = log
        f.write(' RMSMISFIT '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))  

    return data_f, error_f

def chagap_coarse(data, error, log=None):
    """
    Filtering for gaps in the histogram of chargeabilities.
    
    Iterative process. Filtering detects gaps in the histograms and removes
    values after gap. The number of bins is set dynamically.
    """  
    data_f = data.copy()
    error_f = error.copy()
      
    for it in range(5):
    
        mx = np.max(data_f[:,5])*1.25
        nbins = np.ceil(mx)
        
        v, ed = np.histogram(data_f[:,5], nbins)
        tresh = 2
                
        for idx, call in enumerate(v):
            if call < tresh:
                if ed[idx]<1.5 and it>0:
                    break
                if ed[idx]<np.median(data_f[:,5]):
                    continue
                else:  
                    rule =  data_f[:,5]<ed[idx]
                    data_f = data_f[rule]
                    error_f = error_f[rule]
    
    #plot filtered percentage
    print(' CHAGAP - COARSE '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))               

    if log is not None:
        f = log
        f.write(' CHAGAP - COARSE '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100)) 

    return data_f, error_f
    
def chagap_dense(data, error, log=None):
    """
    Filtering for gaps in the histogram of chargeabilities.
    
    Iterative process. Filtering detects gaps in the histograms and removes
    values after gap. The number of bins is set dynamically.
    """  
    data_f = data.copy()
    error_f = error.copy()
      
    for it in range(10):
    
        mx = len(data_f[:,5])*0.05
        nbins = np.ceil(mx)
        #~ print(nbins)
        
        #~ plt.figure(it)
        #~ plt.hist(data_f[:,5], nbins)
        #~ plt.show()
        
        v, ed = np.histogram(data_f[:,5], nbins)
        #~ print(v)
        tresh = 1
                
        for idx, call in enumerate(v):
            if call < tresh:
                if ed[idx]<np.median(data_f[:,5]):
                    continue
                else:  
                    rule =  data_f[:,5]<ed[idx]
                    data_f = data_f[rule]
                    error_f = error_f[rule]
    
    #plot filtered percentage
    print(' CHAGAP - DENSE '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))               
    
    if log is not None:
        f = log
        f.write(' CHAGAP - DENSE '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))  

    return data_f, error_f

def chagap_rev(data, error, log=None):
    """
    Filtering for gaps in the histogram of chargeabilities.
    
    Iterative process. Filtering detects gaps in the histograms and removes
    values after gap. The number of bins is set dynamically.
    """  
    data_f = data.copy()
    error_f = error.copy()
      
    for it in range(0,2):
    
        mx = len(data_f[:,5])*0.05
        nbins = np.ceil(mx)
        
        v, ed = np.histogram(data_f[:,5], nbins)
        v = np.fliplr([v])[0]
        ed = np.fliplr([ed])[0]

        tresh = 1
                
        for idx, call in enumerate(v):
            if call < tresh:
                if ed[idx]>np.median(data_f[:,5]):
                    continue
                else:
                    rule =  data_f[:,5]>ed[idx]
                    data_f = data_f[rule]
                    error_f = error_f[rule]

    #plot filtered percentage
    print(' CHAGAP - REV '.center(40,'='))
    print('# Measurements: %d' % (len(data)))
    print('# Remaining measurements: %d' % (len(data_f)))
    print('# Measurements filtered: %d' % (len(data) - len(data_f)))
    perc = 1 - len(data_f)/len(data)
    print('Percentage filtered: %2.2f\n' % (perc*100))               

    if log is not None:
        f = log
        f.write(' CHAGAP - REV '.center(40,'='))
        f.write('\n')
        f.write('# Measurements: %d\n' % (len(data)))
        f.write('# Remaining measurements: %d\n' % (len(data_f)))
        f.write('# Measurements filtered: %d\n' % (len(data) - len(data_f)))
        f.write('Percentage filtered: %2.2f\n' % (perc*100))  

    return data_f, error_f
    
########################################################################
#    FILTERING                                                         #
########################################################################

def filt(path, path_analysis, folder_out, filtering, path_figures):
    
    frags = path.split('/')
    lid = frags[-1][:-4]
    
    data = np.genfromtxt(path, delimiter='\t',skiprows=1)
    error = np.genfromtxt(path_analysis + lid +'_error_t.dat',delimiter='\t',skiprows=1)

    mdelay = data[0,29]
    ipw = (mdelay+np.cumsum(data[0,30:50]))-data[0,30:50]/2
    
    path_filt = folder_out
    
    f = open(path_analysis + lid + '_' + filtering + 'b.log', 'w')
    
    data_f1, error_f1 = slope(data, error, f)
    data_f2, error_f2 = negcha(data_f1, error_f1, f)
    data_f3, error_f3 = meanrms(data_f2, error_f2, f)
    data_f4, error_f4 = zerogate(data_f3, error_f3, f)
    data_f5, error_f5 = mastercurves(data_f4, error_f4, f)
    data_f6, error_f6 = rmsmisfit(data_f5, error_f5, f)
    
    #~ data_f6, error_f6 = chagap_coarse(data_f6, error_f6, f)
    data_f7, error_f7 = chagap_dense(data_f5, error_f5, f)
    
    data_f8, error_f8 = chagap_rev(data_f7, error_f7, f)
    
    #remove unused misfit values for writing
    error_f8 = np.delete(error_f8, (3, 4, 5, 8, 9), 1)

########################################################################
#    WRITE FILES                                                       #  
########################################################################

    np.savetxt(path_filt + lid + '_' + filtering + 'b.dat', data_f8, delimiter='\t',
    fmt='%d\t%d\t%d\t%d\t%f\t%.4f\t%.1f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t'+
        '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' + 
        '\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' +
        '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d'
        ,
        header='A	B	M	N	R	    	Y		devR	V		I		Mi  mdelay  IPW stack',
    comments='') 
    
    np.savetxt(path_filt + lid + '_' + filtering + 'b.er', error_f8, delimiter='\t',
    fmt='%3.6f',
    header='rms1	rms2	mean_rms	rmsfm   rms_misfit  d_res   d_pha_sf    d_pha_af' +
    '   rmsmm  ipw_misfit ipw_misfit_mm',
    comments='') 

########################################################################
#    PLOT DECAYS AND DISTRIBUTION (INT CHARG)                          #  
########################################################################

    print(' SUMMARY '.center(40,'='))
    print('Number of measurements: %d' % (len(data)))
    print('Number of remaining measurements: %d' % (len(data_f8)))
    print('Number of filtered measurements: %d' % (len(data) - len(data_f8)))
    perc = 1 - len(data_f8)/len(data)
    print('Percentage filtered: %2.2f' % (perc*100))
    
    #~ with open(path_analysis + lid + '_' + filtering + 'b.log', 'w') as f:
        #~ f.write(' SUMMARY '.center(40,'='))
        #~ f.write('\n')
        #~ f.write('Number of measurements: %d\n' % (len(data)))
        #~ f.write('Number of remaining measurements: %d\n' % (len(data_f7)))
        #~ f.write('Number of filtered measurements: %d\n' % (len(data) - len(data_f7)))
        #~ perc = 1 - len(data_f7)/len(data)
        #~ f.write('Percentage filtered: %2.2f\n' % (perc*100))
        
    
    f.write(' SUMMARY '.center(40,'='))
    f.write('\n')
    f.write('Number of measurements: %d\n' % (len(data)))
    f.write('Number of remaining measurements: %d\n' % (len(data_f8)))
    f.write('Number of filtered measurements: %d\n' % (len(data) - len(data_f8)))
    perc = 1 - len(data_f8)/len(data)
    f.write('Percentage filtered: %2.2f\n' % (perc*100))
    f.close()
    
    plt.ioff() 
    plt.figure(0)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized() 

    gs = gridspec.GridSpec(2, 3)
    plt.subplot(gs[0:2,0:2])
    for line in range(len(data_f8)):
        plt.plot(ipw, data_f8[line,9:29])
    plt.grid('on')
    plt.xlabel('Time [ms]')
    plt.ylabel('Chargeability [mV/V]')
    #~ plt.xlim([200, 1900])
    
    plt.subplot(gs[:,2])
    plt.hist(data_f8[:,5], bins=15)
    plt.ylabel('# of Measurements')
    plt.xlabel('Chargeability [mV/V]')
    plt.tight_layout()
    pylab.savefig(path_figures + lid + '_' + filtering + '_sum.png', bbox_inches='tight')
    plt.close()
    #~ plt.show()

