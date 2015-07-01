# -*- coding: utf-8 -*-
"""
Script to display raw and filtered data in pseudosections.
"""

#import of needed modules
import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys
from matplotlib import rcParams
from matplotlib import ticker

#close open figures
plt.close('all')

def plot_pseudo(path, path_filt, filtering, path_figures, elec_spacing=1):
    
    elec_sp = elec_spacing
    
    frags = path.split('/')
    lid = frags[-1][:-4]
    ident = frags[-1].find('mg')
   
    
    if ident == -1: # dipole-dipole data set - plot dd pseudo section
        
        ########################################################################
        #    LOADING PARAMETERS AND DATA                                       #
        ########################################################################
        
        
        #raw 
        raw = np.genfromtxt(path, delimiter='\t', skip_header=1)
        #filtered
        filt = np.genfromtxt(path_filt + lid + '_' + filtering + 'b.dat',
        delimiter='\t', skip_header=1)
                               
        #colorbar limits
        #data limits (chargeability)
        ymin, ymax = 0, 15
        #marker size -> needs to be fitted to screen size
        ms = 80
        
        ########################################################################
        #    GET PLOT POSITIONS                                                #
        ########################################################################
        
        raw[:,0:4] = (raw[:,0:4]-1)*elec_sp
        filt[:,0:4] = (filt[:,0:4]-1)*elec_sp
        
        #### RAW ####
        
        #vector for z-component
        zd_raw = np.zeros(len(raw))
        x_raw = np.zeros(len(raw))
        
        #compute x-component
        x_raw = ((raw[:,0]-1 + raw[:,1]-1)/2 + (raw[:,2]-1 + raw[:,3]-1)/2) / 2
        #compute z-component 
        zd_raw = -(x_raw-(raw[:,0]-1 + raw[:,1]-1)/2) 
        
               
        #### FILTERED ####
        
        #vector for z-component
        zd_filt = np.zeros(len(filt))
        x_filt = np.zeros(len(filt))
        
        #compute x-component
        x_filt = ((filt[:,0]-1 + filt[:,1]-1)/2 + (filt[:,2]-1 + filt[:,3]-1)/2) / 2
        #~ #compute z-component 
        zd_filt = -(x_filt-(filt[:,0]-1 + filt[:,1]-1)/2) 

        zd_raw = np.abs(zd_raw)*-1
        zd_filt = np.abs(zd_filt)*-1
        
             
        ########################################################################
        #    PLOTTING                                                          #
        ########################################################################
        
        #### RAW AND FILTERED - SUBPLOT 1 ####
        
        #maximize plot window
        plt.ioff()
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=17)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        #~ plt.tick_params(axis='y', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + ' - raw data', fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(x_raw, zd_raw,
        s=ms, #size of marker
        c=raw[:,5], #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=17) #label of colorbar
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=18)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + '_' + filtering + 'b' + ' - filtered data', fontsize=14)
        #~ plt.title('decay curve analysis filtered', fontsize=18)
        
        #filtered data - measurements
        scat2 = plt.scatter(x_filt, zd_filt,
        s=ms, #size of marker
        c=filt[:,5], #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax, #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=13) #label of colorbar
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b.png', bbox_inches='tight', dpi=200)
        
        
        plt.close()
        #~ plt.show()

        
        
    else: #multiple gradient data set - plot mg pseudo section
        
        ########################################################################
        #    LOADING PARAMETERS AND DATA                                       #
        ########################################################################

        #raw 
        raw = np.genfromtxt(path, delimiter='\t', skip_header=1)
        #filtered
        filt = np.genfromtxt(path_filt + lid + '_' + filtering + 'b.dat',
        delimiter='\t', skiprows=1)

        #colorbar limits
        #data limits (chargeability)
        ymin, ymax = 0, 15
        #size of marker
        sizem = 80
        size2 = 40
        size3 = 22
        size4 = 10

        ########################################################################
        #    GET PLOT POSITIONS                                                #
        ########################################################################
        
        raw[:,0:4] = (raw[:,0:4]-1)*elec_sp
        filt[:,0:4] = (filt[:,0:4]-1)*elec_sp
        
        #### RAW  ####
        
        #vector for z-component
        zd_raw = np.zeros(len(raw))
        
        #compute z-component
        x_raw = (raw[:,3] + raw[:,2])/2
        zmna = x_raw - raw[:,0]
        zmnb = raw[:,1] - x_raw
        
        #derive z-component 
        I = (zmna<zmnb)
        zd_raw[I] = -zmna[I]/3
        I = (zd_raw==0)
        zd_raw[I] = -zmnb[I]/3
        
        zd_raw = np.abs(zd_raw)*-1
        
        x_raw = np.expand_dims(x_raw, axis=1)
        zd_raw = np.expand_dims(zd_raw, axis=1)
        raws = np.expand_dims(raw[:,5], axis=1)
        B = np.concatenate((x_raw, zd_raw, raws), axis=1)
        
        B = B[np.lexsort((B[:,1], B[:,0]))]
        
        
        ms = np.ones(len(raw))*sizem
        
        for l in range(1,len(raw)):
            if ms[l]==size3 or ms[l]==size4:
                continue
            if B[l,0]==B[l-1,0] and B[l,1]==B[l-1,1]:
                ms[l]=size2
                if l <= len(raw)-2 and B[l+1,0]==B[l,0] and B[l+1,1]==B[l,1]:
                    ms[l+1]=size3
                    if l <= len(raw)-3 and B[l+2,0]==B[l,0] and B[l+2,1]==B[l,1]:
                        ms[l+2]=size4
        
        #### FILTERED ####
        
        #vector for z-component
        zd_filt = np.zeros(len(filt))
        
        #compute z-component
        x_filt = (filt[:,3] + filt[:,2])/2
        zmna2 = x_filt - filt[:,0]
        zmnb2 = filt[:,1] - x_filt
        
        #derive z-component 
        I = (zmna2<zmnb2)
        zd_filt[I] = -zmna2[I]/3
        I = (zd_filt==0)
        zd_filt[I] = -zmnb2[I]/3
        
        zd_filt = np.abs(zd_filt)*-1
        
        x_filt = np.expand_dims(x_filt, axis=1)
        zd_filt = np.expand_dims(zd_filt, axis=1)
        filts = np.expand_dims(filt[:,5], axis=1)
        B2 = np.concatenate((x_filt, zd_filt, filts), axis=1)
        
        B2 = B2[np.lexsort((B2[:,1], B2[:,0]))]
        
        
        msf = np.ones(len(filt))*sizem
        
        for l in range(1,len(filt)):
            if msf[l]==size3 or msf[l]==size4:
                continue
            if B2[l,0]==B2[l-1,0] and B2[l,1]==B2[l-1,1]:
                msf[l]=size2
                if l <= len(filt)-2 and B2[l+1,0]==B2[l,0] and B2[l+1,1]==B2[l,1]:
                    msf[l+1]=size3
                    if l <= len(filt)-3 and B2[l+2,0]==B2[l,0] and B2[l+2,1]==B2[l,1]:
                        msf[l+2]=size4
        

        ########################################################################
        #    PLOTTING                                                          #
        ########################################################################

        #### RAW AND FILTERED - SUBPLOT 1 ####
        
        #maximize plot window
        plt.figure(0)
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - raw data: ' + lid, fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(B[:,0], B[:,1],
        s=ms, #size of marker
        c=B[:,2],
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - DCA filtered data: ' + lid + '_' + filtering + 'b', fontsize=14)
        
        #filtered data - measurements
        scat2 = plt.scatter(B2[:,0], B2[:,1],
        s=msf, #size of marker
        c=B2[:,2],
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax) #set colormap limits
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b.png', bbox_inches='tight', dpi=200)
        plt.close('all')






def plot_pseudo_plus_resistivity(path, path_filt, filtering, path_figures, elec_spacing=1):
    
    elec_sp = elec_spacing
    
    frags = path.split('/')
    lid = frags[-1][:-4]
    ident = frags[-1].find('mg')
   
    
    if ident == -1: # dipole-dipole data set - plot dd pseudo section
        
        ########################################################################
        #    LOADING PARAMETERS AND DATA                                       #
        ########################################################################
        
        
        #raw 
        raw = np.genfromtxt(path, delimiter='\t', skip_header=1)
        #filtered
        filt = np.genfromtxt(path_filt + lid + '_' + filtering + 'b.dat',
        delimiter='\t', skip_header=1)
                               
        #colorbar limits
        #data limits (chargeability)
        ymin, ymax = 0, 15
        #marker size -> needs to be fitted to screen size
        ms = 80
        
        ########################################################################
        #    GET PLOT POSITIONS                                                #
        ########################################################################
        
        raw[:,0:4] = (raw[:,0:4]-1)*elec_sp
        filt[:,0:4] = (filt[:,0:4]-1)*elec_sp
        
        #### RAW ####
        
        #vector for z-component
        zd_raw = np.zeros(len(raw))
        x_raw = np.zeros(len(raw))
        
        #compute x-component
        x_raw = ((raw[:,0]-1 + raw[:,1]-1)/2 + (raw[:,2]-1 + raw[:,3]-1)/2) / 2
        #compute z-component 
        zd_raw = -(x_raw-(raw[:,0]-1 + raw[:,1]-1)/2) 
        
        a = raw[:,0]
        b = raw[:,1]
        m = raw[:,2]
        n = raw[:,3]
        
        am = np.sqrt((m-a)**2)
        an = np.sqrt((n-a)**2)
        bm = np.sqrt((m-b)**2)
        bn = np.sqrt((n-b)**2)

        k_raw = 2*np.pi*(1/(1/am-1/an-1/bm+1/bn))
        
               
        #### FILTERED ####
        
        #vector for z-component
        zd_filt = np.zeros(len(filt))
        x_filt = np.zeros(len(filt))
        
        #compute x-component
        x_filt = ((filt[:,0]-1 + filt[:,1]-1)/2 + (filt[:,2]-1 + filt[:,3]-1)/2) / 2
        #~ #compute z-component 
        zd_filt = -(x_filt-(filt[:,0]-1 + filt[:,1]-1)/2) 

        zd_raw = np.abs(zd_raw)*-1
        zd_filt = np.abs(zd_filt)*-1
        
        a = filt[:,0]
        b = filt[:,1]
        m = filt[:,2]
        n = filt[:,3]
        
        am = np.sqrt((m-a)**2)
        an = np.sqrt((n-a)**2)
        bm = np.sqrt((m-b)**2)
        bn = np.sqrt((n-b)**2)

        k_filt = 2*np.pi*(1/(1/am-1/an-1/bm+1/bn))
             
        ########################################################################
        #    PLOTTING                                                          #
        ########################################################################
        
        #### RAW AND FILTERED - SUBPLOT 1 ####
        
        #maximize plot window
        plt.ioff()
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=17)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        #~ plt.tick_params(axis='y', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + ' - raw data', fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(x_raw, zd_raw,
        s=ms, #size of marker
        c=raw[:,5], #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=17) #label of colorbar
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=18)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + '_' + filtering + 'b' + ' - filtered data', fontsize=14)
        #~ plt.title('decay curve analysis filtered', fontsize=18)
        
        #filtered data - measurements
        scat2 = plt.scatter(x_filt, zd_filt,
        s=ms, #size of marker
        c=filt[:,5], #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax, #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=13) #label of colorbar
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b.png', bbox_inches='tight', dpi=200)
        
        
        plt.close()
        #~ plt.show()
        
        #maximize plot window
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=17)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        #~ plt.tick_params(axis='y', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + ' - raw data', fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(x_raw, zd_raw,
        s=ms, #size of marker
        c=np.log10(np.abs(raw[:,4])*k_raw), #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        #~ vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=17) #label of colorbar
        c.set_label('apparent rho [Ohmm]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=18)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Pseudosection - Dipole-Dipole: ' + lid + '_' + filtering + 'b' + ' - filtered data', fontsize=14)
        #~ plt.title('decay curve analysis filtered', fontsize=18)
        
        #filtered data - measurements
        scat2 = plt.scatter(x_filt, zd_filt,
        s=ms, #size of marker
        c=np.log10(np.abs(filt[:,4])*k_filt), #values for marker (chargeability values)
        marker='s', #type of marker -> square
        linewidth=0.5,
        #~ edgecolors='none',
        cmap='jet', #set colormap
        #~ vmin=ymin, vmax=ymax, #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        #~ c.set_label('Apparent Chargeability m [mV/V]', fontsize=13) #label of colorbar
        c.set_label('apparent rho [Ohmm]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b_rho.png', bbox_inches='tight', dpi=200)
        
        
        plt.close()
        
        
    else: #multiple gradient data set - plot mg pseudo section
        
        ########################################################################
        #    LOADING PARAMETERS AND DATA                                       #
        ########################################################################

        #raw 
        raw = np.genfromtxt(path, delimiter='\t', skip_header=1)
        #filtered
        filt = np.genfromtxt(path_filt + lid + '_' + filtering + 'b.dat',
        delimiter='\t', skiprows=1)

        #colorbar limits
        #data limits (chargeability)
        ymin, ymax = 0, 15
        #size of marker
        sizem = 80
        size2 = 40
        size3 = 22
        size4 = 10

        ########################################################################
        #    GET PLOT POSITIONS                                                #
        ########################################################################
        
        raw[:,0:4] = (raw[:,0:4]-1)*elec_sp
        filt[:,0:4] = (filt[:,0:4]-1)*elec_sp
        
        a = raw[:,0]
        b = raw[:,1]
        m = raw[:,2]
        n = raw[:,3]
        
        am = np.sqrt((m-a)**2)
        an = np.sqrt((n-a)**2)
        bm = np.sqrt((m-b)**2)
        bn = np.sqrt((n-b)**2)

        k_raw = 2*np.pi*(1/(1/am-1/an-1/bm+1/bn))
        
        rho_raw = np.abs(raw[:,4])*k_raw
        
        #### RAW  ####
        
        #vector for z-component
        zd_raw = np.zeros(len(raw))
        
        #compute z-component
        x_raw = (raw[:,3] + raw[:,2])/2
        zmna = x_raw - raw[:,0]
        zmnb = raw[:,1] - x_raw
        
        #derive z-component 
        I = (zmna<zmnb)
        zd_raw[I] = -zmna[I]/3
        I = (zd_raw==0)
        zd_raw[I] = -zmnb[I]/3
        
        zd_raw = np.abs(zd_raw)*-1
        
        x_raw = np.expand_dims(x_raw, axis=1)
        zd_raw = np.expand_dims(zd_raw, axis=1)
        raws = np.expand_dims(raw[:,5], axis=1)
        raws_rho = np.expand_dims(rho_raw, axis=1)
        B = np.concatenate((x_raw, zd_raw, raws, raws_rho), axis=1)
        
        B = B[np.lexsort((B[:,1], B[:,0]))]
        
        
        ms = np.ones(len(raw))*sizem
        
        for l in range(1,len(raw)):
            if ms[l]==size3 or ms[l]==size4:
                continue
            if B[l,0]==B[l-1,0] and B[l,1]==B[l-1,1]:
                ms[l]=size2
                if l <= len(raw)-2 and B[l+1,0]==B[l,0] and B[l+1,1]==B[l,1]:
                    ms[l+1]=size3
                    if l <= len(raw)-3 and B[l+2,0]==B[l,0] and B[l+2,1]==B[l,1]:
                        ms[l+2]=size4
        
        #### FILTERED ####
        
        a = filt[:,0]
        b = filt[:,1]
        m = filt[:,2]
        n = filt[:,3]
        
        am = np.sqrt((m-a)**2)
        an = np.sqrt((n-a)**2)
        bm = np.sqrt((m-b)**2)
        bn = np.sqrt((n-b)**2)

        k_filt = 2*np.pi*(1/(1/am-1/an-1/bm+1/bn))
        
        rho_filt = np.abs(filt[:,4])*k_filt
        
        #vector for z-component
        zd_filt = np.zeros(len(filt))
        
        #compute z-component
        x_filt = (filt[:,3] + filt[:,2])/2
        zmna2 = x_filt - filt[:,0]
        zmnb2 = filt[:,1] - x_filt
        
        #derive z-component 
        I = (zmna2<zmnb2)
        zd_filt[I] = -zmna2[I]/3
        I = (zd_filt==0)
        zd_filt[I] = -zmnb2[I]/3
        
        zd_filt = np.abs(zd_filt)*-1
        
        x_filt = np.expand_dims(x_filt, axis=1)
        zd_filt = np.expand_dims(zd_filt, axis=1)
        filts = np.expand_dims(filt[:,5], axis=1)
        filts_rho = np.expand_dims(rho_filt, axis=1)
        B2 = np.concatenate((x_filt, zd_filt, filts, filts_rho), axis=1)
        
        B2 = B2[np.lexsort((B2[:,1], B2[:,0]))]
        
        
        msf = np.ones(len(filt))*sizem
        
        for l in range(1,len(filt)):
            if msf[l]==size3 or msf[l]==size4:
                continue
            if B2[l,0]==B2[l-1,0] and B2[l,1]==B2[l-1,1]:
                msf[l]=size2
                if l <= len(filt)-2 and B2[l+1,0]==B2[l,0] and B2[l+1,1]==B2[l,1]:
                    msf[l+1]=size3
                    if l <= len(filt)-3 and B2[l+2,0]==B2[l,0] and B2[l+2,1]==B2[l,1]:
                        msf[l+2]=size4
        

        ########################################################################
        #    PLOTTING                                                          #
        ########################################################################

        #### RAW AND FILTERED - SUBPLOT 1 ####
        
        #maximize plot window
        plt.figure(0)
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - raw data: ' + lid, fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(B[:,0], B[:,1],
        s=ms, #size of marker
        c=B[:,2],
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - DCA filtered data: ' + lid + '_' + filtering + 'b', fontsize=14)
        
        #filtered data - measurements
        scat2 = plt.scatter(B2[:,0], B2[:,1],
        s=msf, #size of marker
        c=B2[:,2],
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet', #set colormap
        vmin=ymin, vmax=ymax) #set colormap limits
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        c.set_label('apparent m [mV/V]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b.png', bbox_inches='tight', dpi=200)
        plt.close('all')
        
        # RESISTIVITY
        
        
        #maximize plot window
        plt.figure(0)
        figure = plt.figure(figsize=(8, 5))
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        
        #raw data - electrodes
        plt.subplot(2,1,1)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - raw data: ' + lid, fontsize=14)
        
        #raw data - measurements
        scat1 = plt.scatter(B[:,0], B[:,1],
        s=ms, #size of marker
        c=np.log10(B[:,3]),
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet', #set colormap
        #~ vmin=ymin, vmax=ymax,  #set colormap limits
        alpha=1) #opacity
        c = plt.colorbar(scat1, orientation='vertical', pad=0.01)
        c.set_label('apparent rho [Ohmm]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        #filtered data - electrodes
        plt.subplot(2,1,2)
        #~ plt.plot(elec[:,1], elec[:,2], 'ko')
        plt.hold('True')
        plt.xlim([min(x_raw)-1*elec_sp, max(x_raw)+1*elec_sp])
        plt.ylim([min(zd_raw)-0.5, max(zd_raw)+0.5])
        plt.xlabel('x [m]', fontsize=14)
        #~ plt.ylabel('z [m]', fontsize=14)
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        plt.tick_params(axis='x', which='both', labelsize=12)
        plt.title('Multiple gradient - DCA filtered data: ' + lid + '_' + filtering + 'b', fontsize=14)
        
        #filtered data - measurements
        scat2 = plt.scatter(B2[:,0], B2[:,1],
        s=msf, #size of marker
        c=np.log10(B2[:,3]),
        marker='o', #type of marker -> circle
        edgecolors='k',
        linewidth=0.3,
        cmap='jet') #set colormap
        #~ vmin=ymin, vmax=ymax) #set colormap limits
        c = plt.colorbar(scat2, orientation='vertical', pad=0.01)
        c.set_label('apparent rho [Ohmm]', fontsize=14) #label of colorbar
        c.ax.tick_params(labelsize=12) 
        tick_locator = ticker.MaxNLocator(nbins=7)
        c.locator = tick_locator
        c.update_ticks()
        
        
        plt.tight_layout()
        figure.savefig(path_figures + lid + '_' + filtering + 'b_rho.png', bbox_inches='tight', dpi=200)
        plt.close('all')
