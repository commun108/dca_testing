# Import the needed modules
import numpy as np

def read_ipwin_func(path, folder_out):
    """
    Read the data from the syscal output file (*.txt),
    -> generate a output file which contains:
    A	B	M	N	Resistance [Ohm]	Integral Chargeability [mV/V]	devR []
    Voltage [mV]	Current [mA]	Window Chargeability [mV/V] mdelay [ms]
    IP Window Size [ms] Stack []
    -> the value names will be abbreviated (# of used columns beneath):
    A	B	M	N	R	Y	devR	V	I	Mi	mdelay  IPW Stack
    1	1	1	1	1	1	1	1	1	20	1   20  1
    Jakob Gallist, modified after Adrian Flores-Orozco
    """
    
    frags = path.split('/')
    identifier = frags[-1].find('r')
    
    if identifier == -1:
        #read as normal
        
        
    
        # Setting some variables
        dist = 1 	# distance between electrodes, can be kept at 1
        corr = 0 	# offset in electrode 1
            
        # Reading the data
        din = np.genfromtxt(path, delimiter=' ', skiprows=1, usecols=range(0,53))
        
        # Removing unused columns
        din = np.delete(din,(4,10),1)
        #~ din = np.delete(din,(-4,-3,-2,-1),1)
        
        # Check if current is in an appropriate range, otherwise remove row
        din = din[np.logical_not(np.logical_and(din[:,8] < 0, din[:,8]))] # current
        din = din[np.logical_not(np.logical_and(din[:,8] > 99999, din[:,8]))]
        
        # Loop to switch the potential electrodes
        for x in range(len(din)):
            if din[x,7] < 0: 						# check for negative voltage
                din[x,0] = 1+din[x,0]/dist-corr 	# A
                din[x,1] = 1+din[x,1]/dist-corr		# B
                din[x,3] = 1+din[x,2]/dist-corr		# M
                din[x,2] = 1+din[x,3]/dist-corr		# N
                din[x,7] = -din[x,7]				# change voltage to positiv
            else:
                din[x,0] = 1+din[x,0]/dist-corr 	# A
                din[x,1] = 1+din[x,1]/dist-corr		# B
                din[x,2] = 1+din[x,2]/dist-corr		# M
                din[x,3] = 1+din[x,3]/dist-corr		# N	
        
        
        din[:,6] = din[:,4]
        
        # Calculating the resistance R = V/I
        din[:,4] = din[:,7]/din[:,8]
        #~ din = din[:,0:29]
        
        #path manipulation to save file in right folder
        path_out = folder_out
        lid = frags[-1][:-4]
        
        
        # Writing the data to a .dat file
        np.savetxt(path_out + lid + '.dat', din, 
        fmt='%d\t%d\t%d\t%d\t%f\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t'+
        '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' + 
        '\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' +
        '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d'
        ,
        header='A	B	M	N	R		    Y		devR	V		I		Mi  mdelay  IPW stack',
        comments='')
        
        
        
    else:
    
        #read as reciprocal
        
        # Setting some variables
        dist = 1 	# distance between electrodes, can be kept at 1
        corr = 0 	# offset in electrode 1
        
    
        
        # Reading the data
        din = np.genfromtxt(path, delimiter=' ', skiprows=1, usecols=range(0,53))
        
        nelec = np.max(din[:,0:4])+1
        
        # Removing unused columns
        din = np.delete(din,(4,10),1)
        #~ din = np.delete(din,(-4,-3,-2,-1),1)
        
        # Check if current is in an appropriate range, otherwise remove row
        din = din[np.logical_not(np.logical_and(din[:,8] < 0, din[:,8]))] # current
        din = din[np.logical_not(np.logical_and(din[:,8] > 99999, din[:,8]))]
        
        #file for new values
        m,n = np.shape(din)
        #~ don = np.zeros((m, n))
        don = din.copy()
        
        # Loop to switch the potential electrodes
        for x in range(len(din)):
            if din[x,7] < 0: 						# check for negative voltage
                don[x,0] = nelec-(din[x,0]/dist-corr) 	# A
                don[x,1] = nelec-(din[x,1]/dist-corr)		# B
                don[x,3] = nelec-(din[x,2].copy()/dist-corr)		# M
                don[x,2] = nelec-(din[x,3].copy()/dist-corr)		# N
                don[x,7] = -din[x,7]				# change voltage to positiv
            else:
                don[x,0] = nelec-(din[x,0]/dist-corr) 	# A
                don[x,1] = nelec-(din[x,1]/dist-corr)		# B
                don[x,2] = nelec-(din[x,2]/dist-corr)		# M
                don[x,3] = nelec-(din[x,3]/dist-corr)		# N	
        
        don[:,6] = din[:,4]

        # Calculating the resistance R = V/I
        
        don[:,4] = don[:,7]/din[:,8]
        #~ don = don[:,0:29]
        
        #~ don = din.copy()
        idx = np.lexsort((don[:,3], don[:,2], don[:,1], don[:,0]))
        don = don[idx]
        #~ don = don[np.lexsort((don[:,1], don[:,0]))]
        
        #path manipulation to save file in right folder
        path_out = folder_out
                
        lid = frags[-1][:-4]
    
        
        # Writing the data to a .dat file
        np.savetxt(path_out + lid + '.dat', don, 
        fmt='%d\t%d\t%d\t%d\t%f\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t'+
        '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' + 
        '\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' +
        '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d'
        ,
        header='A	B	M	N	R	    	Y		devR	V		I		Mi  mdelay  IPW stack',
        comments='')
