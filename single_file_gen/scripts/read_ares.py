

# Import the needed modules
import numpy as np
import linecache

def read_ares(path, folder_out):
    """
    Read the data from the ares output file (*.usd),
    -> generate a output file which contains:
    A	B	M	N	Resistance [Ohm]	Integral Chargeability [mV/V]	devR []
    Voltage [mV]	Current [mA]	Window Chargeability [mV/V] mdelay [ms]
    IP Window Size [ms] Stack []
    -> the value names will be abbreviated (# of used columns beneath):
    A	B	M	N	R	Y	devR	V	I	Mi	mdelay  IPW Stack   Midev
    1	1	1	1	1	1	1	1	1	20	1   20  1   20
    Jakob Gallist, modified after Adrian Flores-Orozco, Matthias Bücker
    """
    
    path_out = folder_out
    
    # Reading the data
    din = np.genfromtxt(path, delimiter='\t', skip_header=29)
    
    # Read the gates and compute ipw
    gates = linecache.getline(path, 20)
    gates = gates.split('\t')
    gates[-1] = gates[-1].strip(' ms\n')
    gates = np.array(gates[1:]).astype(float)
    
    # Mdelay and ipw
    mdelay = 5 
    ipw = (mdelay + np.cumsum(gates))-gates/2
    linecache.clearcache()

    # Read stacking parameters
    stack = linecache.getline(path, 25)
    stack = stack.strip('Stacking: ' + '\n')
    stack = stack.replace(' ', '')
    # only use minimum stack information
    stack = float(stack[0])

    don = din.copy()
    
    don[:,0:2] = np.fliplr(din[:,0:2])+1
    don[:,2:4] = np.fliplr(din[:,2:4])+1
    
    din[:,0:4] = din[:,0:4]+1
    
    for line in range(len(din)):
        if din[line,9]<0:
            don[line,2] = din[line,2]
            don[line,3] = din[line,3]
            don[line,9] = -din[line,9]
    
    #transfer resistance
    don[:,4] = don[:,9]/din[:,8]
    
    #deviation resitance
    don[:,6] = din[:,12]*10
    
    #voltage
    don[:,7] = don[:,9]
    
    #current
    don[:,8] = din[:,8]
    
    #window chargeabilities
    mi = din[:,13::2]*10
    
    #deviation per window
    devmi = din[:,14::2]
    
    #integral chargeability
    for line in range(len(din)):
        don[line,5] = np.sum(np.multiply(mi[line,:], gates))/np.sum(gates)
    
    gts = gates.reshape(1,20)
    gts = gts*np.ones((len(don),1))
    
    #pack all values in one array
    dout = np.concatenate((don[:,0:9], mi, np.ones((len(don),1))*mdelay, gts, np.ones((len(don),1))*stack, devmi), axis=1)
    
    frags = path.split('/')
    lid = frags[-1][:-4]
    
    # Writing the %data to a .dat file
    np.savetxt(path_out + lid + '.dat', dout, 
    fmt='%d\t%d\t%d\t%d\t%f\t%.4f\t%.1f\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t'+
    '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' + 
    '\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' +
    '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d' +
    '\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' +
    '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f',
    header='A	B	M	N	R		Y		devR		V		I		Mi  mdelay  IPW stack devMi',
    comments='')







		
