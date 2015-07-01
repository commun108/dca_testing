# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 15:32:36 2014

@author: commun

Converts filtered data sets to tx2 format used in Aarhusinv. 

CHANGELOG:
    
"""
#import of needed modules
import numpy as np

#~ def to_tx2(path_filt, path_tx2, electrode_spacing, ngates=20, pulselength=2):

#~ path_filt = '../shiprock_filt/v3/l6sk0n_1_1b.dat'
path_filt = '../shiprock_filt/v4/l8_1_mg_4b.dat'
#~ path_tx2 = '../shiprock_tx2/'
path_tx2 = '../'
electrode_spacing = 2
pulselength=2
ngates=20

frags = path_filt.split('/')
lid = frags[-1][:-4]

ngates = ngates
elec_sp = electrode_spacing
pl = pulselength
path_filt = path_filt
path_tx2 = path_tx2

filt = np.genfromtxt(path_filt, delimiter='\t', skip_header=1)
    
########################################################################
#    DATA WRANGLING                                                    #
########################################################################

out = np.zeros((len(filt), 28+4*ngates))

#x-electrode position along measuring tape (m)
out[:,0:4] = (filt[:,0:4]-1)*elec_sp

#x-position of electrodes
#eg use of gps coordinates if available
out[:,4:8] = (filt[:,0:4]-1)*elec_sp

#y-position of electrodes
#~ out[:,8:12]

#z-position of electrodes
#~ out[:,12:16]

#apparent resistivity


a = (filt[:,0]-1)*elec_sp
b = (filt[:,1]-1)*elec_sp
m = (filt[:,2]-1)*elec_sp
n = (filt[:,3]-1)*elec_sp

am = np.sqrt((m-a)**2)
an = np.sqrt((n-a)**2)
bm = np.sqrt((m-b)**2)
bn = np.sqrt((n-b)**2)

k = 2*np.pi*(1/(1/am-1/an-1/bm+1/bn))

#~ out[:,16] = filt[:,4]
out[:,16] = k[:]*filt[:,4]

#deviation resistance 
out[:,17] = filt[:,6]/10

#resistance flag (0 keeps data, 1 removes data)
#~ out[:,18]

#number of ip gates
out[:,19] = ngates

#ip values [mV/V] per gate    
out[:,20:20+ngates] = filt[:,9:29]

#mdelay [ms]
#~ out[:,20+ngates] = filt[0,29]
out[:,20+ngates] = filt[0,29]


#gate lengths [ms]
#~ out[:,21+ngates:21+ngates*2] = filt[:,30:50]
out[:,21+ngates:21+ngates*2] = filt[:,30:50]

#deviation of every window
#for syscal files put 0.0 until proper error model is introduced
#for ares: use values given by device

#check for device
c, d = np.shape(filt)

if d>51:
    out[:,21+ngates*2:21+ngates*3] = filt[:,51:]
else:
    out[:,21+ngates*2:21+ngates*3] = 0.1

#ip flag
#~ out[:,21+ngats*3:21+ngates*4]

#stacking
out[:,21+ngates*4] = filt[:,50]
    
#current = (50% of pulse length)
out[:,22+ngates*4] = pl/2

#wave type
out[:,23+ngates*4] = pl/2

#Ton = pulse length [s]
out[:,24+ngates*4] = pl*1000

#Toff = pulse length [s]
out[:,25+ngates*4] = pl*1000

#Tend = time used to collect decay curve [s]
out[:,26+ngates*4] = (np.sum(filt[0,30:50])+filt[0,29])

#Tstart = time at which the decay curve is starting to be measured
out[:,27+ngates*4] = filt[0,29]
    
#######################################################################
#   WRITE TO FILE                                                     #
#######################################################################

M = []
G = []
STD = []
IP_flag = []
Mfmt = []
Gfmt = []
STDfmt = []
IP_flagfmt = []

for num in range(ngates):
    M.append('M' + str(num+1))
    Mfmt.append('%.3f\t')
    G.append('Gate' + str(num+1))
    Gfmt.append('%.3f\t')
    STD.append('Std' + str(num+1))
    STDfmt.append('%.2f\t')
    IP_flag.append('IP_Flag' + str(num+1))
    IP_flagfmt.append('%d\t')

#~ print(M)

M = '\t'.join(M)
G = '\t'.join(G)
STD = '\t'.join(STD)
IP_flag = '\t'.join(IP_flag)

Mfmt = ''.join(Mfmt)
Gfmt = ''.join(Gfmt)
STDfmt = ''.join(STDfmt)
IP_flagfmt = ''.join(IP_flagfmt)


np.savetxt(path_tx2 + lid + '.tx2',
out, 
#~ fmt='%3.6f',
fmt='%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t' + #electrode positions
'%f\t%.3f\t%d\t%d\t' + #resistance, devR, ResFlag, Ngates
Mfmt + '%.3f\t' + Gfmt + STDfmt + IP_flagfmt + 
'%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f',
delimiter='\t',
header='xA	xB	xM	xN	UTMxA	UTMxB	UTMxM	UTMxN	UTMyA	UTMyB\t' +
'UTMyM	UTMyN	zA	zB	zM	zN	Res	Dev	ResFlag	Ngates\t' + M + '\tMdly\t' + G + '\t' + STD +
'\t' + IP_flag + '\tStack	Current\tWaveType	Ton	Toff	Tend	Tstart',
comments='') 
