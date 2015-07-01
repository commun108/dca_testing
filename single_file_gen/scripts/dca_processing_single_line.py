"""
DCA processing

This routine reads syscal or ares raw data, converts to work format,
computes decay curve analysis, does automatic filtering, plotting of raw
and filtered pseudosections and conversion to tx2 for inversion.
and plotting. 
"""

# Import the needed modules
import glob
import read_ipwin_gen as read_ip
import read_ipwin_gen_full_spreadsheet_input as rd
#~ import read_ares 
import filtering_func as filtering
import filt_func_comb as filtcomb
import filt_func_dense as filtdense
import filt_func_coarse as filtcoarse
import plot_pseudo as pp
import to_tx2_func 


########################################################################
#   READING RAW DATA / CONVERSION
########################################################################

## SYSCAL ##

# Get array with paths to files
lid = 'l6_1_mg'
path_raw = '../raw/new_txt/' + lid + '.txt'
path_output = '../raw/new_dat/'

# Calling read_ipwin_func 
#~ read_ip.read_ipwin_func(path_raw, path_output)
rd.read_ipwin_func(path_raw, path_output)
    
## ARES ##

# Get array with paths to files
#~ path_raw = '../shiprock_raw/new_usd/l1sk0n_1.usd'
#~ path_output = '../shiprock_raw/new_dat/'

# Calling read_ares
#~ read_ares.read_ares(path_raw, path_output)

########################################################################
#   DECAY CURVE ANALYSIS
########################################################################

path_raw = '../raw/new_dat/' + lid + '.dat'
path_figures = '../figures/'
path_output = '../filt/analysis/' # slash at end is important

filtering.filtering(path_raw, path_output, path_figures)

########################################################################
#   OUTLIER REMOVAL
########################################################################

path_analysis = '../filt/analysis/'
path_filt = '../filt/v1/'
path_figures = '../figures/'
filtering = '1' # needs to be a string

filtcomb.filt(path_raw, path_analysis, path_filt, filtering, path_figures)
#~ filtcoarse.filt(path_raw, path_analysis, path_filt, filtering, path_figures)
#~ filtdense.filt(path_raw, path_analysis, path_filt, filtering, path_figures)
    
########################################################################
#   PLOT OF PSEUDOSECTION
########################################################################

#for correct x-positions of electrodes set elec_spacing to desired value
#plot_pseudo(path to raw data, path to filtered data, filter number
#path to figure folder, optional: elec_spacing (real electrode spacing, defaul=1))

#~ pp.plot_pseudo(path_raw, path_filt, filtering, path_figures)
#~ pp.plot_pseudo(path_raw, path_filt, filtering, path_figures, 2)
pp.plot_pseudo_plus_resistivity(path_raw, path_filt, filtering, path_figures, 2)
    

########################################################################
#   CONVERT TO TX2
########################################################################

filtering = '1'
path_filt = '../filt/v1/' + lid + '_' + filtering + 'b.dat'
path_tx2 = '../tx2/'

#to_tx2(path to file, path to output folder, electrode spacing, 
#optional: number of gates (default 20), optional: pulse length (default 2s)))
#example: electrode spacing = 2m

to_tx2_func.to_tx2(path_filt, path_tx2, 2, pulselength=1)
