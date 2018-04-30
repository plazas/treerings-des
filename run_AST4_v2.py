#!/usr/bin/python

# NAME:
#     run_AST4_v2.py
#
# PURPOSE:	
#     Code to plot astrometric residuals vector fields along DECam focal plane (AST4), 
#     according to Gary Bernstein's Validation document (ref ?)
#
# EXPLANATION:
#     For further information see the DES WL Data Acceptance document by
#     Gary Bernstein
#
# CALLING SEQUENCE:
#     ./run_AST4_v2.py --config=<config_file.config> --outfile=<outfile.fits> <ReMatch_cat.out> 
#
#     The order of this calling sequence should be respected to prevent a
#     runtime error.
#
# INPUTS:
#     <config_file>      - This should be the filename of the run_AST4.py
#                          configuration file.  Please see the example AST4.config
#                          for details of its contents and further documentation.
#     
#     <ReMatch_cat.fits>  - This is the input file for this program, in fits format (a binary table).

#                          It's the output catalog of Gary Berstein's ReMatch (a program
#                          that improves on E. Bertin's SCAMP
#                          astrometric solution). 
#     
# OUTPUT:
#     <outfile.pdf>     - Output file with the plots in PDF format. 
#
#
# MODIFICATION HISTORY:
#   - Written by Andres Alejandro Plazas. Modified to work with binary fits tables. OCT2012
#


### From Gary:

"""
----------------------------------------------
Outputs from WCSFit.cpp

Binary FITS tables in output file: (single fits file, with several 'extensions' or tables)

HDU name 'Exposures'
...one row for every exposure used in matching
Name (string)
RA (double, degrees)
Dec (double, degrees)
FieldNumber (int)
InstrumentNumber (int)

HDU name 'Extensions'
...one row for every input catalog that was read.
Filename (string)
HDUNumber (int)
ExposureNumber (int)
DeviceNumber (int)
... more...

HDU name 'Instrument' (could be several)
...one table for each Instrument, each row is one Device of instrument
Name (string)
XMin (double)
XMax (double)
YMin (double)
YMax (double)


HDU name 'WCSOUT':
Sequence (int)
Reserved (bool)
Clipped (bool)
Extension (long)  (which input binary fits extension catalog it comes from)
Object (long)     (row number or other ID in the input catalog)
WtFrac (double)
[xy]pix (double)
[xy]w  (double)
[xy]respix (double)
[xy]resw (double)
sigpix (double)
sigw  (double)

"""


"""
Contents of input file: 

No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU       6   ()           uint8
1    FIELDS      BinTableHDU     15   1R x 3C      ['1D', '1PA(4)', '1D']
2    INSTRUMENT  BinTableHDU     21   62R x 5C     ['1PA(3)', '1D', '1D', '1D', '1D']
3    EXPOSURES   BinTableHDU     19   70R x 5C     [1D, 1J, 1J, 1PA(24), 1D]
4    FILES       BinTableHDU     46   4279R x 18C   [1A, 7A, 7A, 13A, 24A, 1J, 9A, 94A, 5A, 9A, 6A, 33A, 1A, 4A, 18A, 1A, 10A, 10A]
5    EXTENSIONS  BinTableHDU     33   4279R x 12C   ['1J', '1PA(13)', '1J', '1PA(94)', '1J', '1J', '1PA(4)', '1PA(937)', '1PA(18)', '1D', '1PA(10)', '1PA(10)']
6    WCSOUT      BinTableHDU     41   679329R x 16C   ['1L', '1K', '1K', '1L', '1J', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E', '1E']

>>> hdulist['FIELDS'].columns
ColDefs(
    name = 'Dec'; format = '1D'
    name = 'Name'; format = '1PA(4)'
    name = 'RA'; format = '1D'
)
>>> hdulist['INSTRUMENT'].columns
ColDefs(
    name = 'Name'; format = '1PA(3)'
    name = 'XMax'; format = '1D'
    name = 'XMin'; format = '1D'
    name = 'YMax'; format = '1D'
    name = 'YMin'; format = '1D'
)
>>> hdulist['EXPOSURES'].columns
ColDefs(
    name = 'Dec'; format = '1D'
    name = 'FieldNumber'; format = '1J'
    name = 'InstrumentNumber'; format = '1J'
    name = 'Name'; format = '1PA(24)'
    name = 'RA'; format = '1D'
)
>>> hdulist['FILES'].columns
ColDefs(
    name = 'AFFINITY'; format = '1A'
    name = 'DEC'; format = '7A'
    name = 'DEVICE'; format = '7A'
    name = 'ERRORKEY'; format = '13A'
    name = 'EXPOSURE'; format = '24A'
    name = 'EXTENSION'; format = '1J'
    name = 'FIELD'; format = '9A'
    name = 'FILENAME'; format = '94A'
    name = 'IDKEY'; format = '5A'
    name = 'INSTRUMENT'; format = '9A'
    name = 'RA'; format = '6A'
    name = 'SELECT'; format = '33A'
    name = 'STAR_SELECT'; format = '1A'
    name = 'WCSFILE'; format = '4A'
    name = 'WCSOUT'; format = '18A'
    name = 'WEIGHT'; format = '1A'
    name = 'XKEY'; format = '10A'
    name = 'YKEY'; format = '10A'
)
>>> hdulist['EXTENSIONS'].columns
ColDefs(
    name = 'DeviceNumber'; format = '1J'
    name = 'errKey'; format = '1PA(13)'
    name = 'ExposureNumber'; format = '1J'
    name = 'Filename'; format = '1PA(94)'
    name = 'FileNumber'; format = '1J'
    name = 'HDUNumber'; format = '1J'
    name = 'idKey'; format = '1PA(4)'
    name = 'WCS'; format = '1PA(937)'
    name = 'WCSOut'; format = '1PA(18)'
    name = 'Weight'; format = '1D'
    name = 'xKey'; format = '1PA(10)'
    name = 'yKey'; format = '1PA(10)'
)
>>> hdulist['WCSOUT'].columns
ColDefs(
    name = 'Clip'; format = '1L'
    name = 'Extension'; format = '1K'
    name = 'Object'; format = '1K'
    name = 'Reserve'; format = '1L'
    name = 'SequenceNumber'; format = '1J'
    name = 'sigPix'; format = '1E'
    name = 'sigW'; format = '1E'
    name = 'wtFrac'; format = '1E'
    name = 'xPix'; format = '1E'
    name = 'xresPix'; format = '1E'
    name = 'xresW'; format = '1E'
    name = 'xW'; format = '1E'
    name = 'yPix'; format = '1E'
    name = 'yresPix'; format = '1E'
    name = 'yresW'; format = '1E'
    name = 'yW'; format = '1E'

"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Pdf')

import subprocess as S
from sys import argv
from os import system
from scipy import *                # It provides numpy already.
from pylab import *
import matplotlib.pyplot as plt
import pyfits 
from scipy import fftpack          # For periodogram
from AST_h import *                # "Header" with class and function definitons for ASTXX tests. 
from read_config import *          # From Peter Melchior, to read config file and do other stuff. 

from scipy.optimize import leastsq

from matplotlib.backends.backend_pdf import PdfPages


def fit_power_law (pinit, xdata, ydata):
    fitfunc = lambda p,x: p[0]*(x**p[1]) 
    errfunc = lambda p,x,y: y - fitfunc(p,x)

    xdata=array(xdata) # Turn the input list into numpy arrays (allows entry-by-enty calculations). 
    ydata=array(ydata)
    pout =  leastsq (errfunc , pinit , args=(xdata , ydata) , maxfev=20000 )   #routine from scipy.optimize

    xfit = linspace(min(xdata), max(xdata) , 100)
    yfit = fitfunc (pout[0] , xfit)

    return pout, xfit, yfit 


# Doc string ? 
"""
run_AST4.py
    Take output of WCSFits (formerly known as ReMatch.cpp) and calculate AST4 test.  

Usage :  ./run_AST4_v2.py --config=<config_file.config> --outfile=<outfile.fits> <WCSFit.out> 

"""

if len(argv) < 2:
    print " "
    print "run_AST4.py" 
    print "Take *output catalog* of ReMatch.cpp and calculate AST4 test. " 
    print "Usage :  ./run_AST4.py --config=<config_file.config> --outfile=<outfile.pdf> <ReMatch_cat.out>"
    print " "
    exit(1)

#-----------------------------------------------------------------------
# Start of program proper
#-----------------------------------------------------------------------


#P. Melchior's  parser and config file manipulation
(prog , config_file, format, outputfile) = parseCommandLine(argv)  # prog = ifname in this case
#print "prog" , prog
#print type(prog)
if outputfile == None:
    print "./run_AST4_v2.py: Specify an output file name for the results of the test!"
    print "For the moment, this is a pdf with all the plots"
    exit(1)
params_string = getConfigArguments(config_file, format)
print "# DES data acceptance testing:"
print "# running " + prog + " " +params_string
#system(prog + " " + params_string)


# get from parameter file
params = getConfig(config_file)

sizex_reg        =   int (params["SIZE_X_REGION"])   
sizey_reg        =   int (params["SIZE_Y_REGION"])     
sizex_reg_small  =   int (params ["BIN_X"])    
sizey_reg_small  =   int (params ["BIN_Y"]) 
bins_from_edge   =   int (params["BINS_FROM_EDGE"]) 
instrument_name  =   params["INSTRUMENT"]
nccds            =   int (params ["NUMBER_CCDS"] )
SIZE_CCD_X       =   int (params["SIZE_CCD_X"] )
SIZE_CCD_Y       =   int (params ["SIZE_CCD_Y"] )
residuals        =        params["RESIDUALS_UNITS"]
requirement1     =   float (params["REQUIREMENT1"])
filter_to_analize = params["FILTER"]
plate_scale      =   float (params["PLATE_SCALE"])    #plate_scale = 0.238 
use_reserved     =   int (params["USE_RESERVED"])  
wtfrac_threshold =   float (params["WTFRAC_THRESHOLD"])  
if abs(wtfrac_threshold) > 1.:
    print "Parameter %s has invalid value %g" %("WTFRAC_THRESHOLD" , wtfrac_threshold)
    exit(1)
pixels_from_edge =   int (params["PIXELS_FROM_EDGE"])
#min_col          =   int ( params["NUMBER_COLUMNS_IFILE"])
#skip_input_lines =   int (params["SKIP_INPUT_LINES"])
arrow_scale      =   float (params["ARROW_SCALE"])
plot_all_ccds    =   int(params["PLOT_ALL_CCDS"])
fold_data        =   int(params["FOLD_DATA"])
T                =   float(params["FOLDING_PERIOD"])

requirement1_tuple= (requirement1, requirement1/(plate_scale*1000))    

if ( residuals == 'p' ) or (residuals == 'pix'):
    string_residuals_units = 'pixels'
    requirement1=requirement1/(plate_scale*1000)
if (residuals == 'w') or (residuals == 'wcs'):
    string_residuals_units = 'mas'

    
     
#requirement1_tuple= (requirement1, requirement1/(plate_scale*1000))    
print "requirement1_tuple", requirement1_tuple
    
# Create lists of CCDs 
# One "focal plane" list (*not*  a focal_plane object , see AST_h.py) for each type of CCD with different grids, depending on what you want to plot ("ASTa, ASTb , ASTc"). The "focal_plane" class was designed after writing this test, and is used in other AST tests (e.g., AST7)


focal_plane_grid = [CCD(i, sizex_reg , sizey_reg , plate_scale = plate_scale)  for i in range(nccds)]  # CCDs go from 1 to 62, BUT in ReMatch, from 0 to 61.
focal_plane_ver  = [CCD(i, sizex_reg_small , SIZE_CCD_Y , plate_scale = plate_scale)  for i in range(nccds)]  # BCS CCDs 1-8
focal_plane_hor  = [CCD(i, SIZE_CCD_X , sizey_reg_small , plate_scale = plate_scale )  for i in range(nccds)]










#### Read input file

hdulist=pyfits.open(prog)




# Read file
#infile = open(prog, 'r')

#hdulist_input = pyfits.open (prog)  # prog is ifname (#AP)
#table_wcsout = hdulist_input["WCSOUT"].data
#table_exposures = hdulist_input["EXPOSURES"].data
#table_extensions = hdulist_input["EXTENSIONS"].data
#table_instrument = hdulist_input["INSTRUMENT"].data
#hdulis_input.close()


#lines = infile.readlines()   # Don't read everything at once ! 
#infile.close()



mean_title="Mean over all CCD vector fields (%g): %s data \n Input file: %s , Filter: %s " %(nccds, instrument_name, prog.split('/')[-1], filter_to_analize)


#initialize flags of rejected objects
skipped_ref=0
skipped_res=0
skipped_clip=0
skipped_wtfrac=0
skipped_sigw=0

#initialize lists 
x_pix_all = []
y_pix_all = []
x_res_all = []
y_res_all = []

counter_lines=0

#subset of exposures to be plotted
first_exp=1
last_exp=10

ext_vec=[]


min_col=16
#for line in lines[skip_input_lines:]:
for i in range(len(hdulist['WCSOUT'].data)/50):
    #first lines of output catalog of ReMatch don't begin with '#'
    #if (counter_lines < skip_input_lines):
    #    counter_lines+=1
    #    continue
    strline = hdulist['WCSOUT'].data[i]
    #print strline
    #if (strline[0] == '#') or (strline[0][0] == '#'):
    #    continue
    if len(strline) < min_col :
                print "Line does not contain at least %g columns: " %(min_col)
                print line
                exit(1)
    # Put input in local variables


    clip = strline[0]
    extension = strline[1]
    ext = hdulist['EXTENSIONS'].data[extension]['DeviceNumber']
    if ext < 0:   # Do objects with DeviceNumber == -1 belong to ref. catalog?? 
        continue

    # First, get exposure number of object, then get Instrument (i.e., filter). Only get instrument for desired filter (hard coded, for the moment)
    exposure_number = hdulist['EXTENSIONS'].data[extension]['ExposureNumber']
    filter_number   = hdulist['EXPOSURES'].data[exposure_number]['InstrumentNumber']


    if not filter_to_analize == '-1':
        if not filter_number == filter_to_analize:
            continue
    
    obj = strline[2] 
    res = strline[3]
    seq_number = strline[4]  ## ID???
    sigpix = strline[5]
    sigw = strline[6]
    wtfrac = strline[7]
    xpix = strline[8]
    xres_pix = strline[9]
    xres_wcs = strline[10]
    xwcs = strline[11]
    ypix = strline[12]
    yres_pix = strline[13]
    yres_wcs = strline[14]
    ywcs = strline[15]
    

    #Don't use reference catalog objects
    #if (exp < 0):
    #   skipped_ref+=1
    #   continue
    # Don't use reserved objects, unless otherwise especified.
    #if (use_reserved==False):

    #if (res == 1):
    #    skipped_res+=1
    #    continue

    # Don't use clipped objects (how to keep track of them more accurately???)
    #if (clip):
        #skipped_clip+=1
        #continue
    # Don't use objects with wtfrac too high: 
    if (wtfrac > wtfrac_threshold):
        skipped_wtfrac+=1
        continue
    # Don't use objects with large statistical uncertainty
    if sigw > 15: 
        skipped_sigw+=1
        continue

    if (residuals == 'p') or (residuals == 'pix'):
        res_x , res_y = xres_pix , yres_pix
    elif (residuals == 'w') or (residuals == 'wcs'):  # wcs residuals must be flipped for DC5(?): map---->     xpix = ywcs  ; ypix = -xwcs
        res_x , res_y = yres_wcs , -xres_wcs
    else:
        print "residuals parameter format not recognized ('pix','p', 'wcs', 'w')"
        sys.exit(1)

    # Correct residuals for 1/(1-wtfrac) factor. 
    res_x , res_y = res_x/(1-wtfrac) , res_y/(1-wtfrac)  # both with same wtfrac? 

    
    if (xpix > pixels_from_edge) and (xpix < (SIZE_CCD_X - pixels_from_edge) ) and (ypix > pixels_from_edge) and (ypix < (SIZE_CCD_Y - pixels_from_edge) ):
        x_pix_all.append(xpix)
        y_pix_all.append(ypix)
        x_res_all.append(res_x)
        y_res_all.append(res_y)
   
    
    # Identify region in CCDs (focal plane) to where object belongs. 


     #Find out the region in the CCD to where the object belongs.
     # Discard clipped objects 

    i = int ( floor(xpix/ sizex_reg) )
    j = int ( floor(ypix / sizey_reg) )
    #print "i, j in grid ", i, j
    #print "ext " , ext
    if (ext <= len(focal_plane_grid)):   #If number of CCDs < 62
        if (clip):
             focal_plane_grid[ext].grid[i][j].increase_clipped()
             skipped_clip+=1
             continue
        else:
                
             focal_plane_grid[ext-1].grid[i][j].add_residuals( res_x , res_y )  # ext-1 , if data from BCS (ReMatch has ext from 1-8)


    k = int ( floor(xpix/ SIZE_CCD_X) )  
    l = int ( floor(ypix / sizey_reg_small) )
    if (ext <= len(focal_plane_hor)):
        if (clip):
             skipped_clip+=1

             focal_plane_hor[ext].grid[k][l].increase_clipped()
             continue
        else:
              #focal_plane_hor[ext].grid[k][l].add_residuals( deltaxw , deltayw )
              focal_plane_hor[ext-1].grid[k][l].add_residuals( res_x , res_y )


    m = int ( floor(xpix/ sizex_reg_small) )
    n = int ( floor(ypix / SIZE_CCD_Y) )
    if (ext <= len(focal_plane_ver)):
        if (clip):
             skipped_clip+=1
             focal_plane_ver[ext].grid[m][n].increase_clipped()
             continue
        else:
              #focal_plane_ver[ext].grid[m][n].add_residuals( deltaxw, deltayw )
              focal_plane_ver[ext-1].grid[m][n].add_residuals( res_x, res_y )


#Close the file !
hdulist.close()


#print "len(ext_vec), max(ext_vec), min(ext_vec)", len(ext_vec), max(ext_vec), min(ext_vec)


#x_pix_all = array(x_pix_all)
#y_pix_all = array(y_pix_all)
#x_res_all = aropen ray(x_res_all)
#y_res_all = array(y_res_all)

print "Skipped: ref = %g , res = %g , clip = %g , wtfrac = %g, skipped_sigw = %g " %(skipped_ref , skipped_res , skipped_clip , skipped_wtfrac, skipped_sigw)

# Output--->three files ith this format:  x_center | y_center  | mean of x residual | mean 0f y res | instrument | exp | CCD | AST4_flag 

ofile= open ('hola1.dat', 'w')
ofile2 = open ('hola2.dat' , 'w')
ofile3 = open ('hola3.dat' , 'w')




pp = PdfPages(outputfile)



# Create a function that plots all 62 CCDS (ASTa and ASTb) and returns an average of all the CCDs as a CCD object.  
# calculate_mean_CCD : it's in AST_h.py file. 

   
#if len(mean_vec_field_x)!= len (mean_vec_field_y):
#    print " x, y vectors for mean vector field don't have same size (equal to # of CCDS each one)"
#    sys.exit(1)



# Function that calculates 1D vectors (x, y , value_at_x , value_at_y) in a single CCD 
# "vectors_to_plot" (in AST_h.py file).    

#print mean_nobjectsvec

#









####### Bueno, ahora a hacer los benditos "plots"


#ASTa: Vector field


#focal_plane_grid is a list of CCDS. maybe I can do a loop over this, to calculate mean of each extension (getting 62 "mean" CCDs. The mean would be over all exposures). 

#focal_plane_grid=[focal_plane_grid[4]]


mean_CCD_grid = calculate_mean_CCD (focal_plane_grid, pp, plot_all_ccds=plot_all_ccds, nccds=nccds, arrow_scale=arrow_scale, string_residuals_units= string_residuals_units, requirement1_tuple=requirement1_tuple, instrument_name=instrument_name)
#print "mean_CCD_grid: " , mean_CCD_grid
x_grid , y_grid , x_value_grid , y_value_grid , nclipped_vec= vectors_to_plot (mean_CCD_grid, write_ofile=False, ofname='mean_vector_field_AST4.dat')
print "x_grid, y_grid , x_value_grid , y_value_grid: " ,  len(x_grid)  , len(y_grid) , len(x_value_grid), len(y_value_grid)




w, h = figaspect(1.)
fig1=plt.figure(figsize=(w,h))
QP = quiver (x_grid, y_grid, x_value_grid, y_value_grid, color='g', scale_units='inches', scale=arrow_scale)
QK = quiverkey (QP, 0.5, 1.10 , requirement1 ,'%.2g  %s (%.3g mas)' %(requirement1, string_residuals_units, requirement1*mean_CCD_grid.plate_scale*1000.), labelpos='E')
# axhline(y=2000, color='r')
# axvline(x=1000, color='r')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#plt.axes().set_aspect('0.5' , 'datalim')
#title('Astrometric Residuals: Mean of all CCD vector fields.' ) #\n Input file: %s' %(prog))
title(mean_title, size=11)
text(-0.01*mean_CCD_grid.SIZE_CCD_X, 1.1*mean_CCD_grid.SIZE_CCD_Y, 'median: %.3g pix (%.3g mas)' %(mean_CCD_grid.get_median_residual(), mean_CCD_grid.get_median_residual()*mean_CCD_grid.plate_scale*1000.), size=8)
plt.xlabel("Pixels")
plt.ylabel("Pixels")
pp.savefig(fig1)      


#pp.close()

#sys.exit(1)




################


mean_CCD_grid_clipped = plot_clipped_focal_plane (focal_plane_grid, pp, plot_all_ccds=plot_all_ccds, nccds=nccds, arrow_scale=arrow_scale, string_residuals_units= string_residuals_units, requirement1_tuple=requirement1_tuple, instrument_name=instrument_name)
x_grid_clipped , y_grid_clipped , x_value_grid_clipped, y_value_grid_clipped , nclipped_vec = vectors_to_plot (mean_CCD_grid, write_ofile=False, ofname='mean_vector_field_AST4.dat')



scatter (x_grid_clipped, y_grid_clipped, nclipped_vec)
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
title("Clipped objects density", size=11)
#text(-0.01*mean_CCD_grid.SIZE_CCD_X, 1.1*mean_CCD_grid.SIZE_CCD_Y, 'median: %.3g pix (%.3g mas)' %(mean_CCD_grid.get_median_residual(), mean_CCD_grid.get_median_residual()*mean_CCD_grid.plate_scale*1000.), size=8)
plt.xlabel("Pixels")
plt.ylabel("Pixels")
pp.savefig(fig1)      


















## AST4b : 1D plots


#### Vertical binning (in the "x" axis)

mean_CCD_ver = calculate_mean_CCD(focal_plane_ver , pp, plot_all_ccds=plot_all_ccds, nccds=nccds,arrow_scale=arrow_scale, string_residuals_units= string_residuals_units, requirement1_tuple=requirement1_tuple, instrument_name=instrument_name)
x_ver, y_ver , x_ver_val , y_ver_val , nclipped_vec = vectors_to_plot (mean_CCD_ver, write_ofile=True, ofname='vectors_to_plot_AST4.dat')
print "len(x_ver) , min(x_ver),  max(x_ver)" , len(x_ver) , min(x_ver),  max(x_ver)
print "x_ver[0] , x_ver[-1]" , x_ver[0] , x_ver[-1]

figure()
plt.plot(x_ver, x_ver_val)
title('Mean binning in x (average over all CCDs)')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
text(l + 0.7*dx, t , 'Bin size: %g %s' %(sizex_reg_small , string_residuals_units))
ylabel('Residuals (%s)' %(string_residuals_units))
xlabel('Pixels')
pp.savefig()


# Lets fit the part that's close to the edge: GLOWING EDGE effect



#Zoom left

pinit_ver = [1.0, -1.0]
params_ver_left , x_ver_fit_left , x_ver_val_fit_left = fit_power_law (pinit_ver, x_ver[0:bins_from_edge] , x_ver_val [0:bins_from_edge] )

figure()
plt.plot(x_ver[0:bins_from_edge], x_ver_val[0:bins_from_edge] ,'ro')
plt.plot(x_ver_fit_left, x_ver_val_fit_left , 'g') #pout[0] is the parameter array (????)
#l, r, b, t  = axis()
#axis([0, bins_from_edge*sizex_reg_small, -0.20, 0.60])
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
legend( ('Residuals' , ' Exp Fit \n Amplitude: %g \n Exponent: %g' %(params_ver_left[0][0],params_ver_left[0][1]) ) , loc = 'upper right' )
ylabel('Residuals (%s)' %(string_residuals_units))
xlabel('Pixels')
title('Mean binning in x(average over all CCDs): zoom of left side ')
pp.savefig()


#Zoom right


pinit_ver_right = [-1.0, -1.0]
params_ver_right , x_ver_fit_right , x_ver_val_fit_right = fit_power_law (pinit_ver_right, x_ver[len(x_ver)-bins_from_edge:len(x_ver)+1] , x_ver_val [len(x_ver_val)-bins_from_edge:len(x_ver_val)+1] )

figure()
plt.plot(x_ver, x_ver_val ,'ro')
plt.plot(x_ver_fit_right, x_ver_val_fit_right , 'g') #pout[0] is the parameter array (????)
l, r, b, t  = axis()
axis([ (len(x_ver) - bins_from_edge)*sizex_reg_small, SIZE_CCD_X  , 0.20, x_ver_val[-1] ])  # In pixels
legend( ('Residuals' , ' Exp Fit \n Amplitude: %g \n Exponent: %g' %(params_ver_right[0][0],params_ver_right[0][1]) ) , loc = 'upper right',  fancybox=True, numpoints=1 )
ylabel('Residuals (pixels)')
xlabel('Pixels')
title('Mean binning in x (average over all CCDs): zoom of right edge ')
pp.savefig()


########  Horizontal binning (in the "y" axis)


mean_CCD_hor = calculate_mean_CCD (focal_plane_hor, pp, plot_all_ccds=plot_all_ccds, nccds=nccds, arrow_scale=arrow_scale, string_residuals_units= string_residuals_units, requirement1_tuple=requirement1_tuple, instrument_name=instrument_name)
x_hor, y_hor , x_hor_val , y_hor_val , nclipped_vec_hor = vectors_to_plot (mean_CCD_hor)
print "len(y_hor) , min(y_hor),  max(y_hor)" , len(y_hor) , min(y_hor),  max(y_hor)
print "y_hor[0] , y_hor[-1]" , y_hor[0] , y_hor[-1]


figure()
plt.plot(y_hor, y_hor_val)
title('Mean binning in y (average over all CCDs)')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
text(l + 0.7*dx, t  , 'Bin size: %g %s' %(sizey_reg_small, string_residuals_units))
ylabel('Residuals (%s)' %(string_residuals_units))
xlabel('Pixels')
pp.savefig()


### Fit part close to lower edge

pinit_hor_low =[1.0, -1.0]
params_hor_low , y_hor_fit_low , y_hor_val_fit_low = fit_power_law (pinit_hor_low, y_hor[0:bins_from_edge] , y_hor_val [0:bins_from_edge] )


figure()
plt.plot(y_hor[0:bins_from_edge], y_hor_val[0:bins_from_edge], 'ro')
plt.plot(y_hor_fit_low, y_hor_val_fit_low , 'g')
#l, r, b, t  = axis()
#axis([0, bins_from_edge*sizey_reg_small , -0.20, 0.60])
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
legend( ('Mean astrometric residuals' , ' Exp Fit \n Amplitude: %g \n Exponent: %g' %(params_hor_low[0][0],params_hor_low[0][1]) ) , loc = 'upper right',  fancybox=True, numpoints=1 )
ylabel('Astrometric Residuals (%s)' %(string_residuals_units))
xlabel('Pixels')
title('Mean residuals close to upper edge')
pp.savefig()


### Fit part close to upper edge 





## AST4c: Periodogram - Take means every nth column (4 pixel bins, for instance) and then FFT it . Use same vectors as in AST4b


#x

print "Length of mean vector before fft: ",  len(x_ver_val[bins_from_edge:-bins_from_edge])
x_ver_val = array(x_ver_val)
hann_x = hanning(len(x_ver_val))
x_ver_val = x_ver_val*hann_x
X = fft (x_ver_val[bins_from_edge:-bins_from_edge], 2*len(x_ver_val[bins_from_edge:-bins_from_edge]))
nx = len(X)
print nx
powerx = abs(X[1:nx/2])**2
nyqx = 1./(2*sizex_reg_small)   # Nyquist freq 
freqx = array (range(nx/2)) / (nx/2.)*nyqx
periodx = 1./freqx

figure()
plt.plot(freqx[1:len(freqx)] , powerx)
ylabel('Power')
xlabel('1/Pixels')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#annotate('Sawtooth function in x \n with T=14 pixels', xy=(1/14., dy/2.), xytext=(0.75*dx, dy/2.),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
title('Power Spectrum in x' )
pp.savefig()

#0.082 , 0.0815 ----> peak in x

"""
#zoom
figure()
plt.plot(freqx[1:len(freqx)] , powerx)
ylabel('Power')
xlabel('1/Pixels')
axis([0.080, 0.083, 0.00, 0.001])
#l, r, b, t = axis()
#dx , dy = r-l , t-b
#axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#annotate('Sawtooth function in x \n with T=14 pixels', xy=(1/14., 0.08), xytext=(1/50., 0.03),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
title('Power Spectrum in x - Zoom in' )
pp.savefig()
"""
"""
figure()
plt.plot(periodx[1:len(periodx)] , powerx)
#axis([0, 50, -0.0002, 0.0030])
#axis([0, 50, -0.001, 0.004])
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#annotate('Sawtooth function in x \n with T=14 pixels', xy=(14., dy/2.), xytext=(25, dy/2.),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
ylabel('Power')
xlabel('Pixels')
title('Periodogram (x) ')
pp.savefig()
"""

"""
#zoom
figure()
plt.plot(periodx[1:len(periodx)] , powerx)
axis([0, 50, -0.002, 0.040])
#axis([0, 50, -0.01, 0.15])
#l, r, b, t = axis()
#dx , dy = r-l , t-b
#axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
annotate('Sawtooth function in x \n with T=14 pixels', xy=(14., 0.02), xytext=(25., 0.030),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
ylabel('Power')
xlabel('Pixels')
title('Periodogram (x) - Zoom in  ')
pp.savefig()
"""

print "FOLD", fold_data
########## FOLDING 
if (fold_data == True):

    x_ver= array(x_ver)

    #x_ver_val_fold = (4*x_ver_val) - floor ((4*x_ver_val)/14.)

    #print min(x_ver_val_fold) , max(x_ver_val_fold)

### x 
# chop glowing edges data
    #x_pix_all = x_pix_all[20:-20]
    #x_res_all = x_res_all[20:-20]


    for f in arange(0.0816, 0.0823, 0.00005):
        T=1./f
        x_pix_fold = (x_pix_all*f) - floor (x_pix_all*f)
        x_pix_fold_bins = [ [] for i in arange(1/0.05)]
        print "min(x_pix_fold) , max(x_pix_fold)" , min(x_pix_fold) , max(x_pix_fold) 

        print "x_pix_fold[0:10]" , x_pix_fold[0:10] 
        phase_bins = array ( [i for i in arange(0.0, 1.0+0.05, 0.05)])
        print "phase_bins" , phase_bins
        inds_x = digitize( x_pix_fold , phase_bins)
        print "x_pix_fold " , len(x_pix_fold) , x_pix_fold[0:10] , x_pix_fold[-11:-1]
        print "inds_x" , len(inds_x) , inds_x[0:10] , inds_x[-11:-1]

        for i in range(x_pix_fold.size):
            x_pix_fold_bins[inds_x[i]-1].append(x_res_all[i])
        print "len x_pix_fold_bins " , len(x_pix_fold_bins)

        mean_x_pix_fold = [ ]
        std_err_mean_x_pix_fold = [ ]

        for i in range(len(x_pix_fold_bins)):
            print "bin number " , i
            print "min(x_pix_fold_bins[i]) " , min(x_pix_fold_bins[i])
            print "max(x_pix_fold_bins[i]) " , max(x_pix_fold_bins[i])
            mean_x_pix_fold.append( mean (x_pix_fold_bins[i]) )
            std_err_mean_x_pix_fold.append(std(x_pix_fold_bins[i])/sqrt(len(x_pix_fold_bins[i])) )

        print "mean_x_pix_fold"  , len (mean_x_pix_fold) ,  mean_x_pix_fold
        print "std_err_mean_x_pix_fold " , len (std_err_mean_x_pix_fold) , std_err_mean_x_pix_fold

        phase_bins = phase_bins[:-1]
        phase_bins+=(0.05/2)
        print "len (phase_bins)" , len(phase_bins)
        figure()
        errorbar(phase_bins, mean_x_pix_fold, yerr=std_err_mean_x_pix_fold, ecolor='r', label='folded x data' , fmt='.')
        xlabel("f=%.6g , T=%g" %(f , T) )
        ylabel("Mean of residuals (%s)" %(string_residuals_units))
        title("Folding-x")
        pp.savefig()


print "Length of mean vector before fft: ",  len(y_hor_val[bins_from_edge:-bins_from_edge])
y_ver_val = array(y_ver_val)
hann_y = hanning(len(y_ver_val))
y_ver_val = y_ver_val*hann_y 
Y = fft (y_hor_val[bins_from_edge:-bins_from_edge] , 2*len(y_hor_val[bins_from_edge:-bins_from_edge]))
ny = len(Y)
powery = abs(Y[1:ny/2])**2
nyqy = 1./(2*sizey_reg_small)   # ??? Is this the right Nyquist freq ?
freqy = array (range(ny/2)) / (ny/2.)*nyqy
periody = 1./freqy

figure()
plt.plot(freqy[1:len(freqy)] , powery)
ylabel('Power')
xlabel('1/Pixels')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#annotate('Sawtooth function in y \n with T=17 pixels', xy=(1./17, dy/2.), xytext=(0.65*dx, dy/2.),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
title('Power Spectrum in y' )
pp.savefig()

"""
#zoom
figure()
plt.plot(freqy[1:len(freqy)] , powery)
ylabel('Power')
xlabel('1/Pixels')
axis([0.04, 0.043 , -0.01, 0.015])
#l, r, b, t = axis()
#dx , dy = r-l , t-b
#axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
#annotate('Sawtooth function in y \n with T=17 pixels', xy=(1./17, 0.05), xytext=(0.09, 0.06),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
title('Power Spectrum in y - Zoom in-' )
pp.savefig()
"""

"""
figure()
plt.plot(periody[1:len(periody)] , powery)
#axis([0, 50 , -0.001, 0.012])
#axis([0, 50 , -0.01, 0.1])
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
ylabel('Power')
xlabel('Pixels')
#annotate('Sawtooth function in y \n with T=17 pixels', xy=(17, dy/2.), xytext=(0.65*dx, dy/2.),
#            arrowprops=dict(facecolor='black', shrink=0.05),
#            )
title('Periodogram (y) ')
pp.savefig()
"""

if (fold_data==True):
    for f in arange(0.040, 0.042, 0.0005):
        T=1./f
        y_pix_fold = (y_pix_all*f) - floor (y_pix_all*f)
        y_pix_fold_bins = [ [] for i in range(10)]
        print "min(y_pix_fold) , max(y_pix_fold)" , min(y_pix_fold) , max(y_pix_fold) 

        print "y_pix_fold[0:10]" , y_pix_fold[0:10] 
        phase_bins = array ( [i for i in arange(0.0, 1.0+0.1, 0.1)])
        print "phase_bins" , phase_bins
        inds_y = digitize( y_pix_fold , phase_bins)
        print "y_pix_fold " , len(y_pix_fold) , y_pix_fold[0:10] , y_pix_fold[-11:-1]
        print "inds_y" , len(inds_y) , inds_y[0:10] , inds_y[-11:-1]

        for i in range(y_pix_fold.size):
            y_pix_fold_bins[inds_y[i]-1].append(y_res_all[i])
        print "len y_pix_fold_bins " , len(y_pix_fold_bins)

        mean_y_pix_fold = [ ]
        std_err_mean_y_pix_fold = [ ]

        for i in range(len(y_pix_fold_bins)):
            print "bin number " , i
            print "min(y_pix_fold_bins[i]) " , min(y_pix_fold_bins[i])
            print "max(y_pix_fold_bins[i]) " , max(y_pix_fold_bins[i])
            mean_y_pix_fold.append( median (y_pix_fold_bins[i]) )
            std_err_mean_y_pix_fold.append(std(y_pix_fold_bins[i])/sqrt(len(y_pix_fold_bins[i])) )

        print "mean_y_pix_fold"  , len (mean_y_pix_fold) ,  mean_y_pix_fold
        print "std_err_mean_y_pix_fold " , len (std_err_mean_y_pix_fold) , std_err_mean_y_pix_fold

        phase_bins = phase_bins[:-1]
        phase_bins+=0.05
        print "len (phase_bins)" , len(phase_bins)
        figure()
        errorbar(phase_bins, mean_y_pix_fold, yerr=std_err_mean_y_pix_fold, ecolor='r', label='folded y data' , fmt='.')
        xlabel("f=%.6g , T=%g" %(f , T) )
        ylabel("Mean of residuals (%s)" %(string_residuals_units))
        title("Folding-y")
        pp.savefig()



    """
    #zoom
    figure()
    plt.plot(periody[1:len(periody)] , powery)
    axis([0, 50 , -0.001, 0.015])
    #axis([0, 50 , -0.01, 0.18])
    #l, r, b, t = axis()
    #dx , dy = r-l , t-b
    #axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
    ylabel('Power')
    xlabel('Pixels')
    annotate('Sawtooth function in y \n with T=17 pixels', xy=(17., 0.008), xytext=(30, 0.010),
                arrowprops=dict(facecolor='black', shrink=0.05),
                )
    title('Periodogram (y) -Zoom in- ')
    pp.savefig()
"""

pp.close()

