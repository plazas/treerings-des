# coding=utf-8
#!/usr/bin/python

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Pdf')

import subprocess as S
from sys import argv
from os import system
from scipy import weave
from scipy.weave import converters
from pylab import *
import matplotlib.pyplot as plt
import pyfits
#from scipy import fftpack          # For periodogram
#from AST_h import *                # "Header" with class and function definitons for ASTXX tests. 
#from read_config import *          # From Peter Melchior, to read config file and do other stuff. 

#from scipy.optimize import leastsq

import matplotlib.font_manager as fm

from matplotlib.backends.backend_pdf import PdfPages

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import LSQUnivariateSpline

from AST_h import *  

from scipy.interpolate import interp1d

from scipy import optimize

import time


## Function to reject outliers from a list (np.array)
def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

#def reject_outliers(data, m = 2.):
    #start = time.time() 
#    d = np.abs(data - np.median(data))
    #end = time.time()
    #print "Time it took to calculate data minus median in outlier rejection: %g secs" %(end - start)
    #start= time.time()
#    mdev = np.median(d)
    #end=time.time() 
    #print "Time it took to calculate data median of d= data - median(data) in outlier rejection: %g secs" %(end - start)
#    s = d/mdev if mdev else 0
    
    #start=time.time()
#    a=data[s<m]
    #end=time.time()
    #print "Time it took to calculate data[s<m] in outlier rejection: %g secs" %(end - start)
#    return a

if not len(argv) == 11:
    print " "
    print "radial_profile_rings.py <flat_image_name.fits> <astrometric_residuals_table.fits> <center_ring_x>(in ccd coordinates) <center_ring_y> (ccd coordinates) <delta_pix> <filter number> (0:i, 1:g, 2:r, 3:z, 4:Y) <CCD_name> (i.e., N22) <glowing_edge_mask>(pixels)  <ofname.pdf> <output_directory_name> "
    print "  "
    print "  "
    sys.exit(1)





hdulist_flat=pyfits.open(sys.argv[1])
hdulist=pyfits.open(sys.argv[2])
center_x= (float(sys.argv[3])) - 1   #Subtract 1  to match 0-indexed coordinates  of array
center_y= (float(sys.argv[4])) - 1
delta_pix= float (sys.argv[5])
filter = int(sys.argv[6])
ccd_name = sys.argv[7]
ccd_number = CCDNametoNumber_dict [ccd_name]

glowing_edge_mask = int(sys.argv[8])

ofilename= sys.argv[9]

odname= sys.argv[10] 

filter_number_to_string = {0:"i", 1:"g", 2:"r", 3:"z", 4:"Y"}

filter_string = filter_number_to_string[filter]



pp= PdfPages(odname + '/' + ofilename)


image=hdulist_flat['PRIMARY'].data
image=np.transpose(image)    ## 0-indexed. The CCD image is 1-indexed. 
print image.shape
#print type(image[50,50])
image = np.array(image, dtype=np.float64)

#sys.exit(1)
dist_to_center, adu_vec = [], []



print "before for loop. "


#def calculate_distance (image, center_x, center_y, ccd_section='all'):
#    dist_to_center_left, adu_vec_left = [], []
#    dist_to_center_right, adu_vec_right = [], []    


#    for (i,j), value in np.ndenumerate(image):
        # No glowing edge pixels
#        if j < glowing_edge_mask or j > (4095 - glowing_edge_mask) or i < glowing_edge_mask or i > (2047 - glowing_edge_mask):
#            continue
        
#        d= math.sqrt( (center_x - i)**2 + (center_y - j)**2 )
#        if i <= 1023:  
#            dist_to_center_left.append(d)
#            adu_vec_left.append(value)
#        else:
#            dist_to_center_right.append(d)
#            adu_vec_right.append(value)
    
#    return np.array(dist_to_center_left), np.array(adu_vec_left), np.array(dist_to_center_right), np.array(adu_vec_right)

#dist_to_center_left, adu_vec_left, dist_to_center_right, adu_vec_right  = calculate_distance (image, center_x, center_y)
#"""




support = """
#include <math.h> 
#include <iostream>
"""

def call_weave(image_matrix, dist_to_center_left, adu_vec_left, dist_to_center_right, adu_vec_right, image_dim_x, image_dim_y, center_x_var, center_y_var, glowing_edge_mask ):
        #print "image[50,50]", image[50,50]
        #sys.exit(1)
        #print type(image)

        code="""
        for (int i=0; i < image_dim_x; ++i) { 
            for (int j=0; j < image_dim_y; ++j) {        
                double value = IMAGE_MATRIX2(i,j);
                
                //std::cerr << IMAGE_MATRIX2(50,50) << std::endl;  
                
                // Exclude glowing edges 
                if ( (j < glowing_edge_mask) || (j > (4095 - glowing_edge_mask)) || (i < glowing_edge_mask) || (i > (2047 - glowing_edge_mask)) ) {  
                    continue; }
                // Exclude those pixels close to the A/B line
                if ( (i >  (1024 - 50))  && (i < (1025 + 50)) ) {continue;}  
                // Exclude tape bumps
                if ( (i < 125 ) && (j < 80) ) {continue; }   // lower left
                if ( (i < 125) && (j > 2020) && (j < 2110)) {continue; }  // middle left
                if ( (i > 40) && (i < 135) && (j > 4020) && (j < 4096)) { continue; }  // upper left
                
                if ( (i > 1920) && (i < 2046) && (j > 4000) && (j < 4096) ) {continue;}   // upper right
                if ( (i > 1930) && (i < 2046) && (j > 2040) && (j < 2140) ) {continue;}   // middle right
                if ( (i > 1900) && (i < 2046) && (j < 110 ) ) {continue; }   // lower right


                //std::cerr << "value: " << value << std::endl ; 
                double d= sqrt( (center_x_var - i)*(center_x_var - i) + (center_y_var - j)*(center_y_var - j) );
                if (i <= 1023) {
                    dist_to_center_left.append(d);
                    adu_vec_left.append(value);
                }
                else {
                    dist_to_center_right.append(d);
                    adu_vec_right.append(value);
                   } 
    
            }
        }
        """
        weave.inline (code, [ 'image_matrix', 'dist_to_center_left', 'adu_vec_left', 'dist_to_center_right', 'adu_vec_right', 'image_dim_x', 'image_dim_y', 'center_x_var', 'center_y_var', 'glowing_edge_mask'],  support_code = support, libraries=['m'])


def calculate_distance (image_matrix, center_x_var, center_y_var):
    dist_to_center_left, adu_vec_left = [], []
    dist_to_center_right, adu_vec_right = [], []

    image_dim_x = image_matrix.shape[0]
    image_dim_y = image_matrix.shape[1]    
    call_weave (image_matrix, dist_to_center_left, adu_vec_left, dist_to_center_right, adu_vec_right, image_dim_x, image_dim_y, center_x_var, center_y_var, glowing_edge_mask)

    return np.array(dist_to_center_left), np.array(adu_vec_left), np.array(dist_to_center_right), np.array(adu_vec_right)


dist_to_center_left, adu_vec_left, dist_to_center_right, adu_vec_right  = calculate_distance (image, center_x, center_y)





### Now bin resuts by distance from center, i.e, by radius

def calculate_binned_adu (dist_to_center, adu_vec, radial_bins):
    bins_index_vec = np.digitize (dist_to_center, radial_bins)
    #print bins_index_vec
    #print len(dist_to_center)
    #print len(adu_vec)
    #print len(bins_index_vec)
    start = time.time()
    #temp=[]
    #len_radial_bins = len(radial_bins)
    #m=3
    #bins_index_vec = bins_index_vec==1
    #code="""
    #for (int i=0; i <  len_radial_bins ; ++i) { 
    #    temp.append(reject_outliers(adu_vec, m=m)); 
    #}
    #"""
    #weave.inline(code, ['temp', 'adu_vec','bins_index_vec','reject_outliers', 'len_radial_bins', 'm'], support_code = support, libraries=['m'])
    temp=[]
    print "adu_vec", adu_vec[0:30]
    adu_vec=np.array(adu_vec, dtype=np.float64)
    #sys.exit(1)
    for i in range(1, len(radial_bins)):
        #start=time.time()
        #index=bins_index_vec==i
        #end = time.time()
        #print "Time it took to get index vec in reject outliers: %g secs" %(end - start)
       
        #start=time.time()
        #vec=[]
        #print "index[0:1000]", index[:-1000] 
        #sys.exit(1)
        #code="""
        #for (int i =0; i < Nindex[0] ; ++i) {
        #    std::cerr << INDEX1(i) << std::endl; 
        #    if ( INDEX1(i) == 1) { 
        #        std::cerr<< ADU_VEC1(i) << std::endl; 
        #        vec.append(ADU_VEC1(i)); 
        #        }
        #} 
        #"""
        #weave.inline(code, ['index', 'adu_vec', 'vec'], support_code = support, libraries=['m'])
        #end = time.time()
        #print "Time it took to evaluate weave in reject outliers: %g secs" %(end - start)        
        
        #start=time.time()
        #vec= adu_vec[index]    
        #end = time.time()
        #print "Time it took to get evaluate index vec in reject outliers: %g secs" %(end - start)

        #start=time.time() 
        temp.append(reject_outliers(adu_vec[bins_index_vec==i], m=3))
        
        #vec=np.array(vec)
        #print "vec[0:100]", vec[0:100]
        #temp.append(reject_outliers(vec, m=3))
        #print "temp[0:100]", temp[0:100]
        #sys.exit(1)
        #end = time.time()
        #print "Time it took to go through one loop of reject outliers: %g secs" %(end - start)
    
    temp=np.array(temp) 
    end = time.time()
    print "Time it took to go through reject outliers: %g secs" %(end - start)
    start=time.time()
    binned_adu_vec = np.array ([vec.mean() for vec in temp ] )
    end = time.time()
    print "Time it took to go through binned_adu_vec: %g secs" %(end - start)
    #binned_adu_vec_err = np.array( [( vec.std() / math.sqrt(len(vec) ) ) for vec in temp] )
    binned_adu_vec_err=[]
    radial_bins_centers=[]
    
    for i in range(len(radial_bins)-1):
        radial_bins_centers.append((radial_bins[i+1] + radial_bins[i])/2.)
    
    return np.array(radial_bins_centers), np.array(binned_adu_vec), np.array(binned_adu_vec_err)


# Create the radial bins. use the same radial bins for the 
max_distance_l = max (dist_to_center_left)
min_distance_l = min (dist_to_center_left)
print "max_distance_left, min_distance_left", max_distance_l, min_distance_l
max_distance_r = max (dist_to_center_right)
min_distance_r = min (dist_to_center_right)
print "max_distance_right, min_distance_right", max_distance_r, min_distance_r

max_distance = max(max_distance_l, max_distance_r)
min_distance = min(min_distance_l, min_distance_r)

print "max_distance, min_distance", max_distance, min_distance

nbins= math.ceil( (max_distance - min_distance)/ delta_pix )
#print "n bins", nbins
radial_bins = linspace (min_distance, max_distance, nbins)

#print "radial_bins", radial_bins
print "len(radial_bins)", len(radial_bins)

#rbins = open ("radial_bins.dat", 'w')
#for i in radial_bins: 
#    rbins.write("%g \n" %i)
#rbins.close()
#sys.exit(1) 

print "adu_vec_left[0:50]", adu_vec_left[0:50]
binned_distance_left, binned_adu_vec_left, binned_adu_vec_err_left = calculate_binned_adu ( dist_to_center_left,  adu_vec_left, radial_bins)
print "binned_distance_left[0:20]", binned_distance_left[0:100]
print "binned_adu_vec_left[0:20]", binned_adu_vec_left[100:300]

binned_distance_right, binned_adu_vec_right, binned_adu_vec_err_right = calculate_binned_adu ( dist_to_center_right,  adu_vec_right, radial_bins)  
print "binned_distance_right[0:20]", binned_distance_right[0:20]
print "binned_adu_vec_right[0:20]", binned_adu_vec_right[100:300]


 

#if not np.array_equal(binned_distance_left, binned_distance_right):
#    print "Unequal arrays"
#    sys.exit(1)


#binned_distance = binned_distance_left


##### Subtract the function we found from original flat  image to see if pattern went away
# Do it only for 'g' filter to save space and time. 
#binned_adu_vec_av = (binned_adu_vec_right + binned_adu_vec_left) / 2.




#def movingaverage(interval, window_size):
#    window= np.ones(int(window_size))/float(window_size)
#    return np.convolve(interval, window, 'same')





##### Get rid of large-scale oscillations in flat. Use splines (penalizing splines (?))

binned_distance_left, binned_distance_right , binned_adu_vec_left, binned_adu_vec_right = np.nan_to_num (binned_distance_left), np.nan_to_num(binned_distance_right),  np.nan_to_num (binned_adu_vec_left ), np.nan_to_num (binned_adu_vec_right)


index_left = binned_adu_vec_left != 0.
index_right = binned_adu_vec_right != 0.

binned_distance_left = binned_distance_left [index_left]
binned_adu_vec_left = binned_adu_vec_left [index_left]

binned_distance_right = binned_distance_right [index_right]
binned_adu_vec_right = binned_adu_vec_right [index_right]

"""
print "Before creating image with no rings "
print image[50,50]


if filter_string == 'g':
    for (i,j), value in np.ndenumerate(image):
        #if ( (j < glowing_edge_mask) or (j > (4095 - glowing_edge_mask)) or (i < glowing_edge_mask) or (i > (2047 - glowing_edge_mask)) ):
            #continue 
        d= math.sqrt( (center_x - i)**2 + (center_y - j)**2 )
        #print "DISTANCE", d
        if (i <= 1023):
            interp_value = np.interp (d, binned_distance_left, binned_adu_vec_left )
        else:
            interp_value = np.interp (d, binned_distance_right, binned_adu_vec_right )
        #print "a", image[i,j]
        #print "b", interp_value
        image[i,j]/= interp_value
        #print "c", image[i,j]

    output_flat_name = sys.argv[1].split("/")[-1].split(".")[0] + "_no_rings.fits"


    fits_odname = odname + '/' + "FITS_NO_RINGS"

    #cmd = "rm -r %s" %(fits_odname)
    #S.Popen ([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

    cmd = "mkdir " + fits_odname
    S.Popen ([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

    #Remove any existing file with this name 
    cmd = "rm " + fits_odname + '/' + "%s" %(output_flat_name)
    S.Popen ([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
    
    hdu=pyfits.PrimaryHDU (np.transpose(image))
    hdu.writeto(fits_odname + '/' + output_flat_name)
    #hdulist_flat.writeto(fits_odname + '/' + output_flat_name)


print "Final"
#sys.exit()
"""



#smooth_left =0.5*len(binned_adu_vec_left)*np.var(binned_adu_vec_left)
#smooth_right =0.5*len(binned_adu_vec_right)*np.var(binned_adu_vec_right)

#smooth_left=0.03
#smooth_right=0.03


first, last = 0, -1


#a_left, a_right = range (100, len(binned_distance_left)-100, 50) , range (100, len(binned_distance_right)-100, 50)

#knots_left, knots_right = [], []
#for i in a_left: 
#    knots_left.append (binned_distance_left[i])



# Just a few knots for each spline
#knots_left = [ binned_distance_left[100], binned_distance_left[-100]]
#knots_right = [ binned_distance_right[100], binned_distance_right[-100]]

dl=len(binned_distance_left)
dr=len(binned_distance_right)

a = arange (0.15, 1., 0.1)



knots_left = [ binned_distance_left[int(dl*a[0])], binned_distance_left[int(dl*a[1])], binned_distance_left[int(dl*a[2])], binned_distance_left[int(dl*a[3])], binned_distance_left[int(dl*a[4])],  binned_distance_left[int(dl*a[5])], binned_distance_left[int(dl*a[6])]  ]   # This will produce about n+2 knots
knots_right = [ binned_distance_right[int(dr*a[0])], binned_distance_right[int(dr*a[1])], binned_distance_right[int(dr*a[2])], binned_distance_right[int(dr*a[3])], binned_distance_right[int(dr*a[4])],  binned_distance_right[int(dr*a[5])], binned_distance_right[int(dr*a[6])] ]



#y_left_spline = UnivariateSpline ( binned_distance_left[first:last] , binned_adu_vec_left[first:last], s=smooth_left, k=3  )   ### s parameter???
#y_right_spline = UnivariateSpline ( binned_distance_right[first:last] , binned_adu_vec_right[first:last], s=smooth_right, k=3 )


y_left_spline = LSQUnivariateSpline ( binned_distance_left[first:last] , binned_adu_vec_left[first:last], knots_left[0:-1],  k=3  )  
y_right_spline = LSQUnivariateSpline ( binned_distance_right[first:last] , binned_adu_vec_right[first:last], knots_right[0:-1],  k=3 )



#print "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"


y_left = y_left_spline (binned_distance_left)
y_right = y_right_spline (binned_distance_right)


print "y_left", y_left
print "y_right", y_right
#sys.exit(1)

assert (len(y_left) == len(binned_adu_vec_left))
assert (len(y_right) == len(binned_adu_vec_right))





####### this is just the wiggle pattern from the flats without any large-scale oscilation

y_left_diff = (binned_adu_vec_left /  y_left) -1 
y_right_diff =  (binned_adu_vec_right / y_right) -1  


print "Y_TOP_DIFF", y_left_diff
print "Y_BOTTOM_DIFF", y_right_diff


print "y_left_spline.get_knots()", y_left_spline.get_knots()
print "y_right_spline.get_knots()", y_right_spline.get_knots()


#file_nudos=open("nudos.dat", 'w')
#line="knots: %g %g \n" %(len(y_left_spline.get_knots()), len(y_right_spline.get_knots()))
#line2="smooth: %g %g \n" %(smooth_left, smooth_right)
#line3="lens: %g %g %g %g \n " %( len(binned_distance_left), len(y_left_diff),  len(binned_distance_right), len(y_right_diff) )
#file_nudos.write(line)
#file_nudos.write(line2)
#file_nudos.write(line3)
#file_nudos.write ("knots left \n")
#file_nudos.close


"""
##### Print the wiggle patterns of teh left and right channels on 2 separate files. 

outfile_left = "rings_pattern_flat" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_A.dat"  
outfile_right = "rings_pattern_flat" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_B.dat"


f=open(odname + '/' + outfile_left, 'w')
f.write ("#1 r_left (pixels) \n ")
f.write ("#2 rings amplitude from dome flat, left (pixels)  \n ")
for r_var_left, left_var in zip (binned_distance_left, y_left_diff):
    line = "%g %g\n " % (r_var_left, left_var)
    f.write(line)
f.close()


f=open(odname + '/' + outfile_right, 'w')
f.write ("#1 r_right (pixels) \n ")
f.write ("#2 rings amplitude from dome flat, right(pixels)  \n ")
for r_var_right, right_var in zip (binned_distance_right, y_right_diff):
    line = "%g %g \n " % (r_var_right, right_var)
    f.write(line)
f.close()
"""


#sys.exit()



##### Plots about pattern from flat


### Flat profile, left channel

marker_size=4.5
fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot ( binned_distance_left ,  binned_adu_vec_left , 'b.',  label="AMP A", markersize=marker_size)
ax.plot ( binned_distance_left , y_left, 'r', label='large-scale oscillation')
#ax.plot (y_left_spline.get_knots(), y_left_spline(y_left_spline.get_knots()), 'go', label='cubic spline knots')

prop = fm.FontProperties(size=8)
ax.legend(loc='upper right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
## limits
mean_y = median(binned_adu_vec_left)
b,t = mean_y - 0.04*mean_y , mean_y + 0.04*mean_y
plt.ylim([b, t])
ax.set_ylabel ("Average value in radial bin (ADU)  ", size=10)
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


ax=fig.add_subplot(212)
ax.plot ( binned_distance_left,  y_left_diff , 'b.', label='AMP A', markersize=marker_size) 
prop = fm.FontProperties(size=8)
plt.grid(True)
ax.set_xlabel ("Radial distance (pix)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU)  ", size=10)
plt.ylim([-0.017,0.017])
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


fig.suptitle('Tree rings, flat profile (amplifier A) \n  CCD: %s (%g), filter: %s, center of rings: (%g, %g) pix.'  %(ccd_name, ccd_number, filter_string, center_x, center_y), fontsize=12)

pp.savefig()


### Flat profile, right channel

marker_size=4.5
fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot ( binned_distance_right,  binned_adu_vec_right , 'b.',  label="AMP B", markersize=marker_size)
ax.plot ( binned_distance_right, y_right, 'r', label='large-scale oscillation')
#ax.plot (y_right_spline.get_knots(), y_right_spline(y_right_spline.get_knots()), 'go', label='cubic spline knots')
prop = fm.FontProperties(size=8)
ax.legend(loc='upper right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
## limits
mean_y = median(binned_adu_vec_right)
b,t = mean_y - 0.04*mean_y , mean_y + 0.04*mean_y
#plt.xlim([-10, max_distance+10])
#plt.xlim([-10, max_distance+10])
plt.ylim([b, t])
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_ylabel ("Average value in radial bin (ADU)", size=10)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


ax=fig.add_subplot(212)
ax.plot ( binned_distance_right,  y_right_diff , 'b.', label='AMP B', markersize=marker_size)
prop = fm.FontProperties(size=8)
plt.grid(True)
ax.set_xlabel ("Radial distance from center (pix)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU)  ", size=10)
plt.ylim([-0.017,0.017])
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)



fig.suptitle('Tree rings, flat profile (amplifier B) \n  CCD: %s (%g), filter: %s, center of rings: (%g, %g) pix.'  %(ccd_name, ccd_number, filter_string, center_x, center_y), fontsize=12)

pp.savefig()

#pp.close()
#sys.exit()







################## PART2:  Integrate the wiggle pattern to create a prediction of the astrometric error from the flat field image.  



# This function calculates the antderivative
def runtotal(function):
    """Accepts a function, returns (x, running total) pairs --
    analogous to (in) definite integral
    """
    # initialize
    rtotal = 0
    output_x, output_y = [],[]
    lastx,lasty = function[0]   # Each entry of 'function' is a tuple
    #print "lastx", lastx
    #print "lasty", lasty
    #sys.exit(1)
    for x,y in function[1:]:
        #print "x, y", x, y
        avgy   = (lasty + y)/2.0        # average f(x)
        #print "avgy", avgy
        avgx   = (lastx + x)/2.0        # average x
        incre  = x - lastx              # domain increment        
        #print "incre", incre
        rtotal = rtotal + (avgy * incre)    # Accumulates sum
        #print "rtotal", rtotal
        output_x.append(avgx)
        output_y.append(rtotal)
        lastx  = x
        lasty  = y
    return np.array(output_x), np.array(output_y)

#Just checking that the function to calculate the antiderivative works

#xvec=arange(-1,1,0.01)
#x2 = [x*x for x in xvec ]
#test_x, test_y = runtotal ( zip(xvec,x2))

#marker_size=5
#fig=plt.figure()
#ax=fig.add_subplot(111)
#ax.plot (xvec, x2, 'g.', markersize=marker_size,  label="original function")
#ax.plot(test_x, test_y, 'r.', markersize=marker_size, label='antiderivative')
#prop = fm.FontProperties(size=7)
#ax.legend(loc='lower right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
#plt.grid(True)
#ax.set_xlabel ("x", size=11)
#ax.set_ylabel ("y ", size=11)
#plt.xlim([-10, max_distance+10])
#plt.ylim([-0.02, 0.02])
#ax.set_xscale("linear")
#ax.set_title ("Checking code to calculate integral of flat field pattern", size=10)
#pp.savefig()


########## Do integrals
### The function we want to find is   f = (1/r) \int r*w(r) dr, where w(r) is the wiggle pattern. 
#### This function should be the astrometric error. 
first= 0
last = -1
####  Replavce NAN's by zeroes before doing integration, otherwise it will fail completely. 

binned_distance_left = np.nan_to_num (binned_distance_left)
binned_distance_right = np.nan_to_num (binned_distance_right)
y_left_diff = np.nan_to_num (y_left_diff)
y_right_diff = np.nan_to_num (y_right_diff)


print "BEFORE MODEL LEFT*******************************************************************************************************************"

integral_radial_distance_left, integral_wiggle_left  =  runtotal (zip(binned_distance_left[first:last] , binned_distance_left[first:last]*y_left_diff[first:last]))

print "BEFORE MODEL RIGHT*******************************************************************************************************************"
integral_radial_distance_right, integral_wiggle_right  =  runtotal (zip(binned_distance_right[first:last] , binned_distance_right[first:last]*y_right_diff[first:last]))


#if not np.array_equal(integral_radial_distance_left,integral_radial_distance_right):
#    print "Unequal arrays"
#    sys.exit(1)
#integral_radial_distance = integral_radial_distance_left


r_model_left, r_model_right, wiggle_left_model, wiggle_right_model = [],[],[], []
#r_model_file = open ("rings_pattern_int.dat", 'w')

for r_int_left, left_int, r_int_right, right_int in zip( integral_radial_distance_left, integral_wiggle_left, integral_radial_distance_right, integral_wiggle_right):
    #print top[0], top[1], bottom[0], bottom[1]
    if r_int_left == 0:
        wl = 0.
    elif r_int_right == 0:
        wr = 0.
    else:
        wl, wr = -left_int/r_int_left , -right_int / r_int_right

    r_model_left.append(r_int_left)
    r_model_right.append(r_int_right)
    wiggle_left_model.append(wl)
    wiggle_right_model.append(wr)

    ### Print a file with the integrated values: this will be the predicted astrometric signal
    #line = "%g %g %g \n" %(r_int, wl, wr)
    #print line
    #r_model_file.write (line)
#r_model_file.close()



### Remove large-scale oscillations from the predicted model of the astrometric residuals (the integrals of the flat).
### Use a smoothing cubic spline

dl=len(r_model_left)
dr=len(r_model_right)


### smaller steps implies more knots in the spline fit

### DICTIONARY of step_a+left, and step_a_right

##### THESE VALUES ASSUME FIRST VALUE OF ARANGE BELOW IS 0.2


step_a_model_dict = {'N1': [0.08,0.07], 'N2':[0.065,0.07], 'N3':[0.08,0.08], 'N4':[0.07,0.075], 'N5':[0.05,0.05], 'N6':[0.07,0.06], 'N7':[0.07,0.07], 'N8':[0.08,0.08], 'N9':[0.08,0.08], 'N10':[0.085,0.085], 'N11':[0.08,0.05], 'N12':[0.05,0.08], 'N13':[0.08, 0.08], 'N14':[0.08, 0.08], 'N15':[0.08, 0.08], 'N16':[0.065, 0.06], 'N17':[0.05, 0.06], 'N18':[0.08, 0.08], 'N19':[0.07, 0.07], 'N20':[0.08, 0.08], 'N21':[0.08, 0.08], 'N22':[0.07, 0.07], 'N23':[0.08, 0.08], 'N24':[0.08, 0.08], 'N25':[0.08, 0.08], 'N26':[0.08, 0.08], 'N27':[0.08, 0.9], 'N28':[0.08, 0.08], 'N29':[0.09, 0.08], 'N31':[0.08, 0.08], 'S1':[0.06,0.06], 'S2':[0.07, 0.07], 'S3':[0.06, 0.06], 'S4':[0.08, 0.07], 'S5':[0.08, 0.07], 'S6':[0.045, 0.05], 'S7':[0.08, 0.07], 'S8':[0.06, 0.06], 'S9':[0.05, 0.05], 'S10':[0.1, 0.1], 'S11':[0.06, 0.06], 'S12':[0.08, 0.08], 'S13':[0.08, 0.08], 'S14':[0.04, 0.04], 'S15':[0.05, 0.05], 'S16':[0.06, 0.06], 'S17':[0.07, 0.07], 'S18':[0.06, 0.06], 'S19':[0.08, 0.08], 'S20':[0.1, 0.85], 'S21':[0.08, 0.08], 'S22':[0.1, 0.1], 'S23':[0.1, 0.1], 'S24':[0.08, 0.08], 'S25':[0.08, 0.08], 'S26':[0.1, 0.08], 'S27':[0.1, 0.1], 'S28':[0.08, 0.08], 'S29':[0.07, 0.08], 'S30':[0.075, 0.1], 'S31':[0.07, 0.1] }


#N22 = 0.05, 0.055

step_a_left = step_a_model_dict [ccd_name][0]
step_a_right = step_a_model_dict [ccd_name][1]

FIRST_ARANGE= 0.2

a_left = arange  (FIRST_ARANGE, 1.,  step_a_left)
a_right = arange (FIRST_ARANGE, 1., step_a_right)


#a_left = arange  (0.15, 1., 0.3)
#a_right = arange (0.15, 1., 0.3)

#print "a", a

knots_left, knots_right = [], []


#for i in [0,1,2,3,4,5,6,7,9]:
for i in range(len(a_left)):
    knots_left.append(r_model_left[int(dl*a_left[i])])

#for i in [0,1,2,3,4,5,6,7,9]:
for i in range(len(a_right)):
    knots_right.append(r_model_right[int(dr*a_right[i])]) 



y_left_model_spline = LSQUnivariateSpline (r_model_left[first:last], wiggle_left_model[first:last], knots_left[0:-1], k=3 )   ### s parameter???
y_right_model_spline = LSQUnivariateSpline (r_model_right[first:last], wiggle_right_model[first:last], knots_right[0:-1], k=3)


#y_left_model_spline = UnivariateSpline (r_model_left[first:last], wiggle_left_model[first:last], k=3, s=2.5 )   ### s parameter???
#y_right_model_spline = UnivariateSpline (r_model_right[first:last], wiggle_right_model[first:last], k=3, s= 1.5)



print "y_right_model_spline.get_knots()", y_right_model_spline.get_knots()
print "y_left_model_spline.get_knots()", y_left_model_spline.get_knots()



y_left_model = y_left_model_spline (r_model_left)
y_right_model = y_right_model_spline (r_model_right)

print "AFTER SPLINES"
print "len(y_left_model)", len(y_left_model)
print "len(y_right_model)", len(y_right_model)
print "len(r_model_left)", len(r_model_left)
print "len(r_model_right)", len(r_model_right)



wiggle_left_model = np.array(wiggle_left_model)
wiggle_right_model = np.array(wiggle_right_model)

#just subtract, don't divide in this case
y_left_model_diff =  (wiggle_left_model - y_left_model)
y_right_model_diff =  (wiggle_right_model - y_right_model)



##### Print predicted astrometric signal



model_directory_name = "astrometric_models"

cmd = "mkdir " + odname + '/' + model_directory_name
S.Popen ([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()


outfile_left_model = "astrometric_model" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_A.dat"  
outfile_right_model = "astrometric_model" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_B.dat"


f=open(odname + '/' + model_directory_name + '/' + outfile_left_model, 'w')
f.write ("#1 r (distance from ring center - pixels) \n ")
f.write ("#2 radial astrometric value predicted from dome flat  (pixels)  \n ")
for r_var_left_model, left_var_model in zip (r_model_left, -y_left_model_diff):
    line_left_model = "%g %g\n " % (r_var_left_model, left_var_model)
    f.write(line_left_model)
f.close()


f=open(odname + '/' + model_directory_name + '/' + outfile_right_model, 'w')
f.write ("#1 r (distance from ring center - pixels) \n ")
f.write ("#2 radial astrometric value predicted from dome flat  (pixels)  \n ")
for r_var_right_model, right_var_model in zip (r_model_right, -y_right_model_diff):
    line_right_model = "%g %g\n " % (r_var_right_model, right_var_model)
    f.write(line_right_model)
f.close()






###### PLOTS INT???

### Plot model

### Residuals model profile, left channel

marker_size=4.5
fig=plt.figure()
plt.subplots_adjust(hspace=0.1)

ax=fig.add_subplot(211)
ax.plot ( r_model_left, -wiggle_left_model , 'r.',  label="predicted model (AMP A)", markersize=marker_size)
ax.plot ( r_model_left, -y_left_model, 'm', label='large-scale oscillations')
ax.plot (y_left_model_spline.get_knots(), y_left_model_spline(y_left_model_spline.get_knots()), 'go', label='cubic spline knots')
prop = fm.FontProperties(size=7)
ax.legend(loc='lower left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU)", size=10)
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
#plt.xlim([-10, max_distance+10])
#plt.ylim([0.95, 1.05])
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=False)


ax=fig.add_subplot(212)
ax.plot ( r_model_left  ,  -y_left_model_diff , 'r.', label='predicted model without large-scale oscillations (AMP A)', markersize=marker_size)
prop = fm.FontProperties(size=7)
plt.grid(True)
ax.set_xlabel ("Radial distance from center (pixels)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU)", size=10)
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.legend(loc='lower left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
#plt.xlim([-10, max_distance+10])
#plt.ylim([-0.01,0.01])
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


fig.suptitle('Residuals model predicted from flat  (AMP A) \n CCD: %s (%g)' %(ccd_name, ccd_number), fontsize=12)

pp.savefig()



# Residuals model profile, right channel

marker_size=4.5
fig=plt.figure()
plt.subplots_adjust(hspace=0.1)

ax=fig.add_subplot(211)
ax.plot ( r_model_right, -wiggle_right_model , 'r.',  label="predicted model (AMP B)", markersize=marker_size)
ax.plot ( r_model_right, -y_right_model, 'm', label='large-scale oscillation')
ax.plot (y_right_model_spline.get_knots(), y_right_model_spline(y_right_model_spline.get_knots()), 'go', label='cubic spline knots')
prop = fm.FontProperties(size=7)
ax.legend(loc='lower left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU) ", size=10)
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
#plt.xlim([-10, max_distance+10])
#plt.ylim([0.95, 1.05])
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=False)


ax=fig.add_subplot(212)
ax.plot ( r_model_right  ,  -y_right_model_diff , 'r.', label='predicted model without large-scale oscillations (AMP B)', markersize=marker_size)
prop = fm.FontProperties(size=7)
ax.legend(loc='lower left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
ax.set_xlabel ("Radial distance from center (pixels)", size=10)
ax.set_ylabel ("Average value in radial bin (ADU)  ", size=10)
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
#plt.xlim([-10, max_distance+10])
#plt.ylim([-0.01,0.01])
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)

fig.suptitle('Residuals model predicted from dome flat  (AMP B) \n CCD: %s (%g)' %(ccd_name, ccd_number), fontsize=12)


pp.savefig()


#pp.close()
#sys.exit()


############################ Part 3. Read astrometric residuals. Compare them wth predicted flat pattern

####  Read astro residuals from output table produced by Gary's refitting code (WCSFit)
####  We need to read the whole list of astrometric residuals, get the radial componets of the errors, and then bin them in radial bins.


instrument = 2  ## 2,3,4,5,6  : has the names of the CCDs.


#hdulist = pyfits.open (astro_res_fname)

clip = hdulist["WCSOUT"].data.field("Clip")
ext =  hdulist["WCSOUT"].data.field("Extension")

row_in_INSTRUMENT = hdulist["EXTENSIONS"].data.field("DeviceNumber")[ext]
CCD_name = hdulist[instrument].data.field('Name') [row_in_INSTRUMENT]   ##### Need to use hdulist[2] instead of hdulist["INSTRUMENT"], because there are 5 intruments. For the catalog decflat3.resid.fits, this is fine because the INSTRUMENTS are the same. 
CCD_name = np.array([''.join(name) for name in CCD_name])


### Get the filter number: ext-->EXTENSION-->ExposureNumber--->EXPOSURE---> InstrumentNumber
exposure_number = hdulist["EXTENSIONS"].data.field('ExposureNumber')[ext]
instrument_number = hdulist["EXPOSURES"].data.field('InstrumentNumber')[exposure_number]   ## 0=i, 1=g , 2=r , 3=z ,4=Y 

#print instrument_number[:100]
#sys.exit(1)

reserve = hdulist["WCSOUT"].data.field('Reserve')
sig_pix =  hdulist["WCSOUT"].data.field('sigPix')
sig_w  =  hdulist["WCSOUT"].data.field('sigW')
wt_frac =  hdulist["WCSOUT"].data.field('wtFrac')

xpix =  hdulist["WCSOUT"].data.field('xPix')
xres_pix =  hdulist["WCSOUT"].data.field('xresPix')
xres_w =  hdulist["WCSOUT"].data.field('xresW')
xw =  hdulist["WCSOUT"].data.field('xW')

ypix =  hdulist["WCSOUT"].data.field('yPix')
yres_pix =  hdulist["WCSOUT"].data.field('yresPix')
yres_w =  hdulist["WCSOUT"].data.field('yresW')
yw =  hdulist["WCSOUT"].data.field('yW')


### Masks of the data 
wt_cut = 0.6
sig_w_cut = 15


#print CCD_name[0]
print [ccd_name] 

#sys.exit(1)

#   // Exclude glowing edges 
#                if ( (j < glowing_edge_mask) || (j > (4095 - glowing_edge_mask)) || (i < glowing_edge_mask) || (i > (2047 - glowing_edge_mask)) ) {  
#                    continue; }
#                // Exclude those pixels close to the A/B line
#                if ( (i >  (1024 - 50))  && (i < (1025 + 50)) ) {continue;}  
#                // Exclude tape bumps
#                if ( (i < 125 ) && (j < 80) ) {continue; }   // lower left
#                if ( (i < 125) && (j > 2020) && (j < 2110)) {continue; }  // middle left
#                if ( (i > 40) && (i < 135) && (j > 4020) && (j < 4096)) { continue; }  // upper left
#                
#                if ( (i > 1920) && (i < 2046) && (j > 4000) && (j < 4096) ) {continue;}   // upper right
#                if ( (i > 1930) && (i < 2046) && (j > 2040) && (j < 2140) ) {continue;}   // middle right
#                if ( (i > 1900) && (i < 2046) && (j < 110 ) ) {continue; }   // lower right







mask1 =  np.bitwise_and ( np.bitwise_and (  np.bitwise_and( np.bitwise_and(clip == False, wt_frac < wt_cut), sig_w < sig_w_cut), CCD_name == [ccd_name] ), instrument_number > -1)
mask =  np.bitwise_and ( np.bitwise_and ( np.bitwise_and( np.bitwise_and ( mask1, xpix > glowing_edge_mask  ), xpix < (2048 - glowing_edge_mask) ),  ypix > glowing_edge_mask ), ypix < (4096 - glowing_edge_mask) )



#mask = np.bitwise_and(  np.bitwise_and(np.bitwise_and(clip == False, wt_frac < wt_cut), sig_w < sig_w_cut), CCD_name == ['N21'])

#mask = mask1

wt_frac= wt_frac[mask]
CCD_name = CCD_name[mask]

xpix = xpix[mask]
xres_pix = xres_pix[mask]
xres_w = xres_w[mask]
xw = xw[mask]

ypix = ypix[mask]
yres_pix =  yres_pix[mask]
yres_w =  yres_w[mask]
yw =  yw[mask]

intrument_number = instrument_number[mask]


xres_pix = xres_pix/ (1 - wt_frac)


##### Calculate the distances from center, and the radial and orthogonal components of the astrometric residuals

#### Have to do it for left and right channels for separate, to compare from predictions above

index_left = xpix < 1023
index_right = xpix > 1022

xpix_left, ypix_left = xpix[index_left], ypix[index_left]
xpix_right, ypix_right = xpix[index_right], ypix[index_right]


xres_pix_left, yres_pix_left = xres_pix[index_left], yres_pix[index_left]
xres_pix_right, yres_pix_right = xres_pix[index_right], yres_pix[index_right]

instrument_number_left = instrument_number [index_left]
instrument_number_right = instrument_number [index_right] 



radial_component_distance_left = np.sqrt ( (xpix_left - center_x)**2 + (ypix_left - center_y)**2)
radial_component_distance_right = np.sqrt ( (xpix_right - center_x)**2 + (ypix_right - center_y)**2)


print "max(xpix_left), max(radial_component_distance_left)" ,  max(xpix_left), max(radial_component_distance_left)
print "max(xpix_right), max(radial_component_distance_right)", max(xpix_right), max(radial_component_distance_right)


rx_left, ry_left = (xpix_left - center_x) ,  (ypix_left-center_y)
rx_right, ry_right = (xpix_right - center_x) ,  (ypix_right-center_y)


radial_component_res_left = (xres_pix_left*rx_left + yres_pix_left*ry_left) / radial_component_distance_left
orthogonal_component_res_left = (xres_pix_left*ry_left - yres_pix_left*rx_left) / radial_component_distance_left

radial_component_res_right = (xres_pix_right*rx_right + yres_pix_right*ry_right) / radial_component_distance_right
orthogonal_component_res_right = (xres_pix_right*ry_right - yres_pix_right*rx_right) / radial_component_distance_right


### Now bin the astrometric residuals by distance


### Left
max_distance_residuals_left = max (radial_component_distance_left )
min_distance_residuals_left = min (radial_component_distance_left )
print "max_distance_residuals_left, min_distance_residuals_left", max_distance_residuals_left, min_distance_residuals_left
nbins_left= math.ceil( (max_distance_residuals_left - min_distance_residuals_left)/ delta_pix )
print "nbins_left for astrometroc residuals", nbins_left
radial_bins_residuals_left = linspace (min_distance_residuals_left, max_distance_residuals_left, nbins_left)


binned_radial_distance_left, binned_radial_residuals_left, binned_radial_residuals_err_left = calculate_binned_adu ( radial_component_distance_left , radial_component_res_left, radial_bins_residuals_left)
binned_orthogonal_distance_left, binned_orthogonal_residuals_left, binned_orthogonal_residuals_err_left = calculate_binned_adu ( radial_component_distance_left , orthogonal_component_res_left, radial_bins_residuals_left)



###Print the binned astrometric residuals in a file
#"_%s" %ccd_name + "_%s" %(filter_string) + ".dat"

out_residuals_file_left = "binned_astro_res" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_A.dat"

binned_res_file_left = open (odname + '/' + out_residuals_file_left, 'w')
for r_res_var, radial_res_var in  zip (binned_radial_distance_left, binned_radial_residuals_left):
    line = "%g %g \n" %(r_res_var, radial_res_var)
    binned_res_file_left.write (line)
binned_res_file_left.close()



### right

max_distance_residuals_right = max (radial_component_distance_right )
min_distance_residuals_right = min (radial_component_distance_right )
print "max_distance_residuals_right, min_distance_residuals_right", max_distance_residuals_right, min_distance_residuals_right
nbins_right= math.ceil( (max_distance_residuals_right - min_distance_residuals_right)/ delta_pix )
print "nbins_right for astrometroc residuals", nbins_right
radial_bins_residuals_right = linspace (min_distance_residuals_right, max_distance_residuals_right, nbins_right)


binned_radial_distance_right, binned_radial_residuals_right, binned_radial_residuals_err_right = calculate_binned_adu ( radial_component_distance_right , radial_component_res_right, radial_bins_residuals_right)
binned_orthogonal_distance_right, binned_orthogonal_residuals_right, binned_orthogonal_residuals_err_right = calculate_binned_adu ( radial_component_distance_right , orthogonal_component_res_right, radial_bins_residuals_right)



###Print the binned astrometric residuals in a file
#"_%s" %ccd_name + "_%s" %(filter_string) + ".dat"

out_residuals_file_right = "binned_astro_res" + "_%s" %ccd_name + "_%s" %(filter_string) + "_amp_B.dat"

binned_res_file_right = open (odname + '/' + out_residuals_file_right, 'w')
for r_res_var_right, radial_res_var_right in  zip (binned_radial_distance_right, binned_radial_residuals_right):
    line_right = "%g %g \n" %(r_res_var_right, radial_res_var_right)
    binned_res_file_right.write (line_right)
binned_res_file_right.close()


print "max(binned_radial_distance_left), max(binned_radial_distance_right)", max(binned_radial_distance_left), max(binned_radial_distance_right)

# Now we need to get rid of the long-scale oscillations for the astrometric residuals


## Plot radial component of residuals and the predicted residuals calculated from the flats. 

"""

marker_size=4.5

fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot ( binned_radial_distance_left, binned_radial_residuals_left, 'g.', label='radial component of astrometric residuals - left channel', markersize=marker_size)
prop = fm.FontProperties(size=7)
ax.legend(loc='lower right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
ax.set_ylabel ("Astrometric residuals (pix)", size=11)
#plt.xlim([1600, 2800])
plt.ylim([-0.15, 0.15])
ax.set_xscale("linear")

#### othogonal component. should be zero
#ax=fig.add_subplot(222)
#ax.plot ( binned_orthogonal_distance_left, binned_orthogonal_residuals_left, 'g.', label='orthogonal component of astrometric residuals - left channel', markersize=marker_size)
#prop = fm.FontProperties(size=7)
#ax.legend(loc='lower right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
#plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
#ax.set_ylabel ("Astrometric residuals (pix)", size=11)
#plt.xlim([1600, 2800])
#plt.ylim([-0.15, 0.15])
#ax.set_xscale("linear")


ax=fig.add_subplot(212)
ax.plot ( binned_radial_distance_right, binned_radial_residuals_right, 'g.', label='radial component of astrometric residuals - right channel', markersize=marker_size)
prop = fm.FontProperties(size=7)
ax.legend(loc='lower right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
ax.set_ylabel ("Astrometric residuals (pix)", size=11)
ax.set_xlabel ("Radial distance (pixels)", size=10)
#plt.xlim([1600, 2800])
plt.ylim([-0.15, 0.15])
ax.set_xscale("linear")



#### othogonal component. should be zero
#ax=fig.add_subplot(224)
#ax.plot ( binned_orthogonal_distance_right, binned_orthogonal_residuals_right, 'g.', label='orthogonal component of astrometric residuals - right channel', markersize=marker_size)
#prop = fm.FontProperties(size=7)
#ax.legend(loc='lower right' , fancybox=True, ncol=2, numpoints=1, prop = prop)
#plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
#ax.set_ylabel ("Astrometric residuals (pix)", size=11)
#plt.xlim([1600, 2800])
#plt.ylim([-0.15, 0.15])
#ax.set_xscale("linear")



fig.suptitle(" Measured astrometric residuals: radial and orthogonal components \n  CCD: %s" %ccd_name, fontsize=12)
pp.savefig()
"""





######  OVERPLOT ASTRO RESIDUALS AND PREDICTION
fig=plt.figure()
ax=fig.add_subplot(211)
ax.plot ( binned_radial_distance_left, binned_radial_residuals_left, 'g.', label='radial component of astrometric residuals - left channel', markersize=marker_size)
ax.plot ( r_model_left,  y_left_model_diff, 'r.', label='prediction from flat - AMP A', markersize=marker_size)
prop = fm.FontProperties(size=7)
ax.legend(loc='upper left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
#ax.set_xlabel ("Radial distance (pixels)", size=10)
ax.set_ylabel ("Astrometric residuals (pix)", size=10)
plt.ylim([-0.17, 0.17])
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


ax=fig.add_subplot(212)
ax.plot ( binned_radial_distance_right, binned_radial_residuals_right, 'g.', label='radial component of astrometric residuals - right channel', markersize=marker_size)
ax.plot (  r_model_right,  y_right_model_diff, 'r.', label='prediction from flat - AMP B', markersize=marker_size)
ax.legend(loc='upper left' , fancybox=True, ncol=2, numpoints=1, prop = prop)
plt.grid(True)
ax.set_xlabel ("Radial distance from center (pixels)", size=10)
ax.set_ylabel ("Astrometric residuals (pix)", size=10)
plt.ylim([-0.17, 0.17])
ax.set_yticklabels(ax.get_yticks(), size=9, visible=True)
ax.set_xscale("linear")
ax.set_xticklabels(ax.get_xticks(), size=9, visible=True)


#ax.set_title ("Tree rings, flat profile (Top/left channel) ", size=10)

fig.suptitle("Tree rings: astrometric residuals and model of the residuals from flat-field images \n CCD: %s (%g), filter: %s " %(ccd_name, ccd_number, filter_string), fontsize=12)

pp.savefig()

#print len(r_model)
#print len(binned_radial_distance)
#print r_model[0:100]
#print binned_radial_distance [0:100]
######## 

pp.close()
sys.exit(1)



##### Part 4:  Use the model to subtract radial component from measured astrometric residuals.  Make vector plot again and see what remains: 


#radial_component_distance = np.sqrt ( (xpix - center_x)**2 + (ypix-center_y)**2)
#rx, ry = (xpix - center_x) ,  (ypix-center_y)

#radial_component_res = (xres_pix*rx + yres_pix*ry) / radial_component_distance
#orthogonal_component_res = (xres_pix*ry - yres_pix*rx) / radial_component_distance


### x_r = (a_r*r_x + a_orth*r_y) / r 
### y_r = (a_r*r_y - a_orth*r_x) / r


#r_model[-1]+= 1


model_left_function = interp1d (r_model_left, y_left_model_diff, bounds_error=False, fill_value=0.)
model_right_function = interp1d( r_model_right,  y_right_model_diff, bounds_error=False, fill_value=0.)

#min_rmodel, max_rmodel = min(r_model), max (r_model)

#print min_rmodel, max_rmodel
#print min(radial_component_distance)
#print max(radial_component_distance)


### Need to scale model for different bands using factors_left, factors_right

factors_left_dict = {1:1, 2:0.991, 0:0.9512, 3:0.5793, 4:0.4260}
factors_right_dict = {1:1, 2:0.983, 0:0.9435, 3:0.5757, 4:0.4279}


factors_left = np.zeros (len(instrument_number_left))
for k, v in factors_left_dict.iteritems(): 
    #print k, v
    factors_left[ instrument_number_left ==k ] = v

factors_right = np.zeros (len(instrument_number_right))
for k, v in factors_right_dict.iteritems(): 
    #print k, v
    #print "hola", instrument_number_right == k
    factors_right[ instrument_number_right == k ] = v


print "min(instrument_number)", min(instrument_number)
print "instrument_number_left[0:20]", instrument_number_left
print "factors_left[0:20]", factors_left



assert(len(factors_left) == len(radial_component_distance_left))




radial_component_res_left = radial_component_res_left -  factors_left*model_left_function (radial_component_distance_left) 
radial_component_res_right = radial_component_res_right - factors_right*model_right_function (radial_component_distance_right)



xres_pix_new_left  = (radial_component_res_left*rx_left + orthogonal_component_res_left*ry_left) / radial_component_distance_left
yres_pix_new_left  = (radial_component_res_left*ry_left - orthogonal_component_res_left*rx_left) / radial_component_distance_left

xres_pix_new_right  = (radial_component_res_right*rx_right + orthogonal_component_res_right*ry_right) / radial_component_distance_right                                            
yres_pix_new_right  = (radial_component_res_right*ry_right - orthogonal_component_res_right*rx_right) / radial_component_distance_right



sizex_reg = 32
sizey_reg = 32
plate_scale = 0.27




chip = CCD (ccd_name, sizex_reg , sizey_reg , plate_scale = plate_scale)   
chip_new = CCD (ccd_name, sizex_reg , sizey_reg , plate_scale = plate_scale)


for ( xpix_left_var, ypix_left_var,  x_res_old_left_var, y_res_old_left_var, x_res_new_left_var, y_res_new_left_var) in zip( xpix_left, ypix_left, xres_pix_left, yres_pix_left, xres_pix_new_left, yres_pix_new_left ):
    
    i = int ( floor(xpix_left_var/ sizex_reg) )
    j = int ( floor(ypix_left_var/ sizey_reg) )
        
    chip.grid[i][j].add_residuals( x_res_old_left_var , y_res_old_left_var ) 
    chip_new.grid[i][j].add_residuals ( x_res_new_left_var, y_res_new_left_var)


for ( xpix_right_var, ypix_right_var,  x_res_old_right_var, y_res_old_right_var, x_res_new_right_var, y_res_new_right_var) in zip( xpix_right, ypix_right, xres_pix_right, yres_pix_right, xres_pix_new_right, yres_pix_new_right ):

    i = int ( floor(xpix_right_var/ sizex_reg) )
    j = int ( floor(ypix_right_var/ sizey_reg) )
        
    chip.grid[i][j].add_residuals( x_res_old_right_var , y_res_old_right_var)  
    chip_new.grid[i][j].add_residuals ( x_res_new_right_var, y_res_new_right_var)



x_old, y_old, x_value_old, y_value_old, clip_old = vectors_to_plot (chip)
x_new, y_new, x_value_new, y_value_new, clip_new = vectors_to_plot (chip_new)



def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)


font_size=9.5
visible_x, visible_y = True, True
tick_size = 8

#plate_scale = 0.27   # arcsec/pix
requirement1 = 15./ (1000*0.27)   #pix


arrow_scale = 0.5
string_residuals_units = "pixels"


fig=plt.figure()
ax = fig.add_subplot(111)
QP = quiver (x_old, y_old, x_value_old, y_value_old,  color='g', scale_units='inches', scale=arrow_scale)
QK = quiverkey (QP, 1.10, 0.5 , requirement1 ,'%.2g  %s (%.3g mas)' %(requirement1, string_residuals_units, requirement1*plate_scale*1000.), labelpos='E')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
plt.xlim ([-20,2070])
plt.ylim ([-20.,4110])
title(" ", size=10.5)
adjustFigAspect(fig,aspect=.5)
ax.set_xticklabels(ax.get_xticks(), size=tick_size, visible=visible_x)
lx=ax.set_xlabel(r"X (pixels)", visible=visible_x)
lx.set_fontsize(font_size)
ax.set_yticklabels(ax.get_yticks(), size=tick_size, visible=visible_y)
ly=ax.set_ylabel(r"Y (pixels)", visible=visible_y)
ly.set_fontsize(font_size)

fig.suptitle ("Astrometric residuals from star flat (%s) \n All filters. CCD: %s (%g). Median: %g , rms: %g "  %(sys.argv[2].split('/')[-1], ccd_name, ccd_number, chip.get_median_residual(), chip.get_rms_residual()))

pp.savefig(fig)


fig=plt.figure()
ax = fig.add_subplot(111)
QP = quiver (x_new, y_new, x_value_new, y_value_new,  color='g', scale_units='inches', scale=arrow_scale)
QK = quiverkey (QP, 1.10, 0.5 , requirement1 ,'%.2g  %s (%.3g mas)' %(requirement1, string_residuals_units, requirement1*plate_scale*1000.), labelpos='E')
l, r, b, t = axis()
dx , dy = r-l , t-b
axis([l-0.05*dx , r + 0.05*dx , b - 0.05*dy , t + 0.05*dy])
plt.xlim ([-20,2070])
plt.ylim ([-20.,4110])
title(" ", size=10.5)
adjustFigAspect(fig,aspect=.5)
ax.set_xticklabels(ax.get_xticks(), size=tick_size, visible=visible_x)
lx=ax.set_xlabel(r"X (pixels)", visible=visible_x)
lx.set_fontsize(font_size)
ax.set_yticklabels(ax.get_yticks(), size=tick_size, visible=visible_y)
ly=ax.set_ylabel(r"Y (pixels)", visible=visible_y)
ly.set_fontsize(font_size)


fig.suptitle ("Astrometric residuals with prediction from flat (radial component) subtracted \n All filters. CCD: %s (%g). Median: %g, rms: %g "  %(ccd_name, ccd_number, chip_new.get_median_residual(), chip_new.get_rms_residual()))

pp.savefig(fig)

pp.close()


sys.exit()

