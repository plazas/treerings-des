#!/usr/bin/python

# Some classes and function definitions used both by the AST4 and AST7 tests. 
# Written by Andres A. Plazas - May 2011 

from scipy import *                # It provides numpy already.
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

# Taken from Jiangang for central coordinates of
# each CCD
S1=["S1",16.908,191.670]
S2=["S2",16.908,127.780]
S3=["S3",16.908,63.890]
S4=["S4",16.908,0.]
S5=["S5",16.908,-63.890]
S6=["S6",16.908,-127.780]
S7=["S7",16.908,-191.670]

S8=["S8",50.724,159.725]
S9=["S9",50.724,95.835]
S10=["S10",50.724,31.945]
S11=["S11",50.724,-31.945]
S12=["S12",50.724,-95.835]
S13=["S13",50.724,-159.725]

S14=["S14",84.540,159.725]
S15=["S15",84.540,95.835]
S16=["S16",84.540,31.945]
S17=["S17",84.540,-31.945]
S18=["S18",84.540,-95.835]
S19=["S19",84.540,-159.725]
   
S20=["S20",118.356,127.780]
S21=["S21",118.356,63.890]
S22=["S22",118.356,0.]
S23=["S23",118.356,-63.890]
S24=["S24",118.356,-127.780]

S25=["S25",152.172,95.835]
S26=["S26",152.172,31.945]
S27=["S27",152.172,-31.945]
S28=["S28",152.172,-95.835]

S29=["S29",185.988,63.890]
S30=["S30",185.988,0.]
S31=["S31",185.988,-63.890]

N1=["N1",-16.908,191.670]
N2=["N2",-16.908,127.780]
N3=["N3",-16.908,63.890]
N4=["N4",-16.908,0.]
N5=["N5",-16.908,-63.890]
N6=["N6",-16.908,-127.780]
N7=["N7",-16.908,-191.670]

N8=["N8",-50.724,159.725]
N9=["N9",-50.724,95.835]
N10=["N10",-50.724,31.945]
N11=["N11",-50.724,-31.945]
N12=["N12",-50.724,-95.835]
N13=["N13",-50.724,-159.725]

N14=["N14",-84.540,159.725]
N15=["N15",-84.540,95.835]
N16=["N16",-84.540,31.945]
N17=["N17",-84.540,-31.945]
N18=["N18",-84.540,-95.835]
N19=["N19",-84.540,-159.725]

N20=["N20",-118.356,127.780]
N21=["N21",-118.356,63.890]
N22=["N22",-118.356,0.]
N23=["N23",-118.356,-63.890]
N24=["N24",-118.356,-127.780]

N25=["N25",-152.172,95.835]
N26=["N26",-152.172,31.945]
N27=["N27",-152.172,-31.945]
N28=["N28",-152.172,-95.835]

N29=["N29",-185.988,63.890]
N30=["N30",-185.988,0.]
N31=["N31",-185.988,-63.890]


#extidx = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,S26,S27,S28,S29,S30,S31,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31]

extidx = [S29,S30,S31,S25,S26,S27,S28,S20,S21,S22,S23,S24,S14,S15,S16,S17,S18,S19,S8,S9,S10,S11,S12,S13,S1,S2,S3,S4,S5,S6,S7,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31]
# defines the size of a chip.  One pixel=15 microns
xsize=2048*15e-6*1000
ysize=4096*15e-6*1000
corners={}
for i,ext in enumerate(extidx):
    xy=[]
    xy.append(ext[1]-xsize/2)
    xy.append(ext[2]-ysize/2)
    corners[i+1]=xy

def LocaltoFocal(ccd_number,x,y):
    xc=corners[ccd_number][0]
    yc=corners[ccd_number][1]
    return x*15e-6*1000+xc,y*15e-6*1000+yc





class Region:
    """ 
    A 'Region' is a rectangular cell in the CCD. It holds a certain number of residuals (from ReMatch).
    It calculates the mean of these residuals along each Cartesian component.

    """
    def __init__(self):
        # Vectors to hold residuals
        #self.deltax=[ ] 
        #self.deltay=[ ]
        self.sum_w=0
        self.nelements = 0
        self.nclipped=0
        self.sum_deltax, self.sum_deltay =0,0
    def add_residuals (self, dx, dy, w=1):
        self.sum_deltax+=dx*w
        self.sum_deltay+=dy*w
        #self.deltax.append(dx)
        #self.deltay.append(dy)
        self.nelements+=1
        self.sum_w+=1
    def increase_clipped (self):
        self.nclipped+=1
    def add_clipped (self,_nclipped):
        self.nclipped+=_nclipped
    def get_nclipped(self):
        return self.nclipped
    def get_mean_residuals (self):
        if self.nelements == 0:
            return 0,0
	else:
	    return self.sum_deltax/self.sum_w , self.sum_deltay/self.sum_w
            
    

class CCD:
    def __init__(self , ccd_ext , box_x , box_y , ccd_width = 2048 , ccd_height = 4096, plate_scale=0.27):
        self.SIZE_CCD_X = ccd_width
        self.SIZE_CCD_Y= ccd_height
        self.sizex_reg = box_x
        self.sizey_reg = box_y
        self.ext = ccd_ext
        self.plate_scale = plate_scale
        self.reg_nx = int ( floor (self.SIZE_CCD_X / box_x ) )  
        self.reg_ny = int ( floor ( self.SIZE_CCD_Y / box_y ) )
        self.grid = [ [Region() for j in range(self.reg_ny) ] for i in range(self.reg_nx) ]   #Matrix of regions: grid in CCD
        self.means_vec = []

        #print "Instance of CCD created. "
        #print "Size of grid: "
        #print len(self.grid) , len(self.grid[0])
        #print "reg_nx , reg_ny"
        #print self.reg_nx, self.reg_ny

    def get_region_size(self):
        return self.sizex_reg , self.sizey_reg
    def get_ccd_size(self):
        return self.SIZE_CCD_X , self.SIZE_CCD_Y
    def get_median_residual(self):    # Calculates the median of the magnitude of the residuals in each region of CCD
        for i in range(self.reg_nx):
            for j in range(self.reg_ny):
                #print i , j
                meanresx, meanresy = self.grid[i][j].get_mean_residuals() 
                self.means_vec.append(math.sqrt(meanresx**2 + meanresy**2))
        return median(self.means_vec)

    # A CCD should know how to draw its grid with vectors, right? Can I return a "quiver" object?
    # def draw_grid(self):

class focal_plane:
    def __init__(self, ID ,  number_ccds , size_x_ccd , size_y_ccd , sizex_reg_ccd , sizey_reg_ccd):
        self.id = ID 
        self.nccds = number_ccds 
        self.ccd_widthx = size_x_ccd 
        self.ccd_width = size_y_ccd
        self.ccd_box_x = sizex_reg_ccd
        self.ccd_box_y = sizey_reg_ccd
        self.ccd_list = [CCD(i ,sizex_reg_ccd ,sizey_reg_ccd, size_x_ccd , size_y_ccd ) for i in range(number_ccds)]
        # make 'calculate_mean_CCD' part of this class? Then I'd have to modify AST4 to use this focal_plane class...



def plot_clipped_focal_plane ( ccd_list, pp, plot_all_ccds=False, nccds=62, arrow_scale=2, string_residuals_units='pixels', requirement1_tuple=(15,0.055), instrument_name='DES'):
    sizex_reg_temp, sizey_reg_temp = ccd_list[0].get_region_size()    # Use one of the CCDs of the focal_plane vector to get the regions sizes.   
    size_x_ccd , size_y_ccd = ccd_list[0].get_ccd_size()

    mean_CCD = CCD (-3, sizex_reg_temp , sizey_reg_temp)  #Single CCD object: will be the mean of all CCDs
    print "len (ccd_list) ", len (ccd_list)


    camera=instrument_name
    fig1=plt.figure()
    cm = matplotlib.cm.get_cmap('winter')

    for ccd in ccd_list:
     

        x=[]
        y=[]
        meanresvec_xpix=[]
        meanresvec_ypix=[]
        n_objects_vec = []
        clipped_objects_vec = []
        

        for i in range(ccd.reg_nx):
            for j in range(ccd.reg_ny):
                #print i , j
                meanresx, meanresy = ccd.grid[i][j].get_mean_residuals()
                x_local, y_local = i*ccd.sizex_reg + floor(ccd.sizex_reg/2) , j*ccd.sizey_reg + floor(ccd.sizey_reg/2)
                x_focal, y_focal = LocaltoFocal (ccd.ext + 1, x_local, y_local)
                x.append(x_focal)
                y.append(y_focal)
                meanresvec_xpix.append( meanresx )
                meanresvec_ypix.append( meanresy )
                n_objects_vec.append(ccd.grid[i][j].nelements)
                clipped_objects_vec.append(ccd.grid[i][j].get_nclipped())
                
                mean_CCD.grid[i][j].add_clipped(ccd.grid[i][j].get_nclipped()) # Mean of all vector fields; in this case contains all the clipped objects
 
  

        if (plot_all_ccds==True):
            
            if string_residuals_units == 'pixels':
                req_main=requirement1_tuple[1]
                req_sec=requirement1_tuple[0]
                main_label=string_residuals_units
                sec_label='mas'
            elif string_residuals_units == 'mas':
                req_main=requirement1_tuple[0]
                req_sec=requirement1_tuple[1]
                main_label=string_residuals_units
                sec_label='pixels'    
            print "CLIPPED OBJECTS for this CCD", clipped_objects_vec
            #sc=plt.scatter (x, y, c=np.array(clipped_objects_vec), cmap=cm)
            plt.hexbin (x, y, C=np.array(clipped_objects_vec), cmap=cm)
            #plt.colorbar(sc)

            axis([-245,245,-245,245])
            xlabel('Focal Plane X (mm)')
            ylabel('Focal Plane Y (mm)') 
           
            #yticks= [int(i) for i in QP.get_yticks()]
            #xticks= [int(i) for i in QP.get_xticks()]
            #QP.set_yticklabels(yticks, size=7, visible=visible_bool_y)
            #if ccd.ext==0 or ccd.ext==4:
            #    ly=QP.set_ylabel("Pixel")
            #    ly.set_fontsize(10)
            #QP.set_xticklabels(xticks[:-2], size=6.5, visible=visible_bool_x)
            #lx=QP.set_xlabel("Pixel")
            #lx.set_fontsize(10)


    #plt.colorbar(sc)    
    #print "n_objects_vec" , n_objects_vec
    #print "clipped_objects_vec" , clipped_objects_vec
    print "hola clipped"
    if (plot_all_ccds == True):
        fig1.suptitle(' %s clipped objects densityper CCD (mean over all exposures)' %(camera), fontsize=11)
        pp.savefig(fig1)

    return mean_CCD 
             



def calculate_mean_CCD (ccd_list, pp, plot_all_ccds=False, nccds=62, arrow_scale=2, string_residuals_units='pixels', requirement1_tuple=(15,0.055), instrument_name='DES'):
    sizex_reg_temp, sizey_reg_temp = ccd_list[0].get_region_size()    # Use one of the CCDs of the focal_plane vector to get the regions sizes.   
    size_x_ccd , size_y_ccd = ccd_list[0].get_ccd_size()
    mean_CCD = CCD (-3, sizex_reg_temp , sizey_reg_temp)  #Single CCD object: will be the mean of all CCDs
    print "len (ccd_list) ", len (ccd_list)


    camera=instrument_name

    fig1=plt.figure()


    """
    if (plot_all_ccds == True):
       fig1=plt.figure()
       if nccds==62:
           rows_subplot, cols_subplot= 8, 8  #DECam of DES 
       elif nccds==8:
           rows_subplot, cols_subplot= 2, 4   # MOSAIC of BCS, WIF
       else:
           rows_subplot= math.ceil(math.sqrt(nccds))    # provisional
           cols_subplot= math.ceil(nccds/rows_subplot)    
    """
    for ccd in ccd_list:
        #print ccd_list
        #if (plot_all_ccds==True):
            #print ccd.ext
        #    QP=plt.subplot(rows_subplot,cols_subplot, ccd.ext+1)
        #    plt.subplots_adjust(hspace=0.003, wspace=0.003)

        x=[]
        y=[]
        meanresvec_xpix=[]
        meanresvec_ypix=[]
        n_objects_vec = []
        clipped_objects_vec = []
        

        for i in range(ccd.reg_nx):
            for j in range(ccd.reg_ny):
                #print i , j
                meanresx, meanresy = ccd.grid[i][j].get_mean_residuals()
                x_local, y_local = i*ccd.sizex_reg + floor(ccd.sizex_reg/2) , j*ccd.sizey_reg + floor(ccd.sizey_reg/2)
                x_focal, y_focal = LocaltoFocal (ccd.ext + 1, x_local, y_local)
                x.append(x_focal)
                y.append(y_focal)
                meanresvec_xpix.append( meanresx )
                meanresvec_ypix.append( meanresy )
                n_objects_vec.append(ccd.grid[i][j].nelements)
                clipped_objects_vec.append(ccd.grid[i][j].nclipped)
                
                mean_CCD.grid[i][j].add_residuals(meanresx, meanresy) # Mean of all vector fields
 

        ena=1
        print "ccd.ext", ccd.ext
        print "ccd.geat_median_residual()", ccd.get_median_residual()

        

        if (plot_all_ccds==True):
            
            if string_residuals_units == 'pixels':
                req_main=requirement1_tuple[1]
                req_sec=requirement1_tuple[0]
                main_label=string_residuals_units
                sec_label='mas'
            elif string_residuals_units == 'mas':
                req_main=requirement1_tuple[0]
                req_sec=requirement1_tuple[1]
                main_label=string_residuals_units
                sec_label='pixels'    

            if ccd.get_median_residual() < req_main:
                color_arrow='green'
            else:
                color_arrow='red'
            meanresvec_xpix = np.array(meanresvec_xpix)
            meanresvec_ypix = np.array(meanresvec_ypix)
            N= np.sqrt(meanresvec_xpix**2 + meanresvec_ypix**2)
            U, V = meanresvec_xpix / N , meanresvec_ypix / N
            
            q=quiver (x, y, U, V,  units='x', width=0.5, color=color_arrow, scale=arrow_scale, angles='xy') #scale_units ='inches' 
            #type(ccd.get_median_residual())
            text(-0.01*ccd.SIZE_CCD_X, 0.9*ccd.SIZE_CCD_Y, 'median: %.3g pix (%.3g mas)' %(ccd.get_median_residual(), ccd.get_median_residual()*ccd.plate_scale*1000. ), size=6.5)
            
            if ccd.ext == 1:
                QK=quiverkey (q, 0.05, 1.1 , req_main*10,' %.2g %s (%.2g %s)' %(req_main,main_label,req_sec,sec_label) ,  labelpos='E', fontproperties={'size':11})
            
            #l, r, b, t = axis()
            #dx , dy = r-l , t-b
            #q.axis([l-0.05*dx , r + 0.02*dx , b - 0.05*dy , t + 0.02*dy])

            axis([-245,245,-245,245])
            xlabel('Focal Plane X (mm)')
            ylabel('Focal Plane Y (mm)') 
           
            #yticks= [int(i) for i in QP.get_yticks()]
            #xticks= [int(i) for i in QP.get_xticks()]
            #QP.set_yticklabels(yticks, size=7, visible=visible_bool_y)
            #if ccd.ext==0 or ccd.ext==4:
            #    ly=QP.set_ylabel("Pixel")
            #    ly.set_fontsize(10)
            #QP.set_xticklabels(xticks[:-2], size=6.5, visible=visible_bool_x)
            #lx=QP.set_xlabel("Pixel")
            #lx.set_fontsize(10)


        
    #print "n_objects_vec" , n_objects_vec
    #print "clipped_objects_vec" , clipped_objects_vec
    print "hola3"
    if (plot_all_ccds == True):
        fig1.suptitle(' %s astrometric residuals per CCD (mean over all exposures)' %(camera), fontsize=11)
        pp.savefig(fig1)
    return mean_CCD 
    


#if len(mean_vec_field_x)!= len (mean_vec_field_y):
#    print " x, y vectors for mean vector field don't have same size (equal to # of CCDS each one)"
#    sys.exit(1)



# Function that calculates 1D vectors (x, y , value_at_x , value_at_y) in a single CCD 

def vectors_to_plot(ccd , write_ofile=False, ofname='default_vectors_to_plot_AST4.dat'): 
    """
    Function that calculates 1D vectors (x, y , value_at_x , value_at_y , clipped_vec) in a single CCD.
    x, y, value_at_x, value_at_y to be plotted in Quiver
    x, y, clipped_vec to be plotted by scatter (???, for the moment)
    """

    print "begin vectors to plot" 
   #center of regions
    x = []
    y = []
    # Average value per region (single number)
    x_value=[]
    y_value=[]
    number_objects = []
    clipped_vec = []
    counter=0
    if (write_ofile):
        ofile = open(ofname, 'w')
    for i in range(ccd.reg_nx):
        for j in range(ccd.reg_ny):
            mean_meanresx, mean_meanresy = ccd.grid[i][j].get_mean_residuals()    
            x.append(i*ccd.sizex_reg + floor(ccd.sizex_reg/2))    #Should these vecs be calculated inside the class ? 
            y.append(j*ccd.sizey_reg + floor(ccd.sizey_reg/2))
            # Arrows must be flipped DC5(?): map---->     xpix = ywcs  ; ypix = -xwcs    
            x_value.append( mean_meanresx )
            y_value.append( mean_meanresy )
            number_objects.append(ccd.grid[i][j].nelements)  # Should be equal to # CCDS for each region
            clipped_vec.append(ccd.grid[i][j].get_nclipped())
            if (write_ofile):
                oline = '%g %g %g %g \n' %(x[counter] , y[counter] , x_value[counter] , y_value[counter] )
                ofile.write(oline)
            counter+=1    
    #print "vectors inside 'vectors to plot' " , len(x) , len(y) , len(x_value) , len(y_value)        
    if (write_ofile):
        ofile.close()
    return x, y , x_value , y_value , clipped_vec   # Vectors ready to be plotted by Quiver
