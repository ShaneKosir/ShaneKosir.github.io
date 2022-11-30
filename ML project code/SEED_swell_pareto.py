import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import os
from matplotlib import cm
import matplotlib
import numpy as np
import matplotlib.ticker as plticker
from molmass import Formula
import math
from matplotlib.patches import Polygon
from load_data import retPQISData
from scipy.spatial import ConvexHull
from scipy.interpolate import UnivariateSpline
import seaborn as sns
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import matplotlib.cbook as cbook
import matplotlib.image as image
from scipy import stats
out_file='seed_swell_pareto.jpg'

#%% Definitions

def pareto_frontier(Xs, Ys, maxX = True, maxY = True):
    """
    Pareto frontier code
    from:
    http://oco-carbon.com/metrics/find-pareto-frontiers-in-python/
    http://code.activestate.com/recipes/578230-pareto-front/
    
    """
    
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    p_front = [myList[0]]    
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[1] <= p_front[-1][1]:
                p_front.append(pair)
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY
    
    
def molar_mass(formula):
  f=Formula(formula)
  return f.mass

def myround(x,base):
  return int(base*math.ceil(float(x)/base))
 
def function(x,*X):
    m = X[0]
    # print(m)
    # print(X)
    A= X[1]
    B = X[2]
    C = X[3]
    D= X[4]
    return (m*x- (A + B*x+C*x**2+D*x**3))**2

def new_hpf_region_boundary(xx,yy):
    xx= np.sort(xx)
    yy= np.sort(yy)
    yy= yy[::-1]
    spl= UnivariateSpline(xx, yy)
    x_limits = intersect(spl)
    xs = np.linspace(x_limits[0],x_limits[1],10)
    # Plot frontier
    ax1.plot(xs, spl(xs), 'go-', label='Pareto Front, Jet-A, N=%s'%(str(int(len(pf_x)))))

def intersect(spl):
    ABCD = spl.get_coeffs()
    # print(ABCD)
    rho = np.array([0.84,.775])
    xo = np.array([44, 44])
    rho[0] = 0.775
    rho[1] = 0.84
    x = np.zeros(2)
    X = (rho[0],ABCD[0],ABCD[1],ABCD[2],ABCD[3])
    # print(X,type(X))
    x[0] = fsolve(function, xo[0], args=X)
    X = (rho[1],ABCD[0],ABCD[1],ABCD[2],ABCD[3])
    x[1] = fsolve(function, xo[1], args=X)
    print(x)
    return x

def watermark():
    datafile = cbook.get_sample_data('HEAT_watermark_8.png', asfileobj=False)
    print('loading %s' % datafile)
    im = image.imread(datafile)
    #image.thumbnail(datafile,'out.png', scale=0.3, interpolation='bilinear', preview=False)
    im[:, :, -1] = 0.5  # set the alpha channel
    return im
def logoWMOutFigure(fileOut, FN, GS):
    #FOR Image
    # Useful Link
    # https://automatetheboringstuff.com/chapter17/
    from PIL import Image, ImageDraw, ImageFont
    im = Image.open(fileOut)
    # im[:, :, -1] = 0.5 # set the alpha channel
    print(im.info['dpi'])
    print(im.size)
    dpi=im.info['dpi']
    draw = ImageDraw.Draw(im)
    width, height = im.size

    if GS == 1:
        im_LOGO = Image.open(FN).convert("RGBA").convert('LA')
    elif GS == 0:
        im_LOGO = Image.open(FN).convert("RGBA")
    w_LOGO, h_LOGO = im_LOGO.size
    im_LOGO.putalpha(128)
    # im_LOGO[:, :, -1] = 0.5 # set the alpha channel
    ## scaling of logo
    scale = 0.40
    print(int(w_LOGO*scale))
    print(int(h_LOGO*scale))
    im_LOGO = im_LOGO.resize((int(w_LOGO*scale), int(h_LOGO*scale)))
    out_file_wm= out_file[:-6]+'_wm.png'
    # Origin is at the top left
    #im.paste(im_LOGO, (width - int(w_LOGO*scale), height - int(h_LOGO*scale)))
    im.paste(im_LOGO, (int((width - int(w_LOGO*scale))*.65), int((height - int(h_LOGO*scale))*0.85)))
    im.save(fileOut.replace('.png', out_file), dpi = dpi)
    return out_file_wm

def conv_jet_names():
    jet_names=['Jet A-1','Jet A','JP-5','JP-8']
    # Position fuel labels
    for ii in range(len(jet_names)):
      if jet_names[ii]=='Jet A-1':
          ax1.text(43.46, jet_ed[ii]-.05, jet_names[ii], size='x-small',zorder=10, color='k') # jet_se[ii]+.03
      if jet_names[ii]=='Jet A':
          ax1.text(42.85, jet_ed[ii]-.05, jet_names[ii], size='x-small',zorder=10, color='k')# jet_se[ii]-.11
      if jet_names[ii]=='JP-5':
          ax1.text(42.9, jet_ed[ii]-.05, jet_names[ii], size='x-small',zorder=10, color='k')# jet_se[ii]-.095
      if jet_names[ii]=='JP-8':
          ax1.text(43.47, jet_ed[ii]-.05, jet_names[ii], size='x-small',zorder=10, color='k') # jet_se[ii]+.03

def hpf_shading():
    ### HPF Range
    fill_list=[[43.204,.84*(43.20125921)],[50,.84*50],[50,34.9072202], \
               [43.59,34.9072202],[43.54,35.2674],[43.3,35.50599999999999], \
               [43.204,35.85]]
    ax1.add_patch(Polygon(fill_list,closed=True,fill=True,lw=.3,facecolor='green',zorder=0,alpha=.3))
    ax1.text(44.5,36,'HPF Region', size='small', color='blue',horizontalalignment='center',
        verticalalignment='center',zorder=10)

def low_int_hpf (x, ED):
    """
    Functions for FSOLVE.
    Compares PQIS Outer Region to to
    lower floor (SE*0.775) = ED
    """
    return ED - lower_lf(x)

def int_line_pqis(x, ED):
    """
    Returns intercept of outerp PQIS
    line with the PQIS intercept with lower bound
    """
    return ED - fn_pqis_interp(x)
        
def hpf_shading_r2():
    import copy
    
    # Starting bounds from upper right region working around clockwise
    st_list = [[50,0.84*50],[50, 0.775*50]]
    
    # PT right by the lower bound and intercept of Jet-A PQIS Scaled Data
    x_pt_RIGHT = fsolve(low_int_hpf, 43.5, args=(REF_ED))[0]
    y_pt_RIGHT = x_pt_RIGHT*0.775
    st_list = [*st_list, [x_pt_RIGHT, y_pt_RIGHT]]
    
    # POINT of intercept with outliers of PQIS Database.
    global fn_pqis_interp, fn_pqis_interp_rev
    fn_pqis_interp = interp1d(pareto_x, pareto_y, fill_value='extrapolate')
    fn_pqis_interp_rev = interp1d(pareto_y, pareto_x, fill_value='extrapolate')
    
    #int_pt_low_SE = fsolve(int_line_pqis, 43.1, args=(y_pt_RIGHT))[0]
    #int_pt_low_ED = fn_pqis_interp(int_pt_low_SE)
    
    int_pt_low_SE = fn_pqis_interp_rev(y_pt_RIGHT)
    int_pt_low_ED = y_pt_RIGHT
    
    
    st_list = [*st_list, [int_pt_low_SE, int_pt_low_ED]]
    #ax1.plot(int_pt_low_SE, fn_pqis_interp(int_pt_low_SE), color='orchid', marker='*')
    
    # POINTS on PQIS OUTLIER 
    pareto_x_2 = [i for i in pareto_x if i > REF_SE and pareto_y[pareto_x.index(i)]> REF_ED]
    pareto_y_2 = [pareto_y[pareto_x.index(i)] for i in pareto_x_2]
    pf = []
    for i in reversed(range(len(pareto_x_2))):
        pf = [*pf, [pareto_x_2[i], pareto_y_2[i]]]
    st_list = [*st_list, *pf]
    #ax1.plot(pareto_x_2, pareto_y_2, color='r', marker='s')
    
    # POINT where Ref SE intercepts on PQIS Outliers
    st_list = [*st_list, [REF_SE, fn_pqis_interp(REF_SE)]]
    #ax1.plot(REF_SE, fn_pqis_interp(REF_SE), color='g', marker='*')
    
    # Point where PQIS SE Line, and then to the upper limit 
    st_list = [*st_list, [REF_SE, REF_SE*0.84]]
    st_list = [*st_list, [50., 50.*0.84]]
    
    ### HPF Range
    # fill_list=[[43.204,.84*(43.20125921)],[50,.84*50],[50,34.9072202], \
               # [43.59,34.9072202],[43.54,35.2674],[43.3,35.50599999999999], \
               # [43.204,35.85]]
    ax1.add_patch(Polygon(st_list,closed=True,fill=True,lw=.3,facecolor='blue',zorder=0,alpha=.3))
    ax1.text(44.6,36.2,'HPF Region', size='small', color='blue',horizontalalignment='center',
        verticalalignment='center',zorder=10)        

def jet_A_average_lines():
    ## Jet A Lines average energy density and specific energy
    conv_xs1=[43.20125921,43.20125921]
    conv_ys1=[34.9072202,.84*43.20125921]
    ax1.plot(conv_xs1,conv_ys1,'red',linestyle='--',linewidth=1,zorder=1)
    conv_xs2=[43.20125921,50]
    conv_ys2=[34.9072202,34.9072202]
    ax1.plot(conv_xs2,conv_ys2,'red',linestyle='--',linewidth=1,zorder=1)

def color_map():
    # Adding colors based on colorbar to plot
    global cmap_name
    norm=matplotlib.colors.Normalize(vmin=0,vmax=1000) # Normalize the colors from [0,1] based on a min and max value
    cmap_name='plasma_r' # _r to reverse
    cmap=matplotlib.cm.get_cmap(cmap_name)
    #color_use=cmap(norm(xxx))

# COLORBAR
def colorbar_addition():        
    # Add in the color bar to the RHS of the figure...
    # for outside the figure...
    #ax1 = fig.add_axes([0.97, 0.15, 0.025, 0.7]) # L,B,W,H
    # for inside the figure and DCN...
    #ax1 = fig.add_axes([0.2, 0.15, 0.025, 0.5]) # L,B,W,H
    # for inside the figure and RI....
    ax5 = fig.add_axes([0.81, 0.13, 0.015, 0.2]) # L,B,W,H
    # Calculate number of ticks needed
    # Arbitrarily picked as the length of the fuels, but can just be any number.
    # num_ticks = len(x1a) 
    num_ticks = 3

    ### Method 1: Using actual values and using as colorbar vars
    # This metod uses the actual values as color bar ticks..
    #tick_loc = np.zeros(num_ticks)
    #for i in range (len(x1a)):
    #    tick_loc[i] = norm(x1a[i])
    #labels = ['{:.1f}'.format(x) for x in x1a]

    ### Method 2: Generic scale appropriate 
    ### for the axis involving color and their associated values.
    ticks = np.linspace(0, 1, num_ticks)
    # Making sure the axis for the color bar aligns to the normalized maxs and mins above.
    labels = np.linspace(0, 1000, num_ticks)
    #labels = ['{:.d}'.format(x) for x in labels]
    labels = [0,500,1000]
    # Adding in the colorbar...
    cbar = matplotlib.colorbar.ColorbarBase(ax5,cmap=cmap_name,
                                     orientation='vertical')#norm=norm, 
    ## If using method 1, uncomment the proceediing line
    #cbar.set_ticks(tick_loc[::1])
    ## If using method 2, uncomment the proceeding line
    cbar.set_ticks(ticks[::1])
    # Add in appropriate labels as per above
    print(labels[::1])
    cbar.set_ticklabels(labels[::1]) 
    cbar.ax.tick_params(labelsize=10) 
    cbar.ax.tick_params(
        axis='y',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        right=False)
    cbar.ax.minorticks_off()

    # Add label for colorbar. Currently label not rotated, but can be
    # cbar.set_label(r'$DCN$', rotation=00, labelpad =15)#20
    cbar.set_label(r'$\%V$', rotation=00, labelpad =10, size=10)#20
    # Place colorbar label in specific position.
    ax5.yaxis.set_label_position("left")
    ####

def scaling_conv_fuels():
    global jet_se, jet_ed
    # # SCALING
    jet_A_val_PQIS_SE=43.20125921
    jet_A_val_PQIS_ED=34.9072202
    jet_A_val_PQIS_DENS = jet_A_val_PQIS_ED/jet_A_val_PQIS_SE
    
    global A2_SE, A2_ED
    A2_SE = 43.06
    A2_ED = 0.803*43.06
    
    global res_SE
    # Caling first on the residual
    res_SE = jet_A_val_PQIS_SE - A2_SE
    jet_dens = [jet_ed[i]/jet_se[i] for i in range(len(jet_se)) if len(jet_se) == len(jet_ed)]
    print (jet_dens)
    jet_se = [i-res_SE for i in jet_se]
    jet_ed = [jet_se[i]*jet_dens[i] for i in range(len(jet_se))]
    
    
    global REF_SE, REF_ED
    REF_SE = jet_A_val_PQIS_SE - res_SE
    REF_ED = jet_A_val_PQIS_DENS*REF_SE
    
    
def pt1_plot():
    """
    Part 1 of plot
    """
    ## Data for conventional fuel boxes
    jet_names=['Jet A-1','Jet A','JP-5','JP-8']
    global jet_se, jet_ed
    jet_se=[43.26950331,43.20125921,43.14795266,43.25658127] 
    jet_ed=[34.41616224,34.9072202,35.05963253,34.65312478]

    
    scaling_conv_fuels()
    
    #%% SNL Data
    snl_se=[44.3,44.3,45,43.54]
    snl_ed=[35.6,37.8,38.7,37.01]
    snl_names=['limonane','pinane','cis-carane','caryophyllane']

    #%% GT molecules
    gt_se1=44.9
    gt_ed1=36.1
    gt_se2=44.9
    gt_ed2=37.7

    #%% Plot setup
    global fig, ax1, ax2, ax3
    fig,ax1=plt.subplots()
    
    ax2=ax1.twinx()
    ax3=ax1.twiny()

    #conv_fuel_molecular_families()

    #ax1.plot(43.4,.827*43.4, color = 'g', marker='*',mec='k',ls='None', ms=10, label='_nolegend_')#'Shell IH'+'$^2$'

    
def savenopen_file(out_file2):
    plt.savefig(out_file2,bbox_inches='tight',dpi=350,transparent=False)
    plt.show()
    #out_file_wm=logoWMOutFigure(out_file, logoFN, 0)#
    os.startfile(out_file2,'open') 


def plot_conv_pareto_front_PQIS():
    global pareto_x, pareto_y
    # Plot Conventional Pareto Front
    pareto_x=[43.1,43.3,43.54,43.6,43.6,43.6] # 43.1,43.6
    pareto_y=[36.204,35.50599999999999,35.2674,34.88,34.444,34.008] # 35.77,33.79
    pareto_dens = [pareto_y[i]/pareto_x[i] for i in range(len(pareto_x)) if len(pareto_x) == len(pareto_y)]
    
    if "res_SE" in globals():
        print("TESTTESTTEST")
        pareto_x = [i-res_SE for i in pareto_x]
        pareto_y = [pareto_x[i]*pareto_dens[i] for i in range(len(pareto_x))]
        # pareto_y = [i-res_ED for i in pareto_y]
        # for i in range (len(pareto_y)):
            # if pareto_y[i] < pareto_x[i]*0.775:
                # pareto_y[i] = pareto_x[i]*0.775
    
    ax1.plot(pareto_x,pareto_y,color='blue',
             label='"Best-case" Conventional Fuels',
             solid_capstyle='round')

def pt2_plot():
    ### Insert Watermark jsh
    # im = watermark()
    # fig.figimage(im, 2250, 500, zorder=3,resize=False)
    
    ## PQIS Data

    #data_sc = retPQISData('tblAFKERO.xlsx')
    data_sc =np.load('jittered_data.npy')
    print(data_sc)
    
    # Rescaling data to keep density constant from database
    # But to switch the density and recalculate the energy density.
    if "res_SE" in globals():
        for i in range(len(data_sc)):
            dens = data_sc[i,1]/data_sc[i,0]
            data_sc[i,0] = data_sc[i,0] - res_SE
            data_sc[i,1] = data_sc[i,0] * dens
        print(type(data_sc))
    
    #ax1.plot(data_sc['SE'], data_sc['ED'], marker='o', color='b', alpha=0.6, ms=6, ls='None', label='_nolegend_')
    #pf_x1, pf_y1 = pareto_frontier(data_sc['SE'].values, data_sc['ED'].values, maxX = True, maxY = True)
    #ax1.scatter(data_sc['SE'],data_sc['ED'])
    #xx = np.array(pf_x)
    #yy = np.array(pf_y)


    
    plot_conv_pareto_front_PQIS()
    


    # Plot Gradient
    sns.kdeplot(data_sc[:,0],data2=data_sc[:,1],ax=ax1,zorder=0)


    # Set axis limits  
    ## Specialty molecules
    # x_lim=[41.5,48.5]
    # y_lim=[32,55]
    ## Coventional Fuel Molecules
    # x_lim=[40.5,45.1]
    # y_lim=[30,40]
    ## Zoomed in Conventional Fuels
    x_lim=[39.75,45.25] # uniform axes
    y_lim=[30,40.5]
    ax1.set_xlim(x_lim)
    ax1.set_ylim(y_lim)

    
    
    # REFerences for lines
    # Commented out the previous ones, because scaling 
    # code above should fix this.
    if "REF_SE" and "REF_ED" not in globals():
        global REF_SE, REF_ED
        REF_SE = 43.20125921
        REF_ED = 34.9072202
    
    # Reference with JEt A, if not comment that out
    #REF_SE_2 = REF_SE
    #REF_ED_2 = REF_ED
    REF_SE_2 = A2_SE
    REF_ED_2 = A2_ED
    
    # Scaling
    # Set twin axis limits
    V_per_min=(y_lim[0]-(REF_ED_2))/(REF_ED_2)*100
    V_per_max=(y_lim[1]-(REF_ED_2))/(REF_ED_2)*100
    ax2.set_ylim(V_per_min,V_per_max)
    m_per_min=(x_lim[0]-(REF_SE_2))/(REF_SE_2)*100
    m_per_max=(x_lim[1]-(REF_SE_2))/(REF_SE_2)*100
    ax3.set_xlim(m_per_min,m_per_max)

    # Labels
    FS_label='x-large'
    ax1.set_xlabel('Specific Energy, MJ/kg',size=FS_label)
    ax1.set_ylabel('Energy Density, MJ/L',size=FS_label)
    ax2.set_ylabel('Energy Density, % Diff.',size=FS_label)
    ax3.set_xlabel('Specific Energy, % Diff.',size=FS_label)

    # Specific energy spec limit
    xrange_se=np.arange(42.8,x_lim[1]+1,.5)
    se_spec_limit_lower=ax1.plot((42.8,42.8),y_lim,color='red',zorder=0,linewidth=.5)
    #ax1.axvspan(42.8,max_x[1],facecolor='red',alpha=.07,zorder=0)

    ## Spec. Limits
    ax1.text(42.75,37.9,'Specific Energy Limit',
             size='small',color='red',horizontalalignment='right')
    # ax1.text(42.9,32,'Jet A Energy Density \n Spec Range', size='x-small', color='blue')

    # Density limit
    xrange=np.arange(42.8,x_lim[1]+1,.5)
    yrange_low=[xrange[ii]*.775 for ii in range(len(xrange))]
    yrange_high=[xrange[ii]*.84 for ii in range(len(xrange))]
    # xrange2=np.arange(jet_se_max,x_lim[1]+1,.5)
    # yrange_low2=[xrange2[ii]*.775 for ii in range(len(xrange2))]
    # yrange_high2=[xrange2[ii]*.84 for ii in range(len(xrange2))]
    density_lower=ax1.plot(xrange,yrange_low,color='blue',zorder=0,linewidth=.5)
    density_upper=ax1.plot(xrange,yrange_high,color='blue',zorder=0,linewidth=.5)
    ax1.fill_between(xrange,yrange_low,yrange_high,facecolor='blue',alpha=.1,zorder=0)

    ### 5% line
    xs=[43.42185921000733,45.204859210066516] # Values from 5% Limits.py
    ys=[36.47436173640615,35.03376588780155]
    #ax1.plot(xs,ys,'black',linestyle='--',label='5% Cumulative\nSE/ED Increase',zorder=0)

    
  
    ## HPF Shading
    #hpf_shading()
    if HPF_LOG == 1:
        hpf_shading_r2()
    
    ## Jet A Lines average energy density and specific energy
    #jet_A_average_lines()
    ## Insert the names of the Conventional fuels
    #conv_jet_names()
    # Minor Ticks
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', direction='out')
    ax1.minorticks_on()
    ax1.tick_params(axis='x', which='minor', direction='out')
    ax2.minorticks_on()
    ax2.tick_params(axis='y', which='minor', direction='out')
    loc=plticker.MultipleLocator(base=.5) # x
    ax1.xaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=1) # y
    ax1.yaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=2.5) # y2
    ax2.yaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=1) # x2
    ax3.xaxis.set_minor_locator(loc)

    # Molecule legend
    legend2=ax1.legend(loc='lower left',prop={'size':'x-small'},
                        ncol=1,columnspacing=1)
    legend2.legendHandles[3].set_facecolor('darkgrey')
    legend2.legendHandles[4].set_facecolor('darkgrey')
    legend2.legendHandles[5].set_facecolor('darkgrey')
    legend2.legendHandles[6].set_facecolor('darkgrey')
    legend2.legendHandles[7].set_facecolor('darkgrey')
    legend2.legendHandles[8].set_facecolor('darkgrey')
    legend2.get_frame().set_facecolor('none')
    legend2.get_frame().set_edgecolor('none')


def low_int (x):
    """
    Functions for FSOLVE.
    Compares SEED Cumulative SEED LINE to
    lower floor (SE*0.775) = ED
    """
    return SEED_PCT_line(x) - lower_lf(x)
def low_int_floor (x):
    """
    Functions for FSOLVE.
    Compares SEED Cumulative SEED LINE to
    lower floor of HPF Region
    Based on the PQIS Information
    
    """
    return SEED_PCT_line(x) - lower_floor_pt(x)
def up_int (x):
    """
    Functions for FSOLVE.
    Compares SEED Cumulative SEED LINE to
    upper ceiling (SE*0.84) = ED
    """
    return SEED_PCT_line(x) - upper_lf(x)    

def get_CFIT_DENS_RANGE ():
    """Curve fit and interpolations for
    upper and lower regions of Jet Fuel"""
    x_FIT = [42.8, 46.]
    y_FIT_LOW = [i*0.775 for i in x_FIT]
    y_FIT_UP = [i*0.84 for i in x_FIT]
    
    global lower_lf, upper_lf
    lower_lf = interp1d(x_FIT, y_FIT_LOW, fill_value='extrapolate')
    upper_lf = interp1d(x_FIT, y_FIT_UP, fill_value='extrapolate')

def plot_pct_inc(PCT, col_r):
    """
    Plotter for the cumulative SE/ED Line.
    This is based on a specified percentage,
    expressed as a decimal.
    
    References are PQIS Data based on the Jet Fuel
    
    """
    
    ## ROB LINES for PCT
    # Reference from PQIS
    SE_REF = 43.20125921
    ED_REF = 34.9072202
    if "REF_SE" and "REF_ED" in globals():
#        print('ROB, youre the man')
        SE_REF = REF_SE
        ED_REF = REF_ED
    
    # Lower Floor of HPF Region Calculation
    x_FIT = [42.8, 46.]
    y_FIT = [ED_REF, ED_REF]
    global lower_floor_pt
    lower_floor_pt = interp1d(x_FIT, y_FIT, fill_value='extrapolate')
    
    # Calculate line for the SE/ED based on percentages
    PCT = float(PCT) # Expressed as decimal
    UP_PT_SE = SE_REF
    UP_PT_ED = ED_REF*(1+PCT)
    LOW_PT_SE = SE_REF*(1+PCT)
    LOW_PT_ED = ED_REF

    # Create Line function of the SE/ED Cumulative Line
    x_LINE = [UP_PT_SE, LOW_PT_SE]
    y_LINE = [UP_PT_ED, LOW_PT_ED]
    global SEED_PCT_line, SEED_PCT_line_MOD
    SEED_PCT_line = interp1d(x_LINE, y_LINE, fill_value='extrapolate')
    
    # Bound the line by the barriers of Ceiling and Floor of Jet-A range.
    x_pt_SE_high = fsolve(up_int, 43.5)[0]
    x_pt_ED_high = x_pt_SE_high*0.84

    x_pt_SE_low = fsolve(low_int, 44.5)[0]
    x_pt_ED_low = x_pt_SE_low*0.775
    
    if x_pt_ED_low < ED_REF:
        """IF right point of the line is intersecting
        in the area below the HPF region, this will adjust the line"""
        print('Adjusting one of the pts on Cumulative Line for HPF Region...')
        x_pt_SE_low = fsolve(low_int_floor, 44.5)[0]
        x_pt_ED_low = float(SEED_PCT_line(x_pt_SE_low))

    # Assembling Line for plotting
    x_LINE_MOD = [x_pt_SE_high, x_pt_SE_low]
    y_LINE_MOD = [x_pt_ED_high, x_pt_ED_low]
    SEED_PCT_line_MOD = interp1d(x_LINE_MOD, y_LINE_MOD)
    
    print(x_LINE_MOD)
    print(y_LINE_MOD)
    
    
    # ax1.plot(x_LINE,y_LINE,
        # 'r',linestyle='-.',label='{:.1f}% Cumulative\nSE/ED Increase'.format(PCT*100),zorder=-5)
    ax1.plot(x_LINE_MOD,y_LINE_MOD,
        color=col_r,linestyle='--',label='{:.1f}% Cumulative\nSE/ED Increase'.format(PCT*100),zorder=0)

def plot_pct_inc_A2(PCT, col_r):
    """
    Plotter for the cumulative SE/ED Line.
    This is based on a specified percentage,
    expressed as a decimal.
    
    References are Jet Fuel (POSF 10325)
    
    """
    
    # Reference from PQIS
    PQIS_SE_REF = 43.20125921
    PQIS_ED_REF = 34.9072202
    
    if "REF_SE" and "REF_ED" in globals():
#        print('ROB, youre the man')
        PQIS_SE_REF = REF_SE
        PQIS_ED_REF = REF_ED
    
    ## ROB LINES for PCT
    # POSF 10325 Values
    SE_REF = 43.06
    ED_REF = 0.803*43.06
    
    # Lower Floor of HPF Region Calculation
    x_FIT = [42.8, 46.]
    y_FIT = [PQIS_ED_REF, PQIS_ED_REF]
    global lower_floor_pt
    lower_floor_pt = interp1d(x_FIT, y_FIT, fill_value='extrapolate')
    
    # Calculate line for the SE/ED based on percentages
    PCT = float(PCT) # Expressed as decimal
    UP_PT_SE = SE_REF
    UP_PT_ED = ED_REF*(1+PCT)
    LOW_PT_SE = SE_REF*(1+PCT)
    LOW_PT_ED = ED_REF

    # Create Line function of the SE/ED Cumulative Line
    x_LINE = [UP_PT_SE, LOW_PT_SE]
    y_LINE = [UP_PT_ED, LOW_PT_ED]
    global SEED_PCT_line, SEED_PCT_line_MOD
    SEED_PCT_line = interp1d(x_LINE, y_LINE, fill_value='extrapolate')
    
    # Bound the line by the barriers of Ceiling and Floor of Jet-A range.
    x_pt_SE_high = fsolve(up_int, 43.5)[0]
    x_pt_ED_high = x_pt_SE_high*0.84
    
    
    if x_pt_SE_high < SE_REF:
        """
        Adjusting the top point (Left point) and bounding to
        the limit of the reference
        """
        print('Adjusting upper point')
        x_pt_SE_high = SE_REF
        x_pt_ED_high = SEED_PCT_line(SE_REF)

    x_pt_SE_low = fsolve(low_int, 44.5)[0]
    x_pt_ED_low = x_pt_SE_low*0.775
    
    if x_pt_ED_low < PQIS_ED_REF:
        """IF right point of the line is intersecting
        in the area below the HPF region, this will adjust the line"""
        print('Adjusting one of the pts on Cumulative Line for HPF Region...')
        x_pt_SE_low = fsolve(low_int_floor, 44.5)[0]
        x_pt_ED_low = SEED_PCT_line(x_pt_SE_low)

    # Assembling Line for plotting
    x_LINE_MOD = [x_pt_SE_high, x_pt_SE_low]
    y_LINE_MOD = [x_pt_ED_high, x_pt_ED_low]

    SEED_PCT_line_MOD = interp1d(x_LINE_MOD, y_LINE_MOD)
    
    #print(x_LINE_MOD)
    #print(y_LINE_MOD)
    
#    ax1.plot(x_LINE_MOD,y_LINE_MOD,
#        color=col_r,linestyle='--',label='{:.1f}% Cumulative\nSE/ED Increase to Jet A'.format(PCT*100),zorder=0)


#def line_shift():
    # NJFCP_fuels_SE = [43.24, 43.06, 42.88]
    # NJFCP_fuels_ED = [0.780*43.24, 0.803*43.06, 0.827*42.88]

    # slope_njfcp, intercept_njfcp, r_value_njfcp, p_value_njfcp, std_err_njfcp = stats.linregress(NJFCP_fuels_SE, NJFCP_fuels_ED)
    # global xs_njfcp, ys_njfcp
    # xs_njfcp = [max(NJFCP_fuels_SE), min(NJFCP_fuels_SE)]
    # ys_njfcp = [slope_njfcp*i+intercept_njfcp for i in xs_njfcp]



    # jet_se_PQIS_DB=[43.26950331,43.20125921,43.14795266,43.25658127] 
    # jet_ed_PQIS_DB=[34.41616224,34.9072202,35.05963253,34.65312478]           
    # #ax1.scatter(jet_se_PQIS_DB,jet_ed_PQIS_DB,facecolor='orange',edgecolor='black',
    # #            marker='d',label='Jet-Fuels, Test',linewidth=.5,s=60,zorder=10)
    # slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(jet_se_PQIS_DB, jet_ed_PQIS_DB)
    # global xs_1, ys_1
    # xs_1 = [max(jet_se_PQIS_DB), min(jet_se_PQIS_DB)]
    # ys_1 = [slope1*i+intercept1 for i in xs_1]
    # print (xs_1)
    # print (ys_1)
    # ax1.plot(xs_1,
            # ys_1,
            # color='orange',mec='black',
            # marker='d',label='Jet-Fuels, Test',linewidth=.5,ms=5,zorder=50)     
    # ax1.plot(xs_njfcp,
            # ys_njfcp,
            # color='green',mec='black',
            # marker='d',label='Jet-Fuels, NJFCP, Test',linewidth=.5,ms=5,zorder=50)       
def alt_fuels():
    """
    refer to the below definition
    "plot_alternative_fuels"
    
    """
    it_fuels = [
        'Jet-A',
        'GEVO ATJ',
        'FT-SPK (Syntroleum GTL)',
        'HEFA-SPK (Sasol IPK)',
        'Farnesane',
        'Lanzatech ETJ', 
        ] #Shell iH2

    shapes = [
            'o',
            'd',
            'p',
            'h',
            '>',
            'X']
    despt = np.array([100.])
    ms = 40.
    se_fuels = [43.9,44.1,43.7,44,43.9, 43.1]
    density_fuels = []
    for i in range(len(it_fuels)):
        se, ed, vol, fls, lbl = gbd('Shell iH2', it_fuels[i])
        # ax1.scatter(se, ed, s=20, c=cmap(norm(vol)), marker=shapes[i], edgecolor='None', label=it_fuels[i],)# linewidth=0.5,
        
        if i == 0:
            ax1.scatter(se[0], ed[0], s=ms, c='g', marker='*', edgecolor='k', zorder=10, linewidth=0.5, label='Shell IH$^2$')
        mask = np.isin(vol, despt)
        se_msk = se[mask]
        ed_msk = ed[mask]
        vol_msk = vol[mask]
        ax1.scatter(se_msk, ed_msk, s=ms, c=cmap(norm(vol_msk)), marker=shapes[i], edgecolor='k', zorder=10, linewidth=0.5, label=it_fuels[i])

def plot_NJFCP_fuels():
    colors=cm.tab10(np.linspace(0,1,10))
    cs=60
    ax1.scatter(42.88,0.827*42.88,facecolor=colors[3][0:4],edgecolor='black',
                    marker='h',label='Conventional Jet Fuels',linewidth=.5,
                    zorder=10,s=cs) #JP-5 (POSF 10289)
    ax1.scatter(43.06,0.803*43.06,facecolor=colors[3][0:4],edgecolor='black',
                    marker='h',label='__nolegend__',linewidth=.5,
                    zorder=10,s=cs) #Jet A (POSF 10325)
    ax1.scatter(43.24,0.780*43.24,facecolor=colors[3][0:4],edgecolor='black',
                    marker='h',label='__nolegend__',linewidth=.5,
                    zorder=10,s=cs) #JP-8 (POSF 10264)

def plot_alternative_fuels():
    shapes = [
            '*',
            'd',
            'p',
            'h',
            '>',
            'X',
            ]
    df = pd.read_excel('Approved_Fuels_Basic_Info.xlsx')
    print('Shape{}'.format(str(df.shape[0])))
    for i in range(df.shape[0]):
        if df['Fuel'].values[i] == 'Shell IH2':
            fc = 'orange'
        else:
            fc = 'forestgreen'
        ax1.scatter(df['LHV (MJ/kg)'].values[i],
                    df['Density (15C)'].values[i]*df['LHV (MJ/kg)'].values[i],
                    facecolor=fc,edgecolor='black',
                    marker=shapes[i],label=df['Fuel Label'].values[i],
                    linewidth=.5,s=65,zorder=10)
      
        
# 'SEED_pareto_Scenarios1to6_r1.png'




pt1_plot()
color_map()
#colorbar_addition()
get_CFIT_DENS_RANGE()
plot_NJFCP_fuels()

import copy
PREV_VALUE_SE = copy.deepcopy(REF_SE)
PREV_VALUE_ED = copy.deepcopy(REF_ED)
REF_SE = A2_SE
REF_ED = A2_ED

plot_pct_inc_A2(0.04,'black')    
# conv_jet_names()


global HPF_LOG
HPF_LOG = 1
# ax1.axhline(y=REF_ED, color='g', label = '__nolegend__')
# ax1.axvline(x=REF_SE, color='g', label = '__nolegend__')
# ax1.axhline(y=A2_ED, color='b', label = '__nolegend__')
#ax1.axvline(x=A2_SE, color='b', label = '__nolegend__')

#REF_SE = PREV_VALUE_SE
#REF_ED = PREV_VALUE_ED
#plot_alternative_fuels()
#alt_fuels()

### Conventional molecules
# Import neat swell data
df_neat_swell=pd.read_excel('all_molecules_nitrile_downselected.xlsx')
neat_swell=df_neat_swell['Neat Swell']
neat_molecule=df_neat_swell['Molecule']

# Import blend swell data
df_blend_seed_swell=pd.read_csv('Nitrile 10000 mfs all props parsed.csv')
blend_se=df_blend_seed_swell['Specific Energy, MJ/kg']
blend_ed=df_blend_seed_swell['Energy Density, MJ/L']
blend_swell=df_blend_seed_swell['Volume Swell']

# Import conventional fuel dataframes
df_conv=pd.read_excel('all_molecules.xlsx')

# Only include molecules with property values
df_conv=df_conv.dropna(subset=['Specific Energy, MJ/kg','Density @ T1, g/cc'])

# Reset indices
df_conv=df_conv.reset_index(drop=True)

# Remove specialty molecules
df_conv=df_conv.iloc[0:54]

# Reset indices
df_conv=df_conv.reset_index(drop=True)

# Convert pandas dataframe to numpy arrays for plotting
cf_se=df_conv['Specific Energy, MJ/kg'].tolist()
number_of_fuels=len(cf_se)
cf_rho=df_conv['Density @ T1, g/cc'].tolist()
cf_ed=[cf_se[ii]*cf_rho[ii] for ii in range(number_of_fuels)]
colors=cm.tab10(np.linspace(0,1,10))

# Get swell molecule indices
conv_molecule=df_conv['Molecule']
swell_overlap=pd.DataFrame([conv_molecule[conv_molecule == neat_molecule[ii]] for ii in range(len(neat_molecule))]).T
"""Molecules need to be in same order in both spreadsheets, should fix this"""
overlap_index_neat=list(swell_overlap.index.values)

# Molecular group indices
n_index=list(range(0,11))
iso_index=list(range(11,23))
monocyclo_index=list(range(23,33))
dicyclo_index=list(range(33,37))
aro_index=list(range(37,53))

# Legend indices
n_legend=[]
iso_legend=[]
monocyclo_legend=[]
dicyclo_legend=[]
aro_legend=[]

# Combine neat and blend swell
#overlap_index_blend=list(range(len(cf_se),len(cf_se)+len(blend_swell)))
#overlap_index_all=np.array(overlap_index_neat+overlap_index_blend)
#all_swell=np.array(neat_swell.append(blend_swell))

# Normalize colors
xxx=np.array(neat_swell) # molecular weights to normalize
min_norm=min(xxx)
max_norm=max(xxx)
norm=matplotlib.colors.Normalize(vmin=min_norm,vmax=max_norm) # Normalize the colors from [0,1] based on a min and max value
cmap_name='Greys' # _r to reverse
cmap=matplotlib.cm.get_cmap(cmap_name)
color_use=cmap(norm(xxx))

# Reassign indices to numpy array
color_use_pd=pd.DataFrame(color_use)
color_use_pd=color_use_pd.set_index(np.array(overlap_index_neat))

# Add colorbar to plot
ax4=fig.add_axes([.58,.18,.22,0.025]) # L,B,W,H
num_ticks=4
ticks=np.linspace(0,1,num_ticks)
labels=np.linspace(min_norm,max_norm,num_ticks)
labels=[round(x,1) for x in labels]
cbar=matplotlib.colorbar.ColorbarBase(ax4,cmap=cmap_name,orientation='horizontal')
cbar.set_ticks(ticks[::1])
cbar.set_ticklabels(labels[::1])
cbar.ax.tick_params(labelsize='x-small')
cbar.set_label('Volume Swell, % v/v',rotation=00,labelpad=-32,size='small')
ax4.yaxis.set_label_position('right')
ax4.tick_params(axis='x',which='minor',bottom=False)

# Swell limits
def min_max_norm(x,minimum,maximum):
    norm_value=(x-minimum)/(maximum-minimum)
    return norm_value
    
swell_lower_limit=min_max_norm(3.7,min_norm,max_norm)
swell_upper_limit=min_max_norm(17.4,min_norm,max_norm)
ax4.scatter(swell_lower_limit,.5,marker='|',color=colors[3][0:4])

# Plot definition
def molecule_plot(indices,overlap_indices,legend_list,legend_label,se_list,
                  ed_list,facecolor_swell,facecolor_no_swell,edgecolor,
                  marker,linewidth,zord,size):
    for ii in indices:
        if ii in overlap_indices and len(legend_list) == 0:
            ax1.scatter(se_list[ii],ed_list[ii],facecolor=np.array(color_use_pd.loc[[ii]]),
                        edgecolor=edgecolor,marker=marker,label=legend_label,
                        linewidth=linewidth,zorder=zord,s=size)
            legend_list.append(ii)
        elif ii in overlap_indices and len(legend_list) != 0:
            ax1.scatter(se_list[ii],ed_list[ii],facecolor=np.array(color_use_pd.loc[[ii]]),
                        edgecolor=edgecolor,marker=marker,label='__nolegend__',
                        linewidth=linewidth,zorder=zord,s=size)
        elif ii not in overlap_indices and len(legend_list) == 0:
            ax1.scatter(se_list[ii],ed_list[ii],facecolor=facecolor_no_swell,
                        edgecolor=edgecolor,marker=marker,label=legend_label,
                        linewidth=linewidth,zorder=zord,s=size)
            legend_list.append(ii)
        elif ii not in overlap_indices and len(legend_list) != 0:
            ax1.scatter(se_list[ii],ed_list[ii],facecolor=facecolor_no_swell,
                        edgecolor=edgecolor,marker=marker,label='__nolegend__',
                        linewidth=linewidth,zorder=zord,s=size)

ms=60
# n-alkanes
molecule_plot(n_index,overlap_index_neat,n_legend,'$\it{n}$-Alkanes (C7-C18)',cf_se,
              cf_ed,color_use_pd,colors[9][0:4],colors[1][0:4],
              '^',1,10,ms)

# iso-alkanes
molecule_plot(iso_index,overlap_index_neat,iso_legend,'$\it{iso}$-Alkanes (C7-C17)',cf_se,
              cf_ed,color_use_pd,colors[9][0:4],colors[2][0:4],
              's',1,10,ms)

# Monocycloalkanes
molecule_plot(monocyclo_index,overlap_index_neat,monocyclo_legend,'Monocycloalkanes (C7-C16)',cf_se,
              cf_ed,color_use_pd,colors[9][0:4],colors[5][0:4],
              'o',1,10,ms) 

# Dicycloalkanes
molecule_plot(dicyclo_index,overlap_index_neat,dicyclo_legend,'Dicycloalkanes (C10-C14)',cf_se,
              cf_ed,color_use_pd,colors[9][0:4],colors[6][0:4],
              'p',1,10,ms) 

# Aromatics
molecule_plot(aro_index,overlap_index_neat,aro_legend,'Aromatics (C6-C15)',cf_se,
              cf_ed,color_use_pd,colors[9][0:4],colors[8][0:4],
              'd',1,10,ms) 

# Farnesane
ax1.scatter(cf_se[53],cf_ed[53],facecolor=np.array(color_use_pd.loc[[53]]),
            edgecolor=colors[4][0:4],marker='>',label='Farnesane',linewidth=1,
            s=ms)
#ax1.text(cf_se[23]+.2,cf_ed[23]-.2,'Farnesane',size='x-small',color='black',
#         horizontalalignment='right',verticalalignment='top')

# Pareto front
df_pareto=pd.read_csv('Overall_Pareto_Front_Summary_Unique_REV_10_thru_6274.csv')
ax1.plot(df_pareto['Specific Energy, MJ/kg'],df_pareto['Energy Density, MJ/L'],
         color='lime',label='Optimized Pareto Front',zorder=10,
         solid_capstyle='round',linewidth=3)

# Swell blend plot
#blend_legend=[]
#for ii in range(len(blend_swell)):
#    if len(blend_legend) == 0:
#        ax1.scatter(blend_se[ii],blend_ed[ii],color=np.array(color_use_pd.loc[[ii+len(df_conv)]]),
#                    marker='.',label='iso/Monocyclo/Dicyclo Blend',s=10,linewidth=1)
#        blend_legend.append(ii)
#    elif len(blend_legend) != 0:
#        ax1.scatter(blend_se[ii],blend_ed[ii],color=np.array(color_use_pd.loc[[ii+len(df_conv)]]),
#                    marker='.',label='__nolegend__',s=10,linewidth=1)        

pt2_plot()

# Swell legend
swell_patch=mpatches.Patch(color='darkgrey',label='Used for swell analysis')
no_swell_patch=mpatches.Patch(color=colors[9][0:4],label='Not used for swell analysis')
#no_swell_patch=mlines.Line2D([],[],color=colors[9][0:4],marker='p',linestyle='None',
#                          markersize=5,label='Molecule not Used for Swell Blending')
swell_legend=ax2.legend(handles=[swell_patch,no_swell_patch],
                        title='Molecule Fill Color:',title_fontsize='small',
                        prop={'size':'small'},loc='upper left',
                        handletextpad=.4)
swell_legend._legend_box.align='left'
swell_legend.get_frame().set_facecolor('none')
swell_legend.get_frame().set_edgecolor('none')

# Watermark
datafile=cbook.get_sample_data('Heat_Lab_bk.png',asfileobj=False)
im=image.imread(datafile)
im[:,:,-2]=0 # Deletes the white space
fig.figimage(im,1750,1900,alpha=.3,zorder=3)

# Save plot
plt.savefig(out_file,bbox_inches='tight',dpi=1000,transparent=False)
os.startfile(out_file,'open')