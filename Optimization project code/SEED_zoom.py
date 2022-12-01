import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import cm
import matplotlib
import numpy as np
import matplotlib.ticker as plticker
import math
from matplotlib.patches import Polygon
import seaborn as sns
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy import interpolate
import matplotlib.cbook as cbook
import matplotlib.image as image
out_file='SEED_zoom.jpg'

#%% Definitions
def linear_interpolation(x, y, functions):
    f = interpolate.interp1d(x, y, bounds_error = False, fill_value = 'NaN')
    functions.append(f)

# Reading Data and Creating Linear Interpolatino Functions
def reading_data(subdirname):
    loc=os.getcwd()+'\\'+subdirname
    lst = [i for i in os.listdir(loc) if i.startswith('Overall_Pareto_Front_Summary_Unique_SS')]
    overall_se = []
    overall_ed = []
    fed=[]
    for i in lst:
        data = pd.read_csv(loc+'\\'+i)
        se = data['Specific Energy, MJ/kg']
        overall_se.append(se)
        ed = data['Energy Density, MJ/L']
        overall_ed.append(ed)
        linear_interpolation(se, ed, fed)
    return overall_se,overall_ed,fed
        
# Converting Arrays   
def converting_arrays(array):
    array = np.asarray(array)
    array = np.ndarray.tolist(array)
    return array

#%% Read data
se_aro,ed_aro,fed_aro=reading_data('aro')
se_no_aro,ed_no_aro,fed_no_aro=reading_data('no_aro')

#%% test
#test=[new_overall_ed[ii][0] for ii in range(len(new_overall_ed))]

#%% Binning SE new (aro)
se_min_list=[se_aro[ii].min() for ii in range(len(se_aro))]
se_min=min(se_min_list)
se_max_list=[se_aro[ii].max() for ii in range(len(ed_aro))]
se_max=max(se_max_list)
bins_aro=np.linspace(se_min,se_max,100)

#%% Binning SE new (no aro)
se_min_list=[se_no_aro[ii].min() for ii in range(len(se_no_aro))]
se_min=min(se_min_list)
se_max_list=[se_no_aro[ii].max() for ii in range(len(ed_no_aro))]
se_max=max(se_max_list)
bins_no_aro=np.linspace(se_min,se_max,100)

#%% Evaluating Functions of binned se
def using_functions(functions,bins):
    new_overall_y=[]
    for i in functions:
        new_y=[]
        for j in bins:
            y_value = i(j)
            new_y.append(y_value)
        new_y = converting_arrays(new_y)
        new_overall_y.append(new_y)
    return new_overall_y
        
new_overall_ed_aro=using_functions(fed_aro,bins_aro)
new_overall_ed_no_aro=using_functions(fed_no_aro,bins_no_aro)

#%% Calculating the median
#def median(new_overall_y):
#    mu = np.nanmedian(new_overall_y, axis = 0)
#    mu = converting_arrays(mu)
#    return mu

#mu_ed = median(new_overall_ed)

#%% Calculating Sigma        
#def sigma_calc(new_overall_y):
#    sigma = np.nanstd(new_overall_y, axis=0)
#    return sigma

#sigma_ed = sigma_calc(new_overall_ed)

#%% Calculate IQR and outliers
def IQR(new_overall_y):
    two_sigma_low,sigma_low,med,sigma_high,two_sigma_high= \
    np.nanpercentile(new_overall_y,[2.275,15.865,50,84.135,97.725],axis=0)
    return two_sigma_low,sigma_low,med,sigma_high,two_sigma_high

ed_stats_aro=IQR(new_overall_ed_aro)
ed_stats_no_aro=IQR(new_overall_ed_no_aro)

#%% Plot stats
def plot_stats(ax,bins,stats,color,capstyle,label,alpha_dark,alpha_light,
               lw_pf,lw_border):
    ax.plot(bins,stats[2],color=color,solid_capstyle=capstyle,lw=lw_pf,
            label=label)
    ax.plot(bins,stats[0],color=color,solid_capstyle=capstyle,lw=lw_border)
    ax.plot(bins,stats[1],color=color,solid_capstyle=capstyle,lw=lw_border)
    ax.plot(bins,stats[3],color=color,solid_capstyle=capstyle,lw=lw_border)
    ax.plot(bins,stats[4],color=color,solid_capstyle=capstyle,lw=lw_border)
    ax.fill_between(bins,stats[2],stats[3],alpha=alpha_dark,facecolor=color)
    ax.fill_between(bins,stats[2],stats[1],alpha=alpha_dark,facecolor=color)
    ax.fill_between(bins,stats[3],stats[4],alpha=alpha_light,facecolor=color)
    ax.fill_between(bins,stats[1],stats[0],alpha=alpha_light,facecolor=color)

#%%
def myround(x,base):
  return int(base*math.ceil(float(x)/base))

def low_int_hpf (x, ED):
    """
    Functions for FSOLVE.
    Compares PQIS Outer Region to to
    lower floor (SE*0.775) = ED
    """
    return ED - lower_lf(x)
        
def hpf_shading_r2():
    global A2_SE, A2_ED
    A2_SE = 43.06
    A2_ED = 0.803*43.06
    
    # Starting bounds from upper right region working around clockwise
    st_list = [[50,0.84*50],[50, 0.775*50]]
    
    # PT right by the lower bound and intercept of Jet-A PQIS Scaled Data
    x_pt_RIGHT = fsolve(low_int_hpf, 43.5, args=(A2_ED))[0]
    y_pt_RIGHT = x_pt_RIGHT*0.775
    st_list = [*st_list, [x_pt_RIGHT, y_pt_RIGHT]]
    
    # POINT of intercept with outliers of PQIS Database.
    global fn_pqis_interp, fn_pqis_interp_rev
    fn_pqis_interp = interp1d(pareto_x, pareto_y, fill_value='extrapolate')
    fn_pqis_interp_rev = interp1d(pareto_y, pareto_x, fill_value='extrapolate')    
    int_pt_low_SE = fn_pqis_interp_rev(y_pt_RIGHT)
    int_pt_low_ED = y_pt_RIGHT    
    st_list = [*st_list, [int_pt_low_SE, int_pt_low_ED]]
    
    # POINTS on PQIS OUTLIER 
    pareto_x_2 = [i for i in pareto_x if i > A2_SE and pareto_y[pareto_x.index(i)]> A2_ED]
    pareto_y_2 = [pareto_y[pareto_x.index(i)] for i in pareto_x_2]
    pf = []
    for i in reversed(range(len(pareto_x_2))):
        pf = [*pf, [pareto_x_2[i], pareto_y_2[i]]]
    st_list = [*st_list, *pf]
    
    # POINT where A2 SE intercepts on PQIS Outliers
    st_list = [*st_list, [A2_SE, fn_pqis_interp(A2_SE)]]
    
    # Point where PQIS SE Line, and then to the upper limit 
    st_list = [*st_list, [A2_SE, A2_SE*0.84]]
    st_list = [*st_list, [50., 50.*0.84]]
    
    ### HPF Range
    ax1.add_patch(Polygon(st_list,closed=True,fill=True,lw=.3,facecolor='blue',zorder=0,alpha=.3))
    ax1.text(44.1,35.7,'HPF Region', size='medium', color='blue',horizontalalignment='center',
        verticalalignment='center',zorder=10)

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
    global jet_se, jet_ed
    jet_se=[43.26950331,43.20125921,43.14795266,43.25658127] 
    jet_ed=[34.41616224,34.9072202,35.05963253,34.65312478]    
    scaling_conv_fuels()

    # Plot setup
    global fig, ax1, ax2, ax3
    fig,ax1=plt.subplots()    
    ax2=ax1.twinx()
    ax3=ax1.twiny()
    
    # Set tick number
    ax1.locator_params(axis='y',nbins=5)
    ax2.locator_params(axis='y',nbins=4)
    ax3.locator_params(axis='x',nbins=5)
    
    # Minor Ticks
    plt.minorticks_on()
    plt.tick_params(axis='both', which='minor', direction='out')
    ax1.minorticks_on()
    ax1.tick_params(axis='x', which='minor', direction='out')
    ax2.minorticks_on()
    ax2.tick_params(axis='y', which='minor', direction='out')
    loc=plticker.MultipleLocator(base=.1) # x
    ax1.xaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=.25) # y
    ax1.yaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=1) # y2
    ax2.yaxis.set_minor_locator(loc)
    loc=plticker.MultipleLocator(base=.25) # x2
    ax3.xaxis.set_minor_locator(loc) 

def plot_conv_pareto_front_PQIS():
    global pareto_x, pareto_y
    # Plot Conventional Pareto Front
    pareto_x=[43.1,43.3,43.54,43.6,43.6,43.6] # 43.1,43.6
    pareto_y=[36.204,35.50599999999999,35.2674,34.88,34.444,34.008] # 35.77,33.79
    pareto_dens = [pareto_y[i]/pareto_x[i] for i in range(len(pareto_x)) if len(pareto_x) == len(pareto_y)]
    
    if "res_SE" in globals():
        pareto_x = [i-res_SE for i in pareto_x]
        pareto_y = [pareto_x[i]*pareto_dens[i] for i in range(len(pareto_x))]    
    ax1.plot(pareto_x,pareto_y,color='blue',solid_capstyle='round',
             label='Best-case\nConventional Fuels')

def pt2_plot():
    data_sc =np.load('jittered_data.npy')
    
    # Rescaling data to keep density constant from database
    # But to switch the density and recalculate the energy density.
    if "res_SE" in globals():
        for i in range(len(data_sc)):
            dens = data_sc[i,1]/data_sc[i,0]
            data_sc[i,0] = data_sc[i,0] - res_SE
            data_sc[i,1] = data_sc[i,0] * dens
    plot_conv_pareto_front_PQIS()
    
    # Plot Gradient
    sns.kdeplot(data_sc[:,0],data2=data_sc[:,1],ax=ax1,zorder=0)

    # Set axis limits  
    x_lim=[43.35,44.3]
    y_lim=[34.25,36.5]
    ax1.set_xlim(x_lim)
    ax1.set_ylim(y_lim)
    
    # REFerences for lines
    # Commented out the previous ones, because scaling 
    # code above should fix this.
    if "REF_SE" and "REF_ED" not in globals():
        global REF_SE, REF_ED
        REF_SE = 43.20125921
        REF_ED = 34.9072202
    
    # Reference with Jet A, if not comment that out
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
    ax1.set_xlabel('Specific Energy, MJ/kg',size='large')
    ax1.set_ylabel('Energy Density, MJ/L',size='large')
    ax2.set_ylabel('Energy Density, % Diff.',size='large')
    ax3.set_xlabel('Specific Energy, % Diff.',size='large')

    # Specific energy spec limit
    xrange_se=np.arange(42.8,x_lim[1]+1,.5)
    se_spec_limit_lower=ax1.plot((42.8,42.8),y_lim,color='red',zorder=0,linewidth=.5)

    ## Spec. Limits
#    ax1.text(42.75,37.7,'Jet A Specific\nEnergy Limit',size='x-small',
#             color='red',ha='right')

    # Density limit
    xrange=np.arange(42.8,x_lim[1]+1,.5)
    yrange_low=[xrange[ii]*.775 for ii in range(len(xrange))]
    yrange_high=[xrange[ii]*.84 for ii in range(len(xrange))]
    density_lower=ax1.plot(xrange,yrange_low,color='blue',zorder=0,linewidth=.5)
    density_upper=ax1.plot(xrange,yrange_high,color='blue',zorder=0,linewidth=.5)
    ax1.fill_between(xrange,yrange_low,yrange_high,facecolor='blue',alpha=.1,zorder=0)
 
    ## HPF Shading
    if HPF_LOG == 1:
        hpf_shading_r2() 

def get_CFIT_DENS_RANGE ():
    """Curve fit and interpolations for
    upper and lower regions of Jet Fuel"""
    x_FIT = [42.8, 46.]
    y_FIT_LOW = [i*0.775 for i in x_FIT]
    y_FIT_UP = [i*0.84 for i in x_FIT]
    
    global lower_lf, upper_lf
    lower_lf = interp1d(x_FIT, y_FIT_LOW, fill_value='extrapolate')
    upper_lf = interp1d(x_FIT, y_FIT_UP, fill_value='extrapolate')

def plot_NJFCP_fuels():
    colors=cm.tab10(np.linspace(0,1,10))
#    ax1.scatter(42.88,0.827*42.88,facecolor=colors[3][0:4],edgecolor='black',
#                    marker='D',label='JP-5',linewidth=.5,s=35,zorder=10)    
#    ax1.scatter(43.06,0.803*43.06,facecolor=colors[3][0:4],edgecolor='black',
#                    marker='h',label='Jet A',linewidth=.5,s=60,zorder=10)
#    ax1.scatter(43.24,0.780*43.24,facecolor=colors[3][0:4],edgecolor='black',
#                    marker='P',label='JP-8',linewidth=.5,s=70,zorder=10)

pt1_plot()
get_CFIT_DENS_RANGE()
plot_NJFCP_fuels()

global HPF_LOG
HPF_LOG = 1

### Conventional molecules
# Import conventional fuel dataframes
df_conv=pd.read_excel('conv_molecules.xlsx')

# Only include molecules with property values
df_conv=df_conv.dropna(subset=['Specific Energy, MJ/kg','Density @ T1, g/cc'])

# Reset indices
df_conv=df_conv.reset_index(drop=True)

# Reset indices
df_conv=df_conv.reset_index(drop=True)

# Convert pandas dataframe to numpy arrays for plotting
cf_se=df_conv['Specific Energy, MJ/kg'].tolist()
number_of_fuels=len(cf_se)
cf_rho=df_conv['Density @ T1, g/cc'].tolist()
cf_ed=[cf_se[ii]*cf_rho[ii] for ii in range(number_of_fuels)]
cf_mw=df_conv['MW, g/mol'].values
colors=cm.tab10(np.linspace(0,1,10))

# Legend indices
n_legend=[]
iso_legend=[]
monocyclo_legend=[]
dicyclo_legend=[]
aro_legend=[]

# Normalize colors
xxx=cf_mw # molecular weights to normalize
min_norm=min(xxx)
max_norm=max(xxx)
norm=matplotlib.colors.Normalize(vmin=min_norm,vmax=max_norm) # Normalize the colors from [0,1] based on a min and max value
cmap_name='Greys' # _r to reverse
cmap=matplotlib.cm.get_cmap(cmap_name)
color_use=cmap(norm(xxx))

# Add colorbar to plot
ax4=fig.add_axes([.67,.8,.2,0.02]) # L,B,W,H
num_ticks=4
ticks=np.linspace(0,1,num_ticks)
labels=np.linspace(min_norm,max_norm,num_ticks)
labels=[myround(x,1) for x in labels]
cbar=matplotlib.colorbar.ColorbarBase(ax4,cmap=cmap_name,orientation='horizontal')
cbar.set_ticks(ticks[::1])
cbar.set_ticklabels(labels[::1])
cbar.ax.tick_params(labelsize='x-small')
cbar.set_label('Molecular Weight, g/mol',rotation=00,labelpad=-27.5,size='x-small')
ax4.yaxis.set_label_position('right')
ax4.tick_params(size=0,which='minor')

# Plot molecules
ms=50
ax1.scatter(cf_se[0:11],cf_ed[0:11],facecolor=color_use[0:11],
            edgecolor=colors[1],marker='^',label='$\it{n}$'+'-Alkane',
            linewidth=1,s=ms)
ax1.scatter(cf_se[11:23],cf_ed[11:23],facecolor=color_use[11:23],
            edgecolor=colors[2],marker='s',label='$\it{iso}$'+'-Alkane',
            linewidth=1,s=ms)
ax1.scatter(cf_se[23:33],cf_ed[23:33],facecolor=color_use[23:33],
            edgecolor=colors[5],marker='o',label='Monocycloalkane',
            linewidth=1,s=ms)
#ax1.scatter(cf_se[33:37],cf_ed[33:37],facecolor=color_use[33:37],
#            edgecolor=colors[6],marker='p',label='Dicycloalkane',
#            linewidth=1,s=ms)
#ax1.scatter(cf_se[37:53],cf_ed[37:53],facecolor=color_use[37:53],
#            edgecolor=colors[8],marker='d',label='Aromatic',
#            linewidth=1,s=ms)
#ax1.scatter(cf_se[53],cf_ed[53],facecolor=color_use[53],
#            edgecolor=colors[4],marker='D',label='Diaromatic',
#            linewidth=1,s=30)

pt2_plot()

# Plot pareto front with sigma shading
plot_stats(ax1,bins_aro,ed_stats_aro,colors[1],'round',
           'Pareto Front\n(Aro Const.)',.5,.25,1.5,.5)
plot_stats(ax1,bins_no_aro,ed_stats_no_aro,colors[2],'round',
           'Pareto Front\n(no Aro Const.)',.5,.25,1.5,.5)
#ax1.plot(bins,mu_ed,color=colors[1],solid_capstyle='round',
#         label='Optimized Pareto Front')
#ax1.fill_between(bins,mu_ed+sigma_ed,mu_ed-sigma_ed,alpha=.65,
#                 facecolor=colors[1])
#ax1.fill_between(bins,mu_ed+sigma_ed,mu_ed+2*sigma_ed,alpha=.4,
#                 facecolor=colors[1])
#ax1.fill_between(bins,mu_ed-sigma_ed,mu_ed-2*sigma_ed,alpha=.4,
#                 facecolor=colors[1])

# Legend
legend=ax1.legend(loc='lower left',bbox_to_anchor=(.12,0),prop={'size':6},
                    ncol=1,labelspacing=.9,columnspacing=.5,
                    facecolor='none',edgecolor='none')
legend.legendHandles[3].set_facecolor('darkgrey')
legend.legendHandles[4].set_facecolor('darkgrey')
legend.legendHandles[5].set_facecolor('darkgrey')
#legend.legendHandles[9].set_facecolor('darkgrey')
#legend.legendHandles[10].set_facecolor('darkgrey')
#legend.legendHandles[11].set_facecolor('darkgrey')

# Watermark
datafile=cbook.get_sample_data('Heat_Lab_bk.png',asfileobj=False)
im=image.imread(datafile)
im[:,:,-2]=0 # Deletes the white space
fig.figimage(im,2825,1500,alpha=.3,zorder=3)

# Save plot
plt.savefig(out_file,bbox_inches='tight',dpi=1000,transparent=False)
os.startfile(out_file,'open')