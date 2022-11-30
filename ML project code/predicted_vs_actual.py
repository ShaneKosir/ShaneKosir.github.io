import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib.ticker as plticker
import os
import pandas as pd
import matplotlib.image as image
import matplotlib.cbook as cbook
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
import matplotlib
import math
out_file_plot='Swell Actual vs Predicted 83 percent kfold 10.jpeg'
out_file_material='Swell Fit Metrics 83% kfold 10 by Material.xlsx'
out_file_mg='Swell Fit Metrics 83% kfold 10 by Molecular Group.xlsx'

#%% Definitions
def myround(x,base):
  return int(base*math.floor(float(x)/base))

def div_zero(num,den):
    return num/den if den else 0

def plot_molecules(axis,actual_data,predict_data,index,marker,edgecolor,
                   facecolor,label,size,zorder):
    axis.scatter(actual_data,predict_data,marker=marker,
            edgecolor=edgecolor,facecolor=facecolor[index],
            label=label,s=size,zorder=zorder)

def calc_mape(act,pred): # Mean absolute percent error
    percent_error_list=[]
    for ii in range(len(act)):
        if act[ii]!=0:
            percent_error_list.append(np.abs((act[ii]-pred[ii])/act[ii]))
    mape=np.mean(percent_error_list)*100
    return mape

def fit_metrics(train_df,test_df,material_column,actual_train_column,
                actual_test_column,train_column,cross_validate_column,
                test_column):
    actual_train=train_df[actual_train_column].values
    actual_test=test_df[actual_test_column].values
    train=train_df[train_column].values
    validate=train_df[cross_validate_column].values
    test=test_df[test_column].values
    train_matl=np.where(train_df[material_column].values==1)[0]
    test_matl=np.where(test_df[material_column].values==1)[0]
    mape_train=calc_mape(actual_train[train_matl],train[train_matl])
    mae_train=mean_absolute_error(actual_train[train_matl],train[train_matl])
    r2_train=r2_score(actual_train[train_matl],train[train_matl])
    mape_validate=calc_mape(actual_train[train_matl],validate[train_matl])
    mae_validate=mean_absolute_error(actual_train[train_matl],validate[train_matl])
    r2_validate=r2_score(actual_train[train_matl],validate[train_matl])
    mape_test=calc_mape(actual_test[test_matl],test[test_matl])
    mae_test=mean_absolute_error(actual_test[test_matl],test[test_matl])
    r2_test=r2_score(actual_test[test_matl],test[test_matl])
    return mape_train,mae_train,r2_train,mape_validate,mae_validate,r2_validate, \
           mape_test,mae_test,r2_test
           
def fit_metrics_mg(train_df,test_df,mg_column,actual_train_column,
                actual_test_column,train_column,cross_validate_column,
                test_column,mg):
    actual_train=train_df[actual_train_column].values
    actual_test=test_df[actual_test_column].values
    train=train_df[train_column].values
    validate=train_df[cross_validate_column].values
    test=test_df[test_column].values
    train_matl=np.where(train_df[mg_column].values==mg)[0]
    test_matl=np.where(test_df[mg_column].values==mg)[0]
    mape_train=calc_mape(actual_train[train_matl],train[train_matl])
    mae_train=mean_absolute_error(actual_train[train_matl],train[train_matl])
    r2_train=r2_score(actual_train[train_matl],train[train_matl])
    mape_validate=calc_mape(actual_train[train_matl],validate[train_matl])
    mae_validate=mean_absolute_error(actual_train[train_matl],validate[train_matl])
    r2_validate=r2_score(actual_train[train_matl],validate[train_matl])
    mape_test=calc_mape(actual_test[test_matl],test[test_matl])
    mae_test=mean_absolute_error(actual_test[test_matl],test[test_matl])
    r2_test=r2_score(actual_test[test_matl],test[test_matl])
    return mape_train,mae_train,r2_train,mape_validate,mae_validate,r2_validate, \
           mape_test,mae_test,r2_test

#%% Data
train_df=pd.read_csv('Train and validate.csv')
test_df=pd.read_csv('Test.csv')
actual_train=train_df['Actual'].values
actual_test=test_df['Actual'].values
train=train_df['Train'].values
validate=train_df['Cross Validate'].values
test=test_df['Predict'].values
mv_train=train_df['Molar Volume, mL/mol'].values
mg_train=train_df['Molecular Group'].values
mv_test=test_df['Molar Volume, mL/mol'].values
mg_test=test_df['Molecular Group'].values

#%% Combine
actual=np.concatenate([actual_train,actual_test])
mv_actual=np.concatenate([mv_train,mv_test])
mg_actual=np.concatenate([mg_train,mg_test])

#%% Calculate mape and r2 all
mape_train=calc_mape(actual_train,train)
mae_train=mean_absolute_error(actual_train,train)
r2_train=r2_score(actual_train,train)
mape_validate=calc_mape(actual_train,validate)
mae_validate=mean_absolute_error(actual_train,validate)
r2_validate=r2_score(actual_train,validate)
mape_test=calc_mape(actual_test,test)
mae_test=mean_absolute_error(actual_test,test)
r2_test=r2_score(actual_test,test)

#%% Calculate fit metrics for material groups
extracted_fit=fit_metrics(train_df,test_df,'N0602e','Actual',
                'Actual','Train','Cross Validate','Predict')
nitrile_fit=fit_metrics(train_df,test_df,'N0602','Actual',
                'Actual','Train','Cross Validate','Predict')
fluorosilicone_fit=fit_metrics(train_df,test_df,'L1120','Actual',
                'Actual','Train','Cross Validate','Predict')
low_temp_fluorocarbon_fit=fit_metrics(train_df,test_df,'V0835','Actual',
                'Actual','Train','Cross Validate','Predict')
lightweight_polysulfide_fit=fit_metrics(train_df,test_df,'PR-1776','Actual',
                'Actual','Train','Cross Validate','Predict')
polythioether_fit=fit_metrics(train_df,test_df,'PR-1828','Actual',
                'Actual','Train','Cross Validate','Predict')
epoxy2_fit=fit_metrics(train_df,test_df,'BMS 10-20','Actual',
                'Actual','Train','Cross Validate','Predict')
epoxy04_fit=fit_metrics(train_df,test_df,'BMS 10-123','Actual',
                'Actual','Train','Cross Validate','Predict')
nylon_fit=fit_metrics(train_df,test_df,'Nylon','Actual',
                'Actual','Train','Cross Validate','Predict')
kapton_fit=fit_metrics(train_df,test_df,'Kapton','Actual',
                'Actual','Train','Cross Validate','Predict')

#%% Export to excel
df_out=pd.DataFrame({'Extracted Nitrile Rubber':extracted_fit,
                     'Nitrile Rubber':nitrile_fit,
                     'Fluorosilicone':fluorosilicone_fit,
                     'Low Temp Fluorocarbon':low_temp_fluorocarbon_fit,
                     'Lightweight Polysulfide':lightweight_polysulfide_fit,
                     'Polythioether':polythioether_fit,
                     'Epoxy 0.2mm':epoxy2_fit,
                     'Epoxy 0.04mm':epoxy04_fit,
                     'Nylon':nylon_fit,
                     'Kapton':kapton_fit}).T
df_out.columns=['Train MAPE','Train MAE','Train R2','Validate MAPE',
                'Validate MAE','Validate R2','Test MAPE','Test MAE','Test R2']
df_out.to_excel(out_file_material)

#%% Calculate fit metrics for molecular groups
iso_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','iso-alkane')
monocyclo_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Monocycloalkane')
dicyclo_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Dicycloalkane')
aro_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Aromatic')
diaro_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Diaromatic')
cycloaro_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Cycloaromatic')
diamond_fit=fit_metrics_mg(train_df,test_df,'Molecular Group','Actual',
                'Actual','Train','Cross Validate','Predict','Diamondoid')

#%% Export to excel
df_out=pd.DataFrame({'iso-alkane':iso_fit,
                     'Monocycloalkane':monocyclo_fit,
                     'Dicycloalkane':dicyclo_fit,
                     'Aromatic':aro_fit,
                     'Diaromatic':diaro_fit,
                     'Cycloaromatic':cycloaro_fit,
                     'diamondoid':diamond_fit}).T
df_out.columns=['Train MAPE','Train MAE','Train R2','Validate MAPE',
                'Validate MAE','Validate R2','Test MAPE','Test MAE','Test R2']
df_out.to_excel(out_file_mg)

#%% Plot
# Setup
fig,ax1=plt.subplots()
colors=cm.tab10(np.linspace(0,1,10))

# Adding colors based on colorbar to plot
xxx=mv_actual # molar volumes to normalize
if min(xxx)==0:
    min_norm=.1e-9
else:
    min_norm=min(xxx)
max_norm=max(xxx)
#norm=matplotlib.colors.Normalize(vmin=min_norm,vmax=max_norm) # Normalize the colors from [0,1] based on a min and max value
norm=matplotlib.colors.LogNorm(vmin=min_norm,vmax=max_norm) # Log10 normalize the colors from [0,1] based on a min and max value
cmap_name='Greys' # _r to reverse
cmap=matplotlib.cm.get_cmap(cmap_name)
color_use=cmap(norm(xxx))

# Add colorbar to plot
ax2=fig.add_axes([.6,.19,.25,0.025]) # L,B,W,H
num_ticks=4
ticks=np.linspace(0, 1, num_ticks)
labels=np.linspace(min_norm,max_norm,num_ticks)
labels=[myround(x,1) for x in labels]
cbar=matplotlib.colorbar.ColorbarBase(ax2,cmap=cmap_name,
                                      orientation='horizontal')
cbar.set_ticks(ticks[::1])
cbar.set_ticklabels(labels[::1])
cbar.ax.tick_params(labelsize='small')
cbar.set_label('Molar Volume, mL/mol',rotation=00,
               labelpad=-31,size='small')
ax2.yaxis.set_label_position('right')

# Plot data
s=50
train_index=list(range(len(actual_train)))
test_index=list(range(len(actual_train),len(actual_train)+len(actual_test)))
plot_molecules(ax1,actual_train,train,train_index,'.',colors[0][0:4],
               color_use,'Train, n=%i' %(len(actual_train)),s,1)
plot_molecules(ax1,actual_train,validate,train_index,'.',colors[1][0:4],
               color_use,'Validate, n=%i' %(len(actual_train)),s,1)
plot_molecules(ax1,actual_test,test,test_index,'.',colors[2][0:4],
               color_use,'Test, n=%i' %(len(actual_test)),s,1)

# Label axes
ax1.set_xlabel('Extrapolated Neat Volume Swell, % v/v',size='x-large')
ax1.set_ylabel('Predicted Volume Swell, % v/v',size='x-large')

# Axis limits
min_val=min(np.concatenate([actual,train,validate,test]))-5
max_val=max(np.concatenate([actual,train,validate,test]))+5
x_lim=y_lim=[min_val,max_val]
ax1.set_xlim(x_lim)
ax1.set_ylim(y_lim)

# 1:1 line
x=np.linspace(x_lim[0],x_lim[1],100)
y=np.linspace(y_lim[0],y_lim[1],100)
ax1.plot(x,y,color='black',linestyle='--',linewidth=.75,zorder=0)

# MAPE lines
#y_above=[y[ii]+y[ii]*mape/100 for ii in range(len(y))]
#y_below=[y[ii]-y[ii]*mape/100 for ii in range(len(y))]
#ax1.plot(x,y_above,color='blue',linestyle='--',linewidth=.5,label='MAPE',
#         zorder=0)
#ax1.plot(x,y_below,color='blue',linestyle='--',linewidth=.5,zorder=0)

# Major ticks
ax1.locator_params(axis='x',nbins=7)
ax1.locator_params(axis='y',nbins=7)

# Minor Ticks
ax1.minorticks_on()
ax1.tick_params(axis='both', which='minor', direction='out')
loc=plticker.MultipleLocator(base=25) # x
ax1.xaxis.set_minor_locator(loc)
loc=plticker.MultipleLocator(base=25) # y
ax1.yaxis.set_minor_locator(loc)

# Label
ax1.text(60,185,
         'Train R-squared='+str(round(r2_train,2)) \
         + '\nTrain MAE='+str(round(mae_train,1))+'% v/v',
         fontsize='small',horizontalalignment='left',color=colors[0][0:4],
         weight='bold')
ax1.text(60,163,
         'Validate R-squared='+str(round(r2_validate,2)) \
         + '\nValidate MAE='+str(round(mae_validate,1))+'% v/v',
         fontsize='small',horizontalalignment='left',color=colors[1][0:4],
         weight='bold')
ax1.text(60,141,
         'Test R-squared='+str(round(r2_test,2)) \
         + '\nTest MAE='+str(round(mae_test,1))+'% v/v',
         fontsize='small',horizontalalignment='left',color=colors[2][0:4],
         weight='bold')

# Legend
legend=ax1.legend(loc='upper left',prop={'size':'small'},ncol=1,
                  handletextpad=-.2)

legend.legendHandles[0].set_facecolor('lightgrey')
legend.legendHandles[1].set_facecolor('lightgrey')
legend.legendHandles[2].set_facecolor('lightgrey')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')

# Outliers
# 1
x_head=actual_train[16]
y_head=validate[16]
x_text=actual_train[16]-15
y_text=validate[16]
ax1.annotate("",xy=(x_head-2.5,y_head),xytext=(x_text,y_text),
             arrowprops=dict(headlength=10,headwidth=10,width=4,
                             facecolor='white',edgecolor='black',
                             linewidth=1))
ax1.text(x_text-1,y_text+4,train_df['Molecule'][16]+'\n(Extracted Nitrile\nRubber)',
         size='small',horizontalalignment='right',verticalalignment='top',
         color='black')
# 2
x_head=actual_train[111]
y_head=validate[111]
x_text=actual_train[111]-10
y_text=validate[111]-17.5
ax1.annotate("",xy=(x_head-1.5,y_head-3.5),xytext=(x_text,y_text),
             arrowprops=dict(headlength=10,headwidth=10,width=4,
                             facecolor='white',edgecolor='black',
                             linewidth=1))
ax1.text(x_text,y_text-2.5,train_df['Molecule'][111]+'\n(Polythioether)',
         size='small',horizontalalignment='center',verticalalignment='top',
         color='black')

# 3
x_head=actual_train[54]
y_head=validate[54]
x_text=actual_train[54]-10
y_text=validate[54]+17.5
ax1.annotate("",xy=(x_head-2,y_head+3),xytext=(x_text,y_text),
             arrowprops=dict(headlength=10,headwidth=10,width=4,
                             facecolor='white',edgecolor='black',
                             linewidth=1))
ax1.text(x_text,y_text+1.5,train_df['Molecule'][54]+'\n(Fluorosilicone)',
         size='small',horizontalalignment='center',verticalalignment='bottom',
         color='black')

# Watermark
datafile=cbook.get_sample_data('Heat_Lab_bk.png',asfileobj=False)
im=image.imread(datafile)
im[:,:,-2]=0 # Deletes the white space
fig.figimage(im,1350,700,alpha=.3,zorder=3)

plt.savefig(out_file_plot,bbox_inches='tight',dpi=1000,transparent=False)
os.startfile(out_file_plot,'open')