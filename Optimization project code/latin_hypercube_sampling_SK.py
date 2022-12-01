import numpy as np
from scipy.special import erfinv
import pandas as pd
import copy
in_file='conv_molecules.xlsx'
# Out file for compiled lhs spreadsheet
out_file_all=in_file.replace('.xlsx','')+'_uncert_all.xlsx'
# Set number of lhs samples per molecule/property
N=100
# Set properties of interest (uncert columns pulled as well)
# Uncert columns need to begin with "Uncert"
prop_cols=['Specific Energy, MJ/kg','Energy Density, MJ/L',
           'Density @ T1, g/cc','Viscosity at T2, mm^2/s',
           'Melting Point, C','Flash Point, C','Boiling Point, C',
           'DCN']
# Add columns without uncertainties to initialization spreadsheets
col_add=['Molecular Group','Formula','H/C','MW, g/mol',
         'Aromatic Volume Percent','Decalin Percent']

### Outputs a N dimensional vector for ONE property for ONE fuel
def lhs_rand_sample_out(mu,sigma,INV_ERF):
    ## mu is the average value for the property
    ## sigma is one standard deviation of the property
    ## INV_ERF is generated in a seperage function
    N= len(INV_ERF)
    out = np.zeros(N-1)
    sqrt2 = 2**.5
    x = sigma*sqrt2*INV_ERF+np.ones(N)*mu
    # cdf_intervals = np.linspace(0,1,num=N+1,endpoint=True)
    for i in range(N-1):
        out[i]=np.random.random()*(x[i+1]-x[i])+x[i]
        
#    print(out,len(out),np.average(out),np.std(out))
    np.random.shuffle(out)
#    print(len(out),np.average(out),np.std(out),out)
    return out

### Calculate the CDF intervals to sample from 
def inv_erf(N):
    ## N is the number of samples
    CDF = np.linspace(1e-9,0.999999999,num=N+1)
    INV_ERF = erfinv(CDF*2-1)
    return INV_ERF

#def get_properties():
#    N = 100 ## Equal to the total number of 'spreadsheet' samples or different solutions
#    INV_ERF = inv_erf(N) 
#    # insert the taking of mu and sigma from other fuels and properties
#    # dummy data
#    mu = 1
#    sigma = .15
#    lhs_rand_sample_out(mu,sigma,INV_ERF)
#get_properties()

### Outputs a N dimensional vector for X properties and Y molecules
def get_props_mult(N,in_file,prop_cols,out_file_all):
    # Set number of lhs samples
    N=N
    inverf=inv_erf(N)
    global df
    df=pd.read_excel(in_file)
    # Select property columns
    prop_cols=prop_cols
    uncert_cols=['Uncert '+i for i in prop_cols]
    df_props=df[prop_cols]
    if df_props.isnull().sum().sum() != 0:
        print(df_props.isnull().sum().sum(),'PROPERTY VALUES MISSING!!!')
    df_uncert=df[uncert_cols]
    if df_uncert.isnull().sum().sum() != 0:
        print(df_uncert.isnull().sum().sum(),'UNCERTAINTY VALUES MISSING!!!')
    # Perform lhs for each [molecule,property]
    global lhs_vals
    lhs_vals=pd.DataFrame(index=df['Molecule'],columns=prop_cols)
    for ii in range(len(prop_cols)):
        for jj in range(len(df_props)):
            lhs_vals.iloc[jj,ii]=list(lhs_rand_sample_out(df_props.iloc[jj,ii], \
                         df_uncert.iloc[jj,ii],inverf)) # iloc[row,col]
    # Add columns without uncertainty
    lhs_vals_all=copy.deepcopy(lhs_vals)
    for ii in reversed(range(len(col_add))):
        lhs_vals_all.insert(0,col_add[ii],df[col_add[ii]].values)    
    lhs_vals_all.to_excel(out_file_all)
get_props_mult(N,in_file,prop_cols,out_file_all)

def zero_filt(df,props):
    for ii in range(len(props)):
        col_new=np.where(df[props[ii]]<=0,1e-9,df[props[ii]])
        df[props[ii]]=col_new
    return df

### Propogate N dimensional vectors to N initialization spreadsheets
def init_spreadsheets(N,lhs_vals,col_add):
    for ii in range(N):
        global ss, ss_filt
        ss=lhs_vals.applymap(lambda x: x[ii])
        # Add columns without uncertainties
        for jj in reversed(range(len(col_add))):
            ss.insert(0,col_add[jj],df[col_add[jj]].values)
            # Replace zero values
        ss_filt=zero_filt(ss,['Specific Energy, MJ/kg','Energy Density, MJ/L',
                              'Density @ T1, g/cc','Viscosity at T2, mm^2/s',
                              'DCN'])
        ss_filt.to_excel(in_file.replace('.xlsx','')+'_uncert_'+str(ii)+'.xlsx')
init_spreadsheets(N,lhs_vals,col_add)