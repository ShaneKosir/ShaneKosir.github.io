"""This code assumes that aliphatic double and triple bonds have been parsed"""
from rdkit import Chem
import pandas as pd
import numpy as np
from itertools import groupby
from operator import itemgetter
in_file='Predict Molecules.xlsx'
out_file='Swell Predict Functional Groups.xlsx'

#%% Definitions
def find(s,ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def div_zero(num,den):
    return num/den if den else np.nan

def group_parse(functional_group,molecule):
    smarts_pattern=functional_group
    smarts_retrieve=Chem.MolFromSmarts(smarts_pattern)
    match=Chem.Mol.GetSubstructMatches(molecule,smarts_retrieve,uniquify=True)
    num=len(match)
    return num

def avg_len(lst):
    lengths = [len(i) for i in lst]
    return 0 if len(lengths) == 0 else (float(sum(lengths)) / len(lengths))

def len_ch2(molecule):
    smarts_retrieve=Chem.MolFromSmarts('[CH2;!R]')
    match=Chem.Mol.GetSubstructMatches(molecule,smarts_retrieve,uniquify=True)
    match=[ii[0] for ii in match]
    match_chunks=[]
    for k, g in groupby(enumerate(match), lambda x: x[0]-x[1]):
        match_chunks.append(list(map(itemgetter(1), g)))
    average_length=avg_len(match_chunks)
    return average_length

def ring_info(molecule):
    atom_rings=molecule.GetRingInfo().AtomRings()
    ring_size=[]
    is_aro=[]
    for ii in range(len(atom_rings)):
        ring_size.append(len(atom_rings[ii]))
        is_aro.append(mol.GetAtomWithIdx(atom_rings[ii][0]).GetIsAromatic())
    num_aro=sum(is_aro)
    num_cyclo=len(atom_rings)-num_aro
    cyclo_index=[ii for ii, x in enumerate(is_aro) if not x]
    cyclo_avg_size=np.average([ring_size[ii] for ii in cyclo_index])
#    cyclo_avg_size=zero_avg(ring_size,cyclo_index)
    return num_aro,num_cyclo,cyclo_avg_size

def sub_info(molecule):
    sub_positions=molecule.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
    atom_rings=mol.GetRingInfo().AtomRings()
    is_aro_sub=[]
    is_aro_ring=[]
    for ii in range(len(sub_positions)):
        is_aro_sub.append(mol.GetAtomWithIdx(sub_positions[ii][1]).GetIsAromatic())
    num_aro_sub=sum(is_aro_sub)
    num_cyclo_sub=len(sub_positions)-num_aro_sub
    for ii in range(len(atom_rings)):
        is_aro_ring.append(mol.GetAtomWithIdx(atom_rings[ii][0]).GetIsAromatic())
    num_aro_ring=sum(is_aro_ring)
    num_cyclo_ring=len(atom_rings)-num_aro_ring
    avg_sub_per_cyclo=div_zero(num_cyclo_sub,num_cyclo_ring)
    avg_sub_per_aro=div_zero(num_aro_sub,num_aro_ring)
    return avg_sub_per_aro,avg_sub_per_cyclo

def chain_info(molecule):
    sub_positions=molecule.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
    is_aro_sub=[]
    for ii in range(len(sub_positions)):
        is_aro_sub.append(mol.GetAtomWithIdx(sub_positions[ii][1]).GetIsAromatic())
    aro_sub_index=[i for i, x in enumerate(is_aro_sub) if x]
    cyclo_sub_index=[i for i, x in enumerate(is_aro_sub) if not x]
    atom_rings=mol.GetRingInfo().AtomRings()
    is_aro_ring=[]
    for ii in range(len(atom_rings)):
        is_aro_ring.append(mol.GetAtomWithIdx(atom_rings[ii][0]).GetIsAromatic())
    all_atoms=[]
    for atom in mol.GetAtoms():
        all_atoms.append(atom.GetIdx())
    sub_indices_chain=[x[0] for x in sub_positions] # position of first non-ring carbon
    ring_atoms=mol.GetRingInfo().AtomRings()
    ring_atoms_flat=list(sum(ring_atoms,()))
    chain_atoms=[x for x in all_atoms if x not in ring_atoms_flat]
    if len(sub_positions) == 0:
        avg_len_aro=np.nan
        avg_len_cyclo=np.nan
    if len(sub_positions) == 1:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) == 2:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        len_chains[1]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[1]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) == 3:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        len_chains[1]=len([ii for ii in chain_atoms if ii > sub_indices_chain[0] and ii < sub_indices_chain[-1]])
        len_chains[2]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[-1]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) == 4:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        len_chains[1]=len([ii for ii in chain_atoms if ii > sub_indices_chain[0] and ii < sub_indices_chain[2]])
        len_chains[2]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[2] and ii < sub_indices_chain[-1]])
        len_chains[3]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[-1]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) == 5:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        len_chains[1]=len([ii for ii in chain_atoms if ii > sub_indices_chain[0] and ii < sub_indices_chain[2]])
        len_chains[2]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[2] and ii < sub_indices_chain[3]])
        len_chains[3]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[3] and ii < sub_indices_chain[-1]])
        len_chains[4]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[-1]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) == 6:
        len_chains=[[] for ii in range(len(sub_positions))]
        len_chains[0]=len([ii for ii in chain_atoms if ii <= sub_indices_chain[0]])
        len_chains[1]=len([ii for ii in chain_atoms if ii > sub_indices_chain[0] and ii < sub_indices_chain[2]])
        len_chains[2]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[2] and ii < sub_indices_chain[3]])
        len_chains[3]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[3] and ii < sub_indices_chain[4]])
        len_chains[4]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[4] and ii < sub_indices_chain[-1]])
        len_chains[5]=len([ii for ii in chain_atoms if ii >= sub_indices_chain[-1]])
        avg_len_aro=np.average([len_chains[ii] for ii in aro_sub_index])
        avg_len_cyclo=np.average([len_chains[ii] for ii in cyclo_sub_index])
    if len(sub_positions) > 6:
        print('Too Many Alkyl Chains!!! GO HOME')
    return avg_len_aro,avg_len_cyclo

#%% Load data
original=pd.read_excel(in_file)
smiles=original['Canonical Smiles']

#%% Dataframe to append
funct_groups=pd.DataFrame({'(linear) -CH3':['na']*len(smiles),
                     '(linear) -CH2-':['na']*len(smiles),
                     '(linear) >CH-':['na']*len(smiles),
                     '(linear) >C<':['na']*len(smiles),
                     'CH3/CH2':['na']*len(smiles),
                     'Avg. Length Uninterupted CH2':['na']*len(smiles),
                     '(cycloalkane) >CH2':['na']*len(smiles),
                     '(cycloalkane) >CH-':['na']*len(smiles),
                     '(cycloalkane) >C<':['na']*len(smiles),
                     '(aromatic) =CH-':['na']*len(smiles),
                     '(aromatic) =C<':['na']*len(smiles),
                     '# Aro Rings':['na']*len(smiles),
                     '# Cyclo Rings':['na']*len(smiles),
                     'Avg Cyclo Ring Size':['na']*len(smiles),
                     'Avg Subs per Aro':['na']*len(smiles),
                     'Avg Subs per Cyclo':['na']*len(smiles),
                     'Avg Len Aro Sub':['na']*len(smiles),
                     'Avg Len Cyclo Sub':['na']*len(smiles)})

#%% Get functional groups
for ii in range(len(smiles)):
    try:
        # Molecule from canonical smile
        mol=Chem.MolFromSmiles(smiles[ii])
        
        # -CH3 (linear)
        ch3_linear=group_parse('[CH3;!R]',mol)
        funct_groups['(linear) -CH3'][ii]=ch3_linear
        
        # -CH2- (linear)
        ch2_linear=group_parse('[CH2;!R]',mol)
        funct_groups['(linear) -CH2-'][ii]=ch2_linear
        
        # >CH- (linear)
        ch_linear=group_parse('[CH;!R]',mol)
        funct_groups['(linear) >CH-'][ii]=ch_linear
        
        # >C< (linear)
        c_linear=group_parse('[CH0;!R]',mol)
        funct_groups['(linear) >C<'][ii]=c_linear
        
        # CH3/CH2 (linear)
        ch3_ch2=div_zero(ch3_linear,ch2_linear)
        funct_groups['CH3/CH2'][ii]=ch3_ch2   
     
        # Average uninterupted CH2 chain length
        ch2_length=len_ch2(mol)
        funct_groups['Avg. Length Uninterupted CH2'][ii]=ch2_length
        
        # >CH2 (cycloalkane)   
        ch2_cyclo=group_parse('[CH2;R]',mol)
        funct_groups['(cycloalkane) >CH2'][ii]=ch2_cyclo
        
        # >CH- (cycloalkane)
        ch_cyclo=group_parse('[CH;R]',mol)
        funct_groups['(cycloalkane) >CH-'][ii]=ch_cyclo
        
        # >C< (cycloalkane)
        c_cyclo=group_parse('[CH0;R]',mol)
        funct_groups['(cycloalkane) >C<'][ii]=c_cyclo
        
        # =CH- (aromatic)
        ch_aro=group_parse('[cH]',mol)
        funct_groups['(aromatic) =CH-'][ii]=ch_aro
        
        # =C< (aromatic)
        c_aro=group_parse('[cH0]',mol)
        funct_groups['(aromatic) =C<'][ii]=c_aro
     
        # Ring info
        ring_data=ring_info(mol)
        funct_groups['# Aro Rings'][ii]=ring_data[0]
        funct_groups['# Cyclo Rings'][ii]=ring_data[1]
        funct_groups['Avg Cyclo Ring Size'][ii]=ring_data[2]
        
        # Ring substitution position info
        sub_data=sub_info(mol)
        funct_groups['Avg Subs per Aro'][ii]=sub_data[0]
        funct_groups['Avg Subs per Cyclo'][ii]=sub_data[1]
        
        # Chain length info
        chain_data=chain_info(mol)
        funct_groups['Avg Len Aro Sub'][ii]=chain_data[0]
        funct_groups['Avg Len Cyclo Sub'][ii]=chain_data[1]
        
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message,'ITERATION ',ii)
    
#%% Merge dataframes
df_out=pd.concat([original,funct_groups],axis=1)

#%% Export
df_out.to_excel(out_file,na_rep=0,index=False)