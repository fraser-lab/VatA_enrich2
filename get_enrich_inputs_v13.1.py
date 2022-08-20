#!/usr/bin/env python
# coding: utf-8

"""
Jenna Pellegrino

python3 get_enrich_inputs_v13.1.py

requires having (A) Expected_mutations.csv and (B) VatA_ref.fa

(A) this has two columns:
    1. position of mutated amino acid
    2. successfully mutated-to codon, which is present in the library
(B) this is the DNA sequence of the WT VatA protein

hardcoded are directories for where to find input files and place output files

"""

import pandas as pd
import os
import re #regex
import json

codon_code = {
            'Cys': ['TGT', 'TGC'],
            'Asp': ['GAT', 'GAC'],
            'Ser': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
            'Gln': ['CAA', 'CAG'],
            'Met': ['ATG'],
            'Asn': ['AAC', 'AAT'],
            'Pro': ['CCT', 'CCG', 'CCA', 'CCC'],
            'Lys': ['AAG', 'AAA'],
            'Stop': ['TAG', 'TGA', 'TAA'],
            'Thr': ['ACC', 'ACA', 'ACG', 'ACT'],
            'Phe': ['TTT', 'TTC'],
            'Ala': ['GCA', 'GCC', 'GCG', 'GCT'],
            'Gly': ['GGT', 'GGG', 'GGA', 'GGC'],
            'Ile': ['ATC', 'ATA', 'ATT'],
            'Leu': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
            'His': ['CAT', 'CAC'],
            'Arg': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
            'Trp': ['TGG'],
            'Val': ['GTA', 'GTC', 'GTG', 'GTT'],
            'Glu': ['GAG', 'GAA'],
            'Tyr': ['TAT', 'TAC']}

def reverse_look(codon_code,codon):
    for key, values in codon_code.items():
        for value in values:
            if value == codon:
                AA = key
    return AA

def enrich_t0_files(loc_t0_files,save_csvs):
    """
    GET ENRICH INPUTS OF TIMEPOINT ZERO (t0) FILES FIRST 
    1) Parse given directory with aligned files and identify t0 files
    2) Store (rep, t0_enrich file name) to dictionary called "rep_t0"
    3) Run t0 files through get_enrich_inputs; flag for is_t0 should be True
    4) Return the 'rep_t0' dictionary to be used in enrich_non_t0_files
    
    N.B. In order to skip saving csv's (False), you need to have the t0_enrich csv's
    in the loc_output_dir. Otherwise, enrich_non_t0_files won't have t0's to look at.
    """   
    rep_t0={} #key = rep; value = t0_enrich file

    for file in os.listdir(loc_t0_files): 
        if file.endswith('mutlist.csv') and "t0" in file:
            outname=file[0:-15] #removes the _R#_mutlist.csv portion of the name
            rep = re.search('(rep[0-9])',file).groups()[0]
            rep_t0[rep]=(loc_output_dir+outname+"_enrich.csv")

            if save_csvs==True:
                df=get_enrich_input(loc_t0_files+file,True) #given mutlist.csv, returns a dataframe
                df.to_csv(loc_output_dir+outname+"_all.csv",index=False, sep="\t")
                df=df.drop(columns=["wt_codon","wt_aa","pos","mut_codon","mut_aa"])
                df.to_csv(loc_output_dir+outname+"_enrich.csv",index=False, sep="\t")
                
    return rep_t0

def enrich_non_t0_files(loc_aligned_files,save_csvs):
    """
    GET ENRICH INPUTS OF ALL OTHER TIMEPOINTS
    1) Parse given directory with aligned files and identify non-t0 files
    2) Run non-t0 files through get_enrich_inputs; flag for is_t0 should be False
    """

    for file in os.listdir(loc_aligned_files):    
        if file.endswith('mutlist.csv') and "t0" not in file:
            outname=file[0:-15] #removes the _R#_mutlist.csv portion of the name       

            df=get_enrich_input(loc_aligned_files+file,False) #given mutlist.csv, returns a dataframe
            
            if save_csvs==True:
                df.to_csv(loc_output_dir+outname+"_all.csv",index=False, sep="\t")
                df=df.drop(columns=["wt_codon","wt_aa","pos","mut_codon","mut_aa"])
                df.to_csv(loc_output_dir+outname+"_enrich.csv",index=False, sep="\t")

def get_enrich_input(mutlist,is_t0):
    print("Parsing",mutlist)
    
    counts=pd.read_csv(mutlist,header=None)
    counts.columns=["pos","mut_codon","count"]
    counts["observed"]=counts["pos"].astype(str)+counts["mut_codon"]
    
    df=pd.DataFrame(columns=["wt_codon","wt_aa","pos","mut_codon","mut_aa","variant","count"])

    
    """
    Go through each row of 'counts' and check if the observed mutation is an expected mutation (will appear in 'twist').
    If mutation is expected, add the hgvs and count to the 'df'.
    Meanwhile, keep track of the counts of all wildtypes.
    """
    
    wt_counts=0
    wt_found=0
    for row in range(0,counts.shape[0]):
        pos=counts.at[row,"pos"]
        mut_codon=counts.at[row,"mut_codon"]    
        observed=counts.at[row,"observed"]

        if twist["expected"].str.fullmatch(observed).any(): #check if observed is in the expected list; fullmatch looks for exact matches
            #the codon is expected at the position; the count is included

            count=counts.at[row,"count"]
            wt_codon=wt_codons[pos-1]
            wt_aa=reverse_look(codon_code,wt_codon)
            mut_aa=reverse_look(codon_code, mut_codon)
            hgvs="p."+wt_aa+str(pos)+mut_aa

            #print(wt_codon,wt_aa,pos,mut_codon,mut_aa,hgvs,count)

            if wt_aa == mut_aa: #using codons here isn't correct, since Twist optimized for E. coli codons; only 12 matches by codon
                wt_counts+=count
                wt_found+=1
                #print("WT detected at pos", pos)

            df.loc[df.shape[0]] = [wt_codon,wt_aa,pos,mut_codon,mut_aa,hgvs,count] #add counts to df

            
        #else, the codon is not expected at the position; the count is trashed
    
    if is_t0 == False:
        print("I am mutlist",mutlist)
        rep = re.search('(rep[0-9])',mutlist).groups()[0]
        df=find_missing_expected_muts(df,rep)
        df=remove_hgvs_counts_not_in_t0(df,rep)
        
    
    """
    Renumber the index numbers so that there are no missing row numbers post removing hgvs counts not in t0.
    NB: If you don't do this, adding the last _wt will overwrite the last row.
    """
    df.reset_index(drop=True, inplace=True) # renumber rows beginning at 0 AND get rid of the old index
     
    
    """
    Finally, add the wildtype counts to the end of 'df'.
    """     
    #                        [wt_codon,wt_aa,pos,mut_codon,mut_aa, hgvs,count]  
    df.loc[df.shape[0]] = [        0,    0,  0,        0,     0,"_wt",wt_counts]

    return df

def find_missing_expected_muts(df,rep):
    """
    1) Using the rep#, look up the appropriate corresponding t0_enrich file.
    2) Check if any expected variants from t0 of a certain rep are missing in 'df'.
    3) If so, add the variant's hgvs with a count of 0 to 'df'.
    """
    print("I am t0", rep_t0[rep])
    
    corresponding_t0_enrich=twist = pd.read_csv(rep_t0[rep],header=0,delimiter='\t')
    for row in range(0,corresponding_t0_enrich.shape[0]-1): #-1 bc last row is _wt and that gets added after
        expected=corresponding_t0_enrich.at[row,"variant"]
        
        #check if expected is NOT in the observed variant list in df, then add it to df w a count of 0
        if not df["variant"].str.fullmatch(expected).any(): #fullmatch looks for exact matches
            df.loc[df.shape[0]] = [        0,    0,  0,        0,     0,expected,0]
        
    return df

def remove_hgvs_counts_not_in_t0(df,rep):
    """
    THROW OUT non-t0 hgvs rows that are not in t0 of that rep
    1) Using the rep#, look up the appropriate corresponding t0_enrich file.
    2) Check if any observed variants from non-t0 'df' of the rep are missing in the t0 of that rep.
    3) If so, remove the row from the 'df'.
    """
    corresponding_t0_enrich=pd.read_csv(rep_t0[rep],header=0,delimiter='\t')
        
    for row in range(0,df.shape[0]):
        observed=df.at[row,"variant"]
        
        #check if observed is NOT in the expected t0 variant list, then delete the row containing the observed hgvs
        if not corresponding_t0_enrich["variant"].str.fullmatch(observed).any(): #fullmatch looks for exact matches
            df=df.drop(row)
            print('not found in t0 of',rep,":", observed)
            
    return df

loc_aligned_files='./subsection/'
loc_t0_files='./aligned_files/'
loc_output_dir='./subsection/'
twist_expected_muts='./Expected_mutations.csv'
WT_DNA_seq='./VatA_ref.fa'

"""
EXTRACT WT CODONS
1) Read WT DNA sequence
2) Break into 3-letter chunks (codons)
"""
with open(WT_DNA_seq,'r') as reader: #edited file to delete \n; maybe doesn't matter
    wt = reader.read().upper()

wt_codons=[wt[i:i+3] for i in range(0, len(wt), 3)] # splits input str into chunks 3 characters long

"""
CREATE EXPECTED TWIST POS/CODON TABLE
1) Read in the Expected Mutations based on codons that Twist installed at each position
2) Convert to hgvs nomenclature and add to a column called "expected_hgvs" in 'twist'
"""
twist = pd.read_csv(twist_expected_muts,header=0)
twist["expected"]=twist["AA_position"].astype(str)+twist["variant_codon"]

for row in range(0,twist.shape[0]):
    pos=twist.at[row,"AA_position"]
    wt_codon=wt_codons[pos-1]
    mut_codon=twist.at[row,"variant_codon"]
    wt_aa=reverse_look(codon_code,wt_codon)
    mut_aa=reverse_look(codon_code, mut_codon)
    expected_hgvs="p."+wt_aa+str(pos)+mut_aa
    twist.at[row,"expected_hgvs"]=expected_hgvs

rep_t0=enrich_t0_files(loc_t0_files,False) #True/False to 'save_csvs'
enrich_non_t0_files(loc_aligned_files,True) #True/False to 'save_csvs'; False is good for testing things