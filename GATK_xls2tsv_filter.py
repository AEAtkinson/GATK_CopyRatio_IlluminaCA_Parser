#!/usr/bin/env python3

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 08:29:37 2019 through 2/24/2019

@author: aaronatkinson
"""
#import packages
import pandas as pd
import glob, os

#Assign directory to change or us within directory
directory="/Users/aaronatkinson/Desktop/Avatar_Upload_01_2019/FINAL_CNV_upload/xls2tsv_filter.py"

#create filelist
for filename in glob.glob('**.xls'):
    print(filename)
    open(filename)
    df = pd.read_table(filename)
    print(df)
    
#split value columns for chr, start, end - per illumina tsv format
    new = df["#Chr:Start-End"].str.split(",", n=1, expand = True) 
    df["Start-End"]= new[1]
    df.drop(columns = ["#Chr:Start-End"], inplace = True)
    new = df["Start-End"].str.split(":", n=1, expand = True)
    df["chr"]= new[0]
    df ["start"]= new[1]
    df["chr"] = df["chr"].map(lambda x: x.lstrip('"'))
    new = df["start"].str.split("-", n=1, expand = True)
    df["start"]= new[0]
    df["stop"]= new[1]
    new = df["stop"].str.split('"', n=1, expand = True)
    df["stop"]= new[0]

#new columns/headers - per illumina tsv format
    df["log2_ratio"] = df["Lg2 Mean Tumor CR"]
    df["copy_number_change"] = df["log2_ratio"] - df["Lg2 Mean Normal CR"]
    df["probe_count"] = df["Num CR Points"]

#filter false copy number changes that didn't pass threshold - note ran into difficulty filtering based of booleans hence the mapping
    booleandf = df.select_dtypes(include=[bool])
    booleanDictionary = {True: 'TRUE', False: 'FALSE'}
    df["Pass Thresholds"] = df["Pass Thresholds"].map(booleanDictionary)
    df1 = df[(df['Pass Thresholds']=='TRUE')]
#filter those deletions that have AF point of <6 and that don't lead to an deletion or amplification i.e. are "0"     
    df2 = df1[(df1['Num Tumor AF Points']>6)]
    df3 = df2[(df2['CR Call']!='0')]
#for bonafide deletions that result in a log2 < -2, change to -2 (maximum diploid deletion or -1 for chrY)     
    mask = df3.log2_ratio < -2
    column_name = 'log2_ratio'
    df3.loc[mask, column_name] = -2
    mask2 = df3.copy_number_change < -2
    column_name2 = 'copy_number_change'
    df3.loc[mask2, column_name2] = -2
    chrYmask = (df3.chr == 'chrY') & (df3.log2_ratio < -1)
    column_name3 = 'log2_ratio'
    df3.loc[chrYmask, column_name3] = -1
    chrYmask2 = (df3.chr == 'chrY') & (df3.copy_number_change < -1)
    column_name4 = 'copy_number_change'
    df3.loc[chrYmask2, column_name4] = -1    
        
#drop unused columns - per illumina tsv format
    df3.drop(columns = ["Pass Thresholds", "CR Call", "Num CR Points", "Geometric Mean Tumor CR", "Lg2 Mean Tumor CR", "Lg2 Mean Normal CR", "Lg2 Mean TN CR Ratios", "Num Tumor AF Points", "Mean Tumor AF", "Num Normal AF Points", "Mean TN AFs", "Mean Normal AF", "Genes", "Start-End"], inplace = True)


#print df split from filelist and save new tsv to_csv function ('\tab') - per illumina format
    print(df3)
    csv_file = os.path.splitext(filename)[0]+".tsv"
    df3.to_csv(csv_file, sep='\t', index=False)