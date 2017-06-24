#!/usr/bin/env python

'''This script accepts input data frames of data sets and corresponding models
    and then computes predictive information of each model on each test dataset
    .'''
from __future__ import division
#Our standard Modules
import argparse
import numpy as np
import scipy as sp
import sys
#Our miscellaneous functions
import pandas as pd
import sortseq_garbage.utils as utils
import sortseq_garbage.predictiveinfo_dev as predictiveinfo_dev
import pickle

def main(data_dfs,model_dfs,dataname_columns,dataname_rows,modeltypes,
    dicttype,exptype=None,start=0,end=None,neuralnetlist=None):
    
    seq_dict,inv_dict = utils.choose_dict(dicttype)
    nummodels = len(model_dfs)
    MI = np.zeros([nummodels,len(data_dfs)])
    #Different Models = rows, test sets = columns
    for z, model_df in enumerate(model_dfs):
        for q, data_df in enumerate(data_dfs):
            print neuralnetlist[q]
            data_df2 = data_df.copy()
            print modeltypes[z]
            print data_df2         
            #Calc predictive MI without error calculation
            MI[z,q],std = predictiveinfo_dev.main(
                data_df2,model_df,
                start=start,end=end,err=False,neuralnet=neuralnetlist[z])
    output_df = pd.DataFrame(MI)
    output_df.columns = dataname_columns
    output_df = pd.concat(
        [pd.Series(dataname_rows,name='Training/Test'),output_df],axis=1)    
    return output_df

# Define commandline wrapper
def wrapper(args):
    ds = pd.io.parsers.read_csv(args.datasets,delim_whitespace=True)
    models = pd.io.parsers.read_csv(args.models,delim_whitespace=True)
    dfs = []
    mymodels = []
    neuralnetlist = []
    #Create lists of test data sets and corresponding models
    
    
    
    modeltypes = models.Model_Type
    for z,f in ds.iterrows():
        tempdf = pd.io.parsers.read_csv(f[1],delim_whitespace=True)
        dfs.append(tempdf)
    for z,m in models.iterrows():
        if modeltypes[z] == 'neuralnet':
            neuralnetlist.append(True)
            mymodels.append(pickle.load(open(m[1],'rb')))
        elif modeltypes[z] == 'pair':
            mymodels.append(np.genfromtxt(m[1]))
            neuralnetlist.append(False)
        else:     
            mymodels.append(pd.io.parsers.read_csv(m[1],delim_whitespace=True))
            neuralnetlist.append(False)
    print neuralnetlist
    dataname_columns = ds.exp
    dataname_rows = models.exp
    output_df = main(
        dfs,mymodels,dataname_columns,dataname_rows,modeltypes,
        args.type,exptype=args.exptype,start=args.start,end=args.end,neuralnetlist=neuralnetlist)
    if args.out:
        outloc = open(args.out,'w')
    else:
        outloc = sys.stdout
    pd.set_option('max_colwidth',int(1e8)) # Makes sure columns are not truncated
    output_df.to_string(
        outloc, index=False,col_space=10,float_format=utils.format_string)

# Connects argparse to wrapper
def add_subparser(subparsers):
    p = subparsers.add_parser('compare_predictiveinformation_dev')
    p.add_argument(
        '-ds','--datasets',help='''White space delimited file, where the 
         columns are name and file of data sets.''')
    p.add_argument(
        '-mt','--modeltype',default='LinearEmat', help='''Type of model to be evaluated''')
    p.add_argument('-s','--start',type=int,default=0,help ='''
         Position to start your analyzed region''')
    p.add_argument(
        '-nn','--neuralnet',action='store_true',help='''Flag to use if you have a
        neuralnet model''')
    p.add_argument('-e','--end',type=int,default = None, help='''
         Position to end your analyzed region''')
    p.add_argument('-expt','--exptype',default=None,choices=[None,'sortseq_garbage',
        'selex','dms','mpra'])
    p.add_argument('-m','--models',help='''File names containing models to
        evaluate.''')
    p.add_argument('-t', '--type', choices=['dna','rna','protein'],
        default='dna')            
    p.add_argument('-o', '--out', default=None)
    p.set_defaults(func=wrapper)
