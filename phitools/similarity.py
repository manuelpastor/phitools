#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Add Data to SDFile
##                   
##    Authors:       Inés Martínez (mmartinez4@imim.es)
##                   Manuel Pastor (manuel.pastor@upf.edu)
##
##    Copyright 2015 Manuel Pastor
##
##    This file is part of PhiTools
##
##    PhiTools is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    PhiTools is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhiTools.  If not, see <http://www.gnu.org/licenses/>


from rdkit import DataStructs
from phitools import moleculeHelper as mh
import pandas as pd
import argparse, sys

sep = '\t'

def compareMax(fpA, fpB=None, cutoff=None):
    
    #################################
    ### Get compound similarities ###
    #################################

    simD = {}
    startSim = 0

    namesA = list(fpA.keys())
    nA = len(namesA)

    # Work with only one input file
    if fpB is None:
        for name in namesA:
            simD[name] = [None, None, None, startSim]

        for i in range(nA):
            name1 = namesA[i]
            [fp1, smiles1] = fpA[name1]
            for j in range(i+1, nA):
                name2 = namesA[j]
                [fp2, smiles2] = fpA[name2]
                sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                if cutoff is not None and sim < cutoff:
                    continue
                
                if sim > simD[name1][-1]: 
                    simD[name1] = [name2, smiles1, smiles2, sim]
                if sim > simD[name2][-1]: 
                    simD[name2] = [name1, smiles2, smiles1, sim]

    # Work with two input files
    else:
        namesB = list(fpB.keys())
        nB = len(namesB)

        for name in namesA:
            simD[name] = [None, None, None, startSim]
        for name in namesB:
            simD[name] = [None, None, None, startSim]

        for i in range(nA):
            name1 = namesA[i]
            [fp1, smiles1] = fpA[name1]            
            for j in range(nB):
                name2 = namesB[j]
                [fp2, smiles2] = fpB[name2]

                sim = DataStructs.TanimotoSimilarity(fp1, fp2) #DataStructs.DiceSimilarity(fp1, fp2)
                if cutoff is not None and sim < cutoff:
                    continue

                if sim > simD[name1][-1]: 
                    simD[name1] = [name2, smiles1, smiles2, sim]
    
    # Remove compounds with no neighbours over the cutoff
    d = {x:simD[x] for x in simD if simD[x][0] is not None}

    # Convert dictionary to pandas dataframe
    df = pd.DataFrame.from_dict(d, orient='index')
    df.reset_index(level=0, inplace=True)
    df.columns = ['cmpd1', 'cmpd2', 'smiles1', 'smiles2', 'similarity']

    return df

def compareMin(fpA, fpB=None, cutoff=None):
    
    #################################
    ### Get compound similarities ###
    #################################

    simD = {}
    startSim = 1

    namesA = list(fpA.keys())
    nA = len(namesA)

    # Work with only one input file
    if fpB is None:

        for i in range(nA):
            simD[name] = [None, None, None, startSim]
            name1 = namesA[i]
            [fp1, smiles1] = fpA[name1]
            for j in range(i+1, nA):
                name2 = namesA[j]
                [fp2, smiles2] = fpA[name2]
                sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                if cutoff is not None and sim < cutoff:
                    continue
                
                if sim < simD[name1][-1]: 
                    simD[name1] = [name2, smiles1, smiles2, sim]
                if sim < simD[name2][-1]: 
                    simD[name2] = [name1, smiles2, smiles1, sim]

    # Work with two input files
    else:
        namesB = list(fpB.keys())

        for name1 in namesA:
            simD[name1] = [None, None, None, startSim]
            [fp1, smiles1] = fpA[name1]            
            for name2 in namesB:
                name2 = namesB[j]
                [fp2, smiles2] = fpB[name2]

                sim = DataStructs.TanimotoSimilarity(fp1, fp2) #DataStructs.DiceSimilarity(fp1, fp2)
                if cutoff is not None and sim < cutoff:
                    continue

                if sim < simD[name1][-1]: 
                    simD[name1] = [name2, smiles1, smiles2, sim]

    # Remove compounds with no neighbours over the cutoff
    d = {x:simD[x] for x in simD if simD[x][0] is not None}

    # Convert dictionary to pandas dataframe
    df = pd.DataFrame.from_dict(d, orient='index')
    # reset_index fails in pandas 0.20.1
    #df.reset_index(level=0, inplace=True)
    df['cmpd1'] = df.index
    df.columns = ['cmpd2', 'smiles1', 'smiles2', 'similarity', 'cmpd1']
    df = df[['cmpd1', 'cmpd2', 'smiles1', 'smiles2', 'similarity']]

    return df

def compareAll(fpA_dict, fpB_dict=None, cutoff=None):
    
    #################################
    ### Get compound similarities ###
    #################################

    simD = {}
    namesA = list(fpA_dict.keys())
    nA = len(namesA)

    # Work with only one input file
    if fpB is None:
        for i in range(nA):
            name1 = namesA[i]
            simD[name1] = {}
            [fp1, smiles1] = fpA_dict[name1]
            for j in range(i+1, nA):
                name2 = namesA[j]
                [fp2, smiles2] = fpA_dict[name2]
                sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                if cutoff is not None and sim < cutoff:
                    simD[name1][name2] = None
                    simD[name2][name1] = None
                else:                                    
                    simD[name1][name2] = [smiles1, smiles2, sim]
                    simD[name2][name1] = [smiles2, smiles1, sim]

    # Work with two input files
    else:
        namesB = list(fpB_dict.keys())

        for nameA in namesA:
            simD[nameA] = {}
            [fpA, smilesA] = fpA_dict[nameA]
            for nameB in namesB:
                [fpB, smilesB] = fpB_dict[nameB]
                sim = DataStructs.TanimotoSimilarity(fpA, fpB) #DataStructs.DiceSimilarity(fp1, fp2)

                if cutoff is not None and sim < cutoff:
                    simD[nameA][nameB] = None
                else:
                    simD[nameA][nameB] = [smilesA, smilesB, sim]

    # Remove compounds with no neighbours over the cutoff
    d = {x:simD[x] for x in simD if simD[x] is not None}

    # Convert dictionary to pandas dataframe
    df = pd.DataFrame.from_records([[i, j] + d[i][j] for i in d for j in d[i]])
    df.columns = ['cmpd1', 'cmpd2', 'smiles1', 'smiles2', 'similarity']

    return df


def compare(args):

    ###########################
    ### Store the compounds ###
    ###########################
    fpType = args.descriptor[0:4]
    fpRadius = int(args.descriptor[4:])

    fpA = mh.getFPdict (args.format, args.filea, molID= args.id, smilesI= args.col, header= args.header, fpType= fpType, radius= fpRadius)

    if args.fileb is not None:
        fpB = mh.getFPdict (args.format, args.fileb, molID= args.id, smilesI= args.col, header= args.header, fpType= fpType, radius= fpRadius)

    if args.sim == 'max':
        df = compareMax(fpA, fpB, args.cutoff)
    elif args.sim == 'min':
        df = compareMin(fpA, fpB, args.cutoff)
    elif args.sim == 'all':
        df = compareAll(fpA, fpB, args.cutoff)

    df.to_csv(args.out, sep= '\t', header= True, index=False)
                        


def main ():
    parser = argparse.ArgumentParser(description='Get the smilarity between the compounds in the input file if only one is provided or between files if two are provided.')
    parser.add_argument('-a', '--filea', type=argparse.FileType('rb'), help='Input file.', required=True)
    parser.add_argument('-b', '--fileb', type=argparse.FileType('rb'), help='Optional input file. If it is provided he compounds in this file will be compared to the compounds in the first input file.')
    parser.add_argument('-f', '--format', action='store', dest='format', choices=['smi', 'sdf'], default='smi', help='Specify the input format (smiles strings (default) or SD file).')
    parser.add_argument('-s', '--sim', action='store', choices=['min', 'max', 'all'], default='max', help='Get only the closest compounds (\'max\', default), the most dissimilar (\'min\'), or all v all similarties (\'all\')')
    parser.add_argument('-d', '--descriptor', action='store', dest='descriptor', default='ecfp4', help='Specify the descriptor to be used in the similarity calculation (ecfp4 (default), fcfp4, ecfp2, etc.).')
    parser.add_argument('-c', '--cutoff', type=float, help='If wanted, set a minimum similarity cutoff.')
    parser.add_argument('-i', '--id', type=str, help='Field containing the molecule ID. If it is not provided for the SD file format, the SD file compound name will be used.')
    parser.add_argument('-x', '--col', type=int, default=1, help='If the input file has smiles, indicate which column contains the smiles strings.')
    parser.add_argument('-n', '--noheader', action='store_false', dest='header', help='Smiles input data file doesn\'t have a header line.')
    parser.add_argument('-o', '--out', type=argparse.FileType('w+'), default='output.txt', help='Output file name (default: output.txt)')
    args = parser.parse_args()
    
    if args.format == 'smi':
        if args.col is not None:
            args.col -= 1

    if args.format == 'smi' and args.id is not None:
        try:
            args.id = int(args.id)-1
        except:
            sys.stderr('The ID argument must be a column index if the input file is of smiles format.\n')
            sys.exit()

    compare(args)    

if __name__ == '__main__':    
    main()