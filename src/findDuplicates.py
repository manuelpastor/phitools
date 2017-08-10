#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Tool for finding duplicates in an SDFile
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

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem,Draw,Descriptors
from SDFhelper import *
import os, sys, argparse

def findDuplicates (args): #sdf, name, out):

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    
    suppl = Chem.ForwardSDMolSupplier(args.sdf,removeHs=False, sanitize=False)

    idlist = []
    nmlist = []
    smlist = []

    sys.stdout.write('reading SDFile...\n')
    counter = 0
    for mol in suppl:
        counter+=1
        
        if mol is None: continue
        try:
            inchi = Chem.MolToInchi(mol)
            inkey = Chem.InchiToInchiKey(inchi)
            smiles = Chem.MolToSmiles(mol)
        except:
            continue

        name = getName(mol, counter)

        idlist.append(inkey[:-3])
        nmlist.append(name)
        smlist.append(smiles)
    args.sdf.close()
    
    n = len(idlist)

    sys.stdout.write('analizing duplicates...\n')

    args.out.write('{}\n'.format('\t'.join(['i', 'j', 'namei', 'namej', 'smilesi', 'smilesj'])))  # Header
    duplicates = 0
    for i in range (n):
        for j in range (i+1,n):
            if idlist[i]==idlist[j]:
                args.out.write('{}\n'.format('\t'.join([str(i), str(j), nmlist[i], nmlist[j], smlist[i], smlist[j]])))
                duplicates+=1
    args.out.close()

    sys.stdout.write('\n%d duplicate molecules found\n' %duplicates)

def main ():

    parser = argparse.ArgumentParser(description='Find duplicated molecules. In the output file, the first columns present the properties of the first molecule duplicated, the last columns contain data about the second molecule identified.')
    parser.add_argument('-f', '--sdf', type=argparse.FileType('rb'), help='SD file', required=True)
    parser.add_argument('-i', '--id', type=str, default='database_substance_id', help='moleculeID')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default='duplicates.txt', help='Output file name (default: duplicates.txt)')
    args = parser.parse_args()

    findDuplicates(args)
    
if __name__ == '__main__':    
    main()
