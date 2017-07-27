#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Tool for generating a SDFile from diverse sources
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

import urllib.request, urllib.parse, urllib.error
import os
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def getSDF(args):
    #out,data,iformat,idname,header):

    vals = []
    noms = []
    args.data.readline()  # Header
    counter = 1
    query = [line.rstrip().split('\t')[args.column-1] for line in args.data]
    args.data.close()

    # for inchis use SDWriter
    if args.inchis:
        args.out.close()
        writer = Chem.rdmolfiles.SDWriter(args.out.name)
        for q in query:
            try:
                mol = Chem.inchi.MolFromInchi (q)
                AllChem.Compute2DCoords(mol)
                mol.SetProp(args.field,q)
                writer.write(mol)
            except:
                print('error processing ', v)
                pass
        writer.close()

        return

    # for names convert to SMILES, in either case use MolFromSmiles
    #for v,n in zip(vals,noms):
    for q in query:

        if args.identifier:
            try:
                smi = urllib.request.urlopen('http://cactus.nci.nih.gov/chemical/structure/'+q+'/smiles')
            except:
                continue
            smi1 = smi.readline().decode("utf-8").rstrip()
        elif args.smiles:
            smi1 = q

        if smi1 == None or ('Page not found' in smi1):
            continue
        
        try:
            m = Chem.MolFromSmiles(smi1)
            Chem.AllChem.Compute2DCoords(m)
        except:
            print('error processing', q)
            continue
            
        m.SetProp("_Name",q)
        mb = Chem.MolToMolBlock(m)        

        args.out.write(mb)

        args.out.write('>  <'+args.field+'>\n'+q+'\n\n')
        args.out.write('>  <smiles>\n'+smi1+'\n\n')

        args.out.write('$$$$\n')        
    args.out.close()

def main ():

    parser = argparse.ArgumentParser(description='Get the structure of the compounds in the input data file.')
    parser.add_argument('-d', '--data', type=argparse.FileType('r'), required=True, help='File with molecule identifiers.')
    parser.add_argument('-c', '--column', type=int, default=1, help='Column in the input file that contains the molecule identifiers (default= 1).')
    parser.add_argument('-n', '--noheader', action='store_false', dest='header', help='Input data file doesn\'t have a header line.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-x', '--identifier', action='store_true', help='The input format will be compound indentifiers.')
    group.add_argument('-s', '--smi', action='store_true', help='The input format will be a smiles string.')
    group.add_argument('-i', '--inchis', action='store_true', help='The input format will be molecule InChis.')
    parser.add_argument('-f', '--field', type=str, default='name', help='Field in the output SD file that will contain the unique ID (default= name).')
    parser.add_argument('-o', '--out', type=argparse.FileType('w+'), default='output.sdf', help='Output file name (default: output.sdf)')
    args = parser.parse_args()
    
    getSDF(args)    

if __name__ == '__main__':    
    main()
