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

from rdkit import Chem
import os
import sys
import argparse
import re

sep = '\t'

def addData (args):
    #args.sdf, args.data, args.id, args.out
    #sdf_id, data_id, cmpdID, out_id):

    ###################
    # Read data file
    ###################
    tags = args.data.readline().rstrip().split(sep) 
    vals = []
    for line in args.data:
        line=line.rstrip().split(sep)
        vals.append(line)            
    args.data.close()  

    ###################
    # Search property #
    ###################

    cmpdID = args.id.replace(" ","") # delete whitespaces in the compound identifier 
    
    if not cmpdID in tags:        
        for i in tags:
            match=re.search(cmpdID[:5], i)
            if match:
                sys.stderr.write("Property descriptor not found, please try again. Similar property found:" + i)
                return
            else:
                match=re.search(cmpdID[len(cmpdID)-5:], i)
                if match:
                    sys.stderr.write("Property descriptor not found, please try again. Similar property found:" + i)
                    return                
        #not found 
        if i == tags[-1]:
            sys.stderr.write("property descriptor not found")
            return
    else:
        ind=tags.index(cmpdID)
        
    ###################
    # Process SDFile  #
    ###################
  
    suppl = Chem.SDMolSupplier(sdf_id, removeHs=False, sanitize=False)
    print("Input file has",len(suppl),"molecules")
    fo = open (out_id,'w+')
    
    for mol in suppl:
    
        if mol is None: continue

        # Add MolBlock
        fo.write(Chem.MolToMolBlock(mol, kekulize=False))

        l = []
        for i in mol.GetPropNames():
            l.append(i)
        if not cmpdID in l: continue

        db_id=mol.GetProp(mol.GetPropNames()[l.index(cmpdID)])

        match = False
        for i in vals:
            if i[ind]==db_id:
                match = True
                for j in range(len(i)):                          
                    fo.write('>  <'+tags[j]+'>\n'+i[j]+'\n\n')               
        if not match:
            for t in tags:                          
                fo.write('>  <'+t+'>\nna\n\n')  

        fo.write('$$$$\n')
        
    fo.close()
    
    suppl1 = Chem.SDMolSupplier(out_id)
    print("Output file has",len(suppl1),"molecules")

def main ():

    parser = argparse.ArgumentParser(description='Add data from an input table into SD file fields. The data file must be a tab separated file with a single line header. One of the columns must contain a unique id, present also in the SDFile, which is used for the matching. This field can be specified using the parameter -i | --id.')
    parser.add_argument('-f', '--sdf', type=argparse.FileType('r'), help='SD file', required=True)
    parser.add_argument('-d', '--data', type=argparse.FileType('r'), help='Data file', required=True)
    parser.add_argument('-i', '--id', type=str, help='moleculeID')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default='output.sdf', help='Output file name (default: output.sdf)')
    args = parser.parse_args()

    addData(args)    

if __name__ == '__main__':    
    main()
