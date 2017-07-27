#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Add molecular information to a data file
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

def addData (args): #sdf_id, data_id, prop_id, out_id):
    #args.sdf, args.data, args.id, args.out

    ###################
    # Read csv file
    ###################
    #f = open (data_id)
    tags = args.data.readline().rstrip().split(sep) 
    vals = []
    for line in args.data:
        line=line.rstrip().split(sep)
        vals.append(line)            
    args.data.close()  

    ###################
    # Search property #
    ###################
    
    prop_id=args.id.replace(" ","") # delete whitespaces
    if not prop_id in tags:        
        for i in tags:
            match=re.search(prop_id[:5], i)
            if match:
                print("Property descriptor not found, please try again. Similar property found:" + i)
                return
            else:
                match=re.search(prop_id[len(prop_id)-5:], i)
                if match:
                    print("Property descriptor not found, please try again. Similar property found:" + i)
                    return                
        #not found 
        if i == tags[-1]:
            print("property descriptor not found")
            return
    else:
        ind=tags.index(prop_id)
        
    ###################
    # Process SDFile  #
    ###################
  
    suppl = Chem.SDMolSupplier(args.sdf.name, removeHs=False, sanitize=False)
    print("Input file has",len(suppl),"molecule(s)")
    
    for mol in suppl:
    
        if mol is None: continue        

        l = []
        for i in mol.GetPropNames():
            l.append(i)
        if not prop_id in l: continue

        db_id=mol.GetProp(mol.GetPropNames()[l.index(prop_id)])

        for i in vals:
            if i[ind]==db_id:
                # Add MolBlock
                args.out.write(Chem.MolToMolBlock(mol, kekulize=False))
                
                # Add values
                for j in range(len(i)):                          
                    args.out.write('>  <'+tags[j]+'>\n'+i[j]+'\n\n')               
                args.out.write('$$$$\n')
        
    args.out.close()
    
    suppl1 = Chem.SDMolSupplier(args.out.name)
    print("Output file has",len(suppl1),"molecule(s)")

def main ():
    parser = argparse.ArgumentParser(description='Add data from an input SD file\'s fields into a table file. One of the fields must contain a unique id which is used to identify the compounds in the output data file. This field can be specified using the parameter -i | --id.')
    parser.add_argument('-f', '--sdf', type=argparse.FileType('r'), help='SD file', required=True)
    parser.add_argument('-d', '--data', type=argparse.FileType('r'), help='Data file', required=True)
    parser.add_argument('-i', '--id', type=str, help='moleculeID', required=True)
    parser.add_argument('-o', '--out', type=argparse.FileType('w+'), default='output.sdf', help='Output file name (default: output.sdf)')
    args = parser.parse_args()
    args.sdf.close()

    addData(args)    

if __name__ == '__main__':    
    main()
