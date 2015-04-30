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
import getopt
import re

def addData (sdf_id, data_id, prop_id, out_id):

    vals = []
    tags = []
    prop_id=prop_id.replace(" ","") # delete whitespaces

    ###################
    # Read csv file
    ###################
    f = open (data_id)
    for line in f:
        line=line.rstrip()
        line=line.split('\t')
        if len(tags) == 0:
            tags.append(line)            
        else:
            vals.append(line)            
    f.close()   

    ###################
    # Search property #
    ###################
    
    if not prop_id in tags[0]:        
        for i in tags[0]:
            match=re.search(prop_id[:5], i)
            if match:
                print "Property descriptor not found, please try again. Similar property found:" + i
                return
            else:
                match=re.search(prop_id[len(prop_id)-5:], i)
                if match:
                    print "Property descriptor not found, please try again. Similar property found:" + i
                    return                
        #not found 
        if i == tags[0][-1]:
            print "property descriptor not found"
            return
    else:
        ind=tags[0].index(prop_id)
        
    ###################
    # Process SDFile  #
    ###################
  
    suppl = Chem.SDMolSupplier(sdf_id, removeHs=False, sanitize=False)
    print "Input file has",len(suppl),"molecules"
    fo = open (out_id,'w+')
    
    for mol in suppl:
    
        if mol is None: continue

        # Add MolBlock
        fo.write(Chem.MolToMolBlock(mol, kekulize=False))

        l = []
        for i in mol.GetPropNames():
            l.append(i)
        if not prop_id in l: continue

        db_id=mol.GetProp(mol.GetPropNames()[l.index(prop_id)])

        match = False
        for i in vals:
            if i[ind]==db_id:
                match = True
                for j in range(len(i)):                          
                    fo.write('>  <'+tags[0][j]+'>\n'+i[j]+'\n\n')               
        if not match:
            for t in tags[0]:                          
                fo.write('>  <'+t+'>\nna\n\n')  

        fo.write('$$$$\n')
        
    fo.close()
    
    suppl1 = Chem.SDMolSupplier(out_id)
    print "Output file has",len(suppl1),"molecules"

def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'addDataToSDF -f file.sdf -d data.csv --id=molecule_id [-o output.sdf]'
    print '\n\t data.csv contains tab separated info and a single line header'
    print '\t one of the columns must contain a unique id, present also in the SDFile, which is used for the matching'
    print '\t this field can be specified using the parameter --id '
    print '\t the output file will include all molecules present in the original SDFile'
    sys.exit(1)

def main ():
    
    sdf = None
    data = None
    prop = None
    out = 'output.sdf'
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:d:o:', ['id='])
    except getopt.GetoptError:
       usage()
       print "Error. Arguments not recognized"

    if args:
       usage()
       print "Error. Arguments not recognized"

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-f':
                sdf = arg
            elif opt in '-d':
                data = arg
            elif opt in '--id':
                prop = arg
            elif opt in '-o':
                out = arg

    if sdf==None or data==None or prop==None:
        usage()
        print "Error. Missing arguments"

    addData(sdf, data, prop, out)    

if __name__ == '__main__':    
    main()
