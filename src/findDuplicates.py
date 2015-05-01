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
import os
import sys
import getopt

def findDuplicates (sdf, name, out):

    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    
    suppl = Chem.SDMolSupplier(sdf,removeHs=False, sanitize=False)

    idlist = []
    nmlist = []
    smlist = []

    print 'reading SDFile...'
    counter = 0
    for mol in suppl:

        counter+=1
        
        if mol is None: continue
        try:
            inchi = Chem.MolToInchi(mol)
            inkey = Chem.InchiToInchiKey(inchi)
            smile = Chem.MolToSmiles(mol)
        except:
            continue

        try:
            ni = mol.GetProp(name)
        except:
            ni = 'mol%0.8d' %counter

        idlist.append(inkey[:-3])
        nmlist.append(ni)
        smlist.append(smile)
    
    n = len(idlist)

    print 'analizing duplicates...'

    fo = open (out,'w+')
    fo.write('i\tj\tnamei\tnamej\tsmilesi\tsmilesj\n')
    duplicates = 0
    for i in range (n):
        for j in range (i+1,n):
            if idlist[i]==idlist[j]:
                line=str(i)+'\t'+str(j)+'\t'+nmlist[i]+'\t'+nmlist[j]+'\t'+smlist[i]+'\t'+smlist[j]
                fo.write(line+'\n')
                duplicates+=1
    fo.close()

    print '\n%d duplicate molecules found' %duplicates

               
def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'findDuplicates -f file.sdf [-n database_substance_id] [-o output.csv]'
    print '\t the output file csvfile.csv contains tab separated info'
    print '\t the first columns present the properties of the first molecule duplicated, the last columns contain data about the second molecule identified.'    
    sys.exit(1)

def main ():
    
    sdf = None
    out = 'output.csv'
    name = 'database_substance_id'

    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:o:n:', [])
    except getopt.GetoptError:
       usage()
       print "False, Arguments not recognized"
       sys.exit(1)

    if args:
       usage()
       print "False, Arguments not recognized"

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-f':
                sdf = arg
            elif opt in '-o':
                out = arg
            elif opt in '-n':
                name = arg
                
                
    if sdf is None:
        usage()

    findDuplicates(sdf, name, out)
    
if __name__ == '__main__':    
    main()
