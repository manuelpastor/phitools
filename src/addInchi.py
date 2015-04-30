#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Add an InChI field to a SDFile 
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
from rdkit.Chem import AllChem
import os
import sys
import getopt

def addInchi (sdf,out):

    # Read SDF
    suppl = Chem.SDMolSupplier(sdf,removeHs=False, sanitize=False)

    # Create header
    fo = open(out,'w+')
    
    for mol in suppl:

        if mol is None: continue
        
        fo.write(Chem.MolToMolBlock(mol,kekulize=False))

        pnames = []
        for i in mol.GetPropNames():
            pnames.append(i)

        pvalues = []
        for i in pnames:
            pvalues.append(mol.GetProp(i))
            
        for pn,pv in zip(pnames,pvalues):        
            fo.write('>  <'+pn+'>\n'+pv+'\n\n')
        
        try:
            inchi = Chem.MolToInchi(mol)
            inkey = Chem.InchiToInchiKey(inchi)
        except:
            inchi = 'na'
            inkey = 'na'

        fo.write('>  <inchi>\n'+inchi+'\n\n')
        fo.write('>  <inchikey>\n'+inkey+'\n\n')

        fo.write('$$$$\n')

    fo.close()
 
def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'addInchi [-f file.sdf] [-o output.sdf]'
    sys.exit(1)

def main ():
    sdf = None
    out = 'output.sdf'
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:o:', [])
    except getopt.GetoptError:
       usage()
       print "False, Arguments not recognized"
       sys.exit(1)

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-f':
                sdf = arg
            elif opt in '-o':
                out = arg
                
    if sdf is None:
        usage()

    addInchi (sdf,out)

if __name__ == '__main__':    
    main()
