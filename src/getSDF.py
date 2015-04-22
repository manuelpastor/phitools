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

import urllib
import os
import sys
import getopt
from rdkit import Chem
from rdkit.Chem import AllChem

def getSDF(out,data,incl_sml):

    vals = []
    
    f = open (data)
    for line in f:
        line=line.rstrip()
        line=line.split('\t')
        vals.append(line[0])
    f.close()

    fo = open (out,'w+')     
    for v in vals:
        if not incl_sml:
            print v, 
            smi = urllib.urlopen('http://cactus.nci.nih.gov/chemical/structure/'+v+'/smiles')
            smi1 = smi.readline().rstrip()
            
        else:
            smi1 = v

        if smi1:
            print smi1
        else:
            print 'error processing ', v
            continue

        try:
            m = Chem.MolFromSmiles(smi1)
            Chem.AllChem.Compute2DCoords(m)
        except:
            print 'error processing ', v
            continue
            
        mb = Chem.MolToMolBlock(m)        

        fo.write(mb)
        
        if not incl_sml:
            fo.write('>  < name >\n'+v+'\n\n')
            
        fo.write('>  < smiles >\n'+smi1+'\n\n')

        fo.write('$$$$\n')        
    fo.close()

def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'getSDF -f|s input.csv [-o output.sdf]'
    print '\n\t -f input.csv (list of compound names)'
    print '\t -s input.csv (list of SMILES)'
    print '\t -o output.sdf (output SDFile)'
    sys.exit(1)

def main ():
    out= 'output_getSDF.sdf'
    data = None
    incl_sml = False
 
    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:s:o:')
    except getopt.GetoptError:
       usage()
       print "False, Arguments not recognized"
    
    if not len(opts):
       usage()
       print "False, Arguments not recognized"

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-o':
                out = arg
            elif opt in '-f':
                data = arg
            elif opt in '-s':
                data = arg
                incl_sml = True

    if not data: usage()
    
    getSDF(out,data,incl_sml)    

if __name__ == '__main__':    
    main()
