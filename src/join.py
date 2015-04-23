#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Inner Join of two CSV tables
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

import os
import sys
import getopt

def join (fileA,fileB,key,out):

    # Open files
    fa = open (fileA,'r')
    fb = open (fileB,'r')
    fo = open (out,'w+')

    #identify key in both files
    ha = fa.readline().strip('\t')
    try:
        indexA = ha.index(key)
    except:
        print 'key not fount in file '+fileA
        sys.exit(1)

    hb = fb.readline().strip('\t')
    try:
        indexB = hb.index(key)
    except:
        print 'key not fount in file '+fileB
        sys.exit(1)

    # read both files
  

    fa.close()
    fb.close()

    # write output
    
    fo.close()
 
def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'join -a fileA.csv -b fileB.csv --id molecule_id [-o output.csv]'
    sys.exit(1)

def main ():
    fileA = None
    fileB = None
    key = None
    out = 'output.csv'
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'a:b:o:', ['id='])
    except getopt.GetoptError:
       usage()
       print "False, Arguments not recognized"
       sys.exit(1)

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-a':
                fileA = arg
            elif opt in '-b':
                fileB = arg
            elif opt in '-o':
                out = arg
            elif opt in '--id':
                key = arg
                
    if fileA is None or fileB is None or key is None:
        usage()

    join (fileA, fileB, key, out)

if __name__ == '__main__':    
    main()
