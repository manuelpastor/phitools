#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Split an SDFile randomly in a two datasets
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu)
##
##    Copyright 2016 Manuel Pastor
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
import os, sys, getopt
import numpy as np

def splitSDF (sdf, prop, seed):

    nmols = 0
    try:
        f = open (sdf,'r')
    except:
        print('unable to open file ',sdf, 'ABORT')
        exit(1)
        
    for line in f:
        if line.startswith('$$$$'):
            nmols = nmols+1
    f.close()

    ntrai = int(np.round(prop*nmols/100.0))
    npred = nmols - ntrai
    
    print(nmols, "compounds found. Creating series of ", ntrai, " for training and ", npred, " for prediction")

    if seed != None :
        npseed = int(seed)
        np.random.seed(npseed)
        
    elements = np.random.choice(nmols, ntrai, False)
    #print elements
    
    f  = open (sdf,'r')
    fp = open ('pr-'+sdf,'w')
    ft = open ('tr-'+sdf,'w')

    i = 0
    for line in f:
        if i in elements :
            ft.write(line)
        else:
            fp.write(line)
        if line.startswith('$$$$'):
            i=i+1

    f.close()
    fp.close()
    ft.close()

def usage ():
    """Prints in the screen the command syntax and argument"""
    print('randomSplitSDF -f file.sdf -p 70 [-s 2356]')
    sys.exit(1)

def main ():
    
    sdf = None
    p = None
    seed = None
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:p:s:')
    except getopt.GetoptError:
       usage()
       print("Error. Arguments not recognized")

    if args:
       usage()
       print("Error. Arguments not recognized")

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-f':
                sdf = arg
            elif opt in '-p':
                p = arg
            elif opt in '-s':
                seed = arg

    if sdf==None or p==None:
        print("Error. Missing arguments")
        usage()
        
    try:
        prop = int(p)
    except:
        print("Error. Proportion not an integer")
        usage()
        
        
    splitSDF(sdf, prop, seed)    

if __name__ == '__main__':    
    main()
