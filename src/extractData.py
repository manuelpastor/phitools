#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Tool for extracting data within an SDFile
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
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
from rdkit import Chem

def extractNames (mol):

    suppl=Chem.SDMolSupplier(mol)

    count = 1
    while True:
        try:
            mi = suppl.next()
        except:
            break

        if not mi:
            break

        name = ''
        if mi.HasProp ('name'):
            name = mi.GetProp('name')
        if not name:
            name = mi.GetProp('_Name')
        if not name:
            name = 'mol%0.8d'%count
        count +=1

        print name

def extractField (mol,field):

    suppl=Chem.SDMolSupplier(mol)

    count = 0
    while True:
        try:
            mi = suppl.next()
        except:
            break

        if not mi:
            break

        name = ''
        if mi.HasProp (field):
            name = mi.GetProp(field)
        else:
            name = 'na'
            
        print name

def extractAll (mol):

    suppl = Chem.SDMolSupplier(mol)
    
    isFirst = True
    header = []
    for m in suppl:

        if m is None: continue
        
        if isFirst:
            for i in m.GetPropNames():
                header.append(i)
                print i+'\t',
            print
            isFirst = False

        for i in header:
            try:
                print m.GetProp(i)+'\t',
            except:
                print 'na\t',
        print


    
def usage ():
    """Prints in the screen the command syntax and argument"""
    
    print 'extractData -f sdfile.sdf [--name|--field=Activ|--table]'

def main ():

    mol = None
    action = None
    data = None
    
    try:
       opts, args = getopt.getopt(sys.argv[1:], 'f:',
                                  ['name','add=','field=','table'])

    except getopt.GetoptError:
       print('Error. Arguments not recognized')
       usage()
       sys.exit(1)

    if args:
       print('Error. Arguments not recognized')
       usage()
       sys.exit(1)
        
    if len( opts ) > 0:
        for opt, arg in opts:
                
            if opt in '-f':
                mol = arg
                
            elif opt in '--name':
                action = 'name'
                
            elif opt in '--field':
                action = 'field'
                field = arg

            elif opt in '--table':
                action = 'table'
                
            elif opt in '-h':
                usage()
                sys.exit(0)

    #print action, mol, data

    if action == None:
        usage()
        sys.exit(1)
    elif action == 'name':
        if mol == None:
            sys.exit(1)
            
        extractNames (mol)

    elif action == 'field':
        if mol== None or field== None:
            usage()
            sys.exit(1)
            
        extractField (mol, field)
        
    elif action == 'table':
        if mol== None:
            usage()
            sys.exit(1)
            
        extractAll (mol)
        
    sys.exit(0)
        
if __name__ == '__main__':
    
    main()
