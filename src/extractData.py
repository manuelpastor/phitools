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
import argparse
from rdkit import Chem

def molName(mol, count):
    name = ''
    if mol.HasProp ('name'):
        name = mol.GetProp('name')
    if not name:
        name = mol.GetProp('_Name')
    if not name:
        name = 'mol%0.8d'%count
    return name

def extractNames (args):
    suppl=Chem.SDMolSupplier(args.sdf.name)
    count = 1
    for m in suppl:
        name = molName(m, count)
        count +=1
        args.out.write('{}\n'.format(name))

def extractField (args):
    suppl=Chem.SDMolSupplier(args.sdf.name)
    for m in suppl:
        if m.HasProp (args.field):
            value = m.GetProp(args.field)
        else:
            value = 'NA'
        args.out.write('{}\n'.format(value))

def extractAll (args):

    fields = set()
    # Cycle through all molecules to make sure all field names are stored
    suppl = Chem.SDMolSupplier(args.sdf.name)
    for m in suppl:
        if m is None: continue
        fields = fields.union(set(m.GetPropNames()))
    fields = list(fields)
    header = ['ID']
    header.extend(fields)
    args.out.write('{}\n'.format('\t'.join(header)))
    
    suppl = Chem.SDMolSupplier(args.sdf.name)
    count = 1
    for m in suppl:
        if m is None: continue
        line = [molName(m, count)]
        for field in fields:
            if m.HasProp (field):
                value = m.GetProp(field)
            else:
                value = 'NA'
            line.append(value)
        args.out.write('{}\n'.format('\t'.join(line)))
        count += 1


def main ():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--sdf', type=argparse.FileType('rb'), help='SD file', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--field', type=str, help='Extract only the data in this field of the SD file.')
    group.add_argument('-a', '--all', action='store_true', help='Extract data in all fields.')
    group.add_argument('-n', '--name', action='store_true', help='Extract molecule names.')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default='output.tsv', help='Output file name (default: output.tsv)')
    args = parser.parse_args()
    args.sdf.close()

    if args.name:
        extractNames (args)
    elif args.field:
        extractField (args)
    else:
        extractAll (args)

    args.out.close()
    
        
if __name__ == '__main__':
    
    main()
