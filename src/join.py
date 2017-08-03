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

import os, sys, argparse

sep = '\t'

def join (args):

    # Open files
    if args.type == 'right':
        fa = args.fileb
        fb = args.filea     
    else:   
        # Left, inner or outer join
        fa = args.filea
        fb = args.fileb 

    # Identify key in both files
    ha = fa.readline().decode("utf-8").rstrip().split('\t')
    try:
        indexA = ha.index(args.key)
    except:
        sys.stderr.write('Key not fount in file '+fa.name)
        sys.exit(1)

    hb = fb.readline().decode("utf-8").rstrip().split('\t')
    try:
        indexB = hb.index(args.key)
    except:
        sys.stderr.write('Key not fount in file '+fb.name)
        sys.exit(1)

    # Write header in the output file
    header = ha[:]
    del hb[indexB]
    header.extend(hb)
    args.out.write('{}\n'.format(sep.join(header)))

    # Read the second file in memory and index key
    vIndex = {}
    for line in fb:
        linelist = line.rstrip().decode("utf-8").split('\t')
        lineraw = []
        for i in range(len(linelist)):
            value = linelist[i]
            if i!=indexB:
                lineraw.append(value)
            else:
                if args.soft: key = value[:-3]
                else: key = value
        vIndex[key] = '{}'.format(sep.join(lineraw))
    fb.close()
    
    # Read the first file
    for line in fa:
        line = line.decode("utf-8").rstrip()
        linelist = line.split('\t')
        k = linelist[indexA]
        if args.soft: k = k[:-3]
        if k not in vIndex:
            if args.type != 'inner':
                args.out.write('{}\n'.format(line))
            continue
        args.out.write('{}\n'.format(sep.join([line, vIndex[k]])))
        del vIndex[k]        
    fa.close()

    if args.type == 'outer':
        for key in vIndex:
            # For lines in B that weren't found in A fill in the fields corresponding to the first line's header.
            line = '{}'.format(sep.join(['' if i != indexA else key for i in range(len(ha))]))
            args.out.write('{}\n'.format(sep.join([line, vIndex[key]])))

    args.out.close()

def main ():

    parser = argparse.ArgumentParser(description='Joins the two input files using the column label indicated by the --id parameter as a key. The --soft parameter is used when InChiKey based comparisons are performed, discarding the last 3 chars. By default it performs a left join, but you can also chose right, inner, or outer join.')
    parser.add_argument('-a', '--filea', type=argparse.FileType('rb'), help='First file to join.', required=True)
    parser.add_argument('-b', '--fileb', type=argparse.FileType('rb'), help='Second file to join.', required=True)
    parser.add_argument('-f', '--field', type=str, dest='key', help='Name of the field to be used as a common key.', required=True)
    parser.add_argument('-s', '--soft', action='store_true', help='When InChiKey based comparisons are performed, discard the last 3 chars.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-l', '--left', action='store_const', dest='type', const='left', default='left', help='Left join (default).')
    group.add_argument('-r', '--right', action='store_const', dest='type', const='right', help='Right join.')
    group.add_argument('-i', '--inner', action='store_const', dest='type', const='inner', help='Inner join.')
    group.add_argument('-x', '--outer', action='store_const', dest='type', const='outer', help='Outer join.')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default='output.tsv', help='Output file name (default: output.tsv)')
    args = parser.parse_args()

    join (args)

if __name__ == '__main__':    
    main()
