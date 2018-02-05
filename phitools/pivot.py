#!/usr/bin/env python

import sys, os, argparse

sep = '\t'  # Column separator

def pivotmany(args, colsN):
  header = args.fn.readline().rstrip().split(sep)
  freezeI = []
  newLine = []
  for i in range(len(header)):
    if i not in colsN:
      freezeI.append(i)
      newLine.append(header[i])
  newLine.extend(['Category', 'Value'])
  args.out.write('{}\n'.format(sep.join(newLine)))

  for line in args.fn:
    fields = line.rstrip().replace('"', '').split(sep)
    newLine = []
    for i in freezeI:
      try:
        newLine.append(fields[i])
      except:
        newLine.append('')

    for i in colsN:
      field = fields[i].strip()
      if field == '':
        # Empty cell
        continue
      while field[-1] == args.sep:
        field = field.strip(args.sep)
      field = field.split(args.sep)
      for fieldValue in field:
        if fieldValue == '':
          continue
        nl = newLine[:]
        nl.extend([header[i], fieldValue.strip()])
        args.out.write('{}\n'.format(sep.join(nl)))
  args.fn.close()
  args.out.close()

def pivotone(args, coli):
  headerLine = args.fn.readline()
  header = headerLine.rstrip().split(sep)
  freezeI = [i for i in range(len(header)) if i != coli]
  newLine = []
  for i in freezeI:
    newLine.append(header[i])
  newLine.append(header[coli])
  args.out.write('{}\n'.format(sep.join(newLine)))

  for line in args.fn:
    fields = line.rstrip().replace('"', '').split(sep)
    newLine = []
    for i in freezeI:
      try:
        newLine.append(fields[i].strip())
      except:
        newLine.append('')
        continue

    field = fields[coli].strip()
    if field == '':
      # Empty cell
      newLine.append('')
      args.out.write('{}\n'.format(sep.join(newLine)))
      continue
    while field[-1] == args.sep:
      field = field.strip(args.sep)
    
    if field == '':
      # Empty cell
      newLine.append('')
      args.out.write('{}\n'.format(sep.join(newLine)))
      continue

    field = field.split(args.sep)
    for fieldValue in field:
      if fieldValue == '':
        continue
      nl = newLine[:]
      nl.append(fieldValue.strip())
      args.out.write('{}\n'.format(sep.join(nl)))
  args.fn.close()
  args.out.close()

def pivot(args):
  colsN = []
  for i in args.cols.strip().split(','):
    try:
      colsN.append(int(i)-1)
    except:
      iRange = i.split('-')
      try:
        colsN.extend(range(int(iRange[0])-1, int(iRange[1])))
      except:
        sys.stderr.write('Wrong column index format: %s\n' %i)
        sys.exit()

  if len(colsN) == 1:
    pivotone(args, colsN.pop())
  else:
    pivotmany(args, colsN)

def main ():

    parser = argparse.ArgumentParser(description='Pivot the columns indicated in argument -c.')
    parser.add_argument('-f', '--fn', type=argparse.FileType('r'), required=True, help='Data file.')
    parser.add_argument('-c', '--cols', type=str, help='Columns with multile values to be split.')
    parser.add_argument('-s', '--sep', type=str, default= '|', help='Subfield separator (default= \'|\').')
    parser.add_argument('-o', '--out', type=argparse.FileType('w+'), default='output.sdf', help='Output file name (default: output.txt)')
    args = parser.parse_args()
    
    pivot(args)    

if __name__ == '__main__':    
    main()