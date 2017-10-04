#!/usr/bin/env python

import os, sys, argparse, re
from rdkit import Chem

def getName(mol, count=1, field=None):
    if field is not None:
        name = mol.GetProp(field)
    else:
        name = mol.GetProp('_Name')
        if name == '':
            name = 'mol%0.8d'%count
    return name
    
def setName(mol, ID):
    mol.SetProp("_Name", ID)
    return mol

def getProperties(mol):
    propD = mol.GetPropsAsDict()
    return propD

def writeProperties(fh, propD):
    for prop in propD:
        fh.write('>  <{}>\n{}\n\n'.format(prop, propD[prop])) 

def writeSDF(mol, fh, propD, ID=None):
    if ID:
        mol = setName(mol, ID)
    fh.write(Chem.MolToMolBlock(mol))
    writeProperties(fh, propD)
    fh.write('$$$$\n')