#!/usr/bin/env python

import os, sys, argparse, re
from rdkit import Chem

def getName(mol):
    return mol.GetProp("_Name")
    
def setName(mol, ID):
    mol.SetProp("_Name", ID)
    return mol

def writeProperties(fh, propD):
    for prop in propD:
        fh.write('>  <{}>\n{}\n\n'.format(prop, propD[prop])) 

def writeSDF(mol, fh, propD, ID=None):
    if ID:
        mol = setName(mol, ID)
    fh.write(Chem.MolToMolBlock(mol))
    writeProperties(fh, propD)
    fh.write('$$$$\n')