#!/usr/bin/env python

import os, sys, argparse, re
from rdkit import Chem

def getPropertyValue(mol, prop):
    l = [i for i in mol.GetPropNames()]
    if not cmpdID in l: continue

    print (cmpdID)
    print (mol.GetPropNames()[l.index(cmpdID)])
    sys.exit()
    val = mol.GetProp(mol.GetPropNames()[l.index(cmpdID)])
    
    return val

def getName(mol):
    return mol.GetProp("_Name")

def writeProperties(fh, propD):
    for prop in propD:
        fh.write('>  <'+prop+'>\n'+propD[prop]+'\n\n') 
    
def setName(mol, ID):
    mol.SetProp("_Name", ID)
    return mol

def writeMol(mol, fh, propD, ID=None):
    if ID:
        mol = setName(mol, ID)
    mb = Chem.MolToMolBlock(mol)
    fh.write(mb)
    writeProperties(fh, propD)
    fh.write('$$$$\n')