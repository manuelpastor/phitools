#!/usr/bin/env python

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, SaltRemover
import sys

def morgan(mol, fpType, r=4):
    # r is the radius parameter. For instance, r=4 will yield ECFP4 or FCFP4 fingerprints
    if fpType == 'fcfp': feat = True
    else: feat = False

    return AllChem.GetMorganFingerprint(mol,r,useFeatures=feat)

def getFP(mol, fpType, r=4):
    if fpType in ('ecfp', 'fcfp'):
        return morgan(mol, fpType, r)

def getFPdict_smi (fh, molID= None, fpType= 'ecfp', radius= 4, smilesI= 1, header= False):

    if header:
        fh.readline()

    fpD = {}
    count = 0
    for line in fh:
        count += 1
        fields = line.decode('utf-8').rstrip().split('\t')
        if len(fields) > smilesI:
            smiles = fields[smilesI]
        else: continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: continue

        if molID is not None and len(fields) > molID: name = fields[molID]
        else: name = 'mol%0.8d'%count

        fp = getFP(mol, fpType, radius)

        fpD[name] = [fp, smiles]
    
    return fpD

def getFPdict_sdf (fh, molID= None):

    suppl = Chem.ForwardSDMolSupplier(fh,removeHs=False)
    fpD = {}
    count = 0
    for mol in suppl:
        count += 1
        if mol is None: continue
        name = getName(mol, count, molID)
        mol.UpdatePropertyCache(strict=False)
        mh=Chem.AddHs(mol, addCoords=True)
        fpD[name] = AllChem.GetMorganFingerprint(mh,4)
    
    return fpD

def getFPdict (inFormat, fh, molID= None, smilesI= None):

    if inFormat == 'smi':
        if molID is not None:
            try:
                i = int(molID)
            except:
                sys.stderr.write('The ID argument must be a column index if the input file is of smiles format.\n')
                sys.exit()
        else: i = None
        return getFPdict_smi (fh, i, smilesI)
    elif inFormat == 'sdf':
        return getFPdict_sdf (fh, molID)

def RemoveSalts(mol, fh):
    f = open('HighQuality.smi')
    o = open('HighQuality.NoSalts.smi', 'w')
    remover = SaltRemover.SaltRemover()
    for line in f:
        cas, smi = line.rstrip().split('\t')
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = remover.StripMol(mol,dontRemoveEverything=True)
        except:
            pass
        else:
            smi = Chem.MolToSmiles(mol)
        fh.write('{}\t{}\n'.format(cas, smi))
    o.close()

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

def readSmi(fname, smiI, nameI):
    with open(fname) as f:
        suppl = []
        for line in fname:
            fields = line.rstrip().split(sep)
            smi = fields[smiI]
            name = fields[nameI]
            mol = MolFromSmiles(smi)
            setName(mol, name)
    return suppl

def writePropertiesSD(fh, propD):
    for prop in propD:
        fh.write('>  <{}>\n{}\n\n'.format(prop, propD[prop])) 

def writeSDF(mol, fh, propD, ID=None):
    if ID:
        mol = setName(mol, ID)
    fh.write(Chem.MolToMolBlock(mol))
    writeProperties(fh, propD)
    fh.write('$$$$\n')