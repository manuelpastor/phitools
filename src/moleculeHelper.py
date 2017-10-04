#!/usr/bin/env python

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, SaltRemover
from standardise import standardise

from pathlib import Path
import sys, tempfile

rand_str = lambda n: ''.join([random.choice(string.ascii_lowercase) for i in range(n)])

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

    name = ''

    if field is not None and mi.HasProp (field):
        name = mol.GetProp(field)
    else:
        name = mol.GetProp('_Name')
        
    if name == '':
        name = 'mol%0.8d'%count

    # get rid of strange characters
    name = name.decode('utf-8')
    name = name.encode('ascii','ignore')  # use 'replace' to insert '?'

    if ' ' in name:
        name = name.replace(' ','_')

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

def writeSDF(mol, fh, propD=None, ID=None):
    if ID:
        mol = setName(mol, ID)
    fh.write(Chem.MolToMolBlock(mol))
    if propD == None:
        getProperties(mol)
    writeProperties(fh, propD)
    fh.write('$$$$\n')

def standardize(mol):
    """
    Wrapper to aply the structure normalization protocol provided by Francis Atkinson (EBI). If no non-salt components can be found in the mixture, the original mixture is returned.

    Returns a tuple containing:
        1) True/False: depending on the result of the method
        2) (if True ) The output molecule
           (if False) The error message
    """
    try:
        parent = standardise.run(Chem.MolToMolBlock(mol))
    except standardise.StandardiseException as e:
        if e.name == "no_non_salt":
            parent = Chem.MolToMolBlock(mol)
        else:
            return (False, e.name)

    return (True, parent)

def protonateMol(mol, pH= 7.4, mokaPath= os.environ.get('MOKA_PATH'), clean= True):
    """Adjusts the ionization state of the molecule.

        In this implementation, it uses blabber_sd from Molecular Discovery
        The result is a tuple containing:
           1) True/False: describes the success of the protonation for this compound
           2) (if True ) The name of the protonated molecules and its formal charge
              (if False) The error message
    """
    stderrf = open (os.devnull, 'w')
    stdoutf = open (os.devnull, 'w')

    if mokaPath == None:
        return (False, 'Moka path not found.', None)

    name = getName(mol)
    fnameI = '.'.join([name, 'sdf'])
    if not fnameI.is_file():
        with open(fnameI, 'w') as f:
            writeSDF(mol, f)

    fnameO = '.'.join([name, 'prot', 'sdf'])
    call = [mokaPath+'/blabber_sd', fnameI,
            '-p',  str(pH),
            '-o',  fnameO]

    try:
        retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
    except:
        return (False, 'Blabber execution error', 0.0)

    stdoutf.close()
    stderrf.close()

    if 'blabber110' in mokaPath: # in old blabber versions, error is reported as '0'
        if retcode == 0:
            return (False, 'Blabber 1.0 execution error', 0.0)
    else:
        if retcode != 0:
            return (False, 'Blabber execution error', 0.0)

    try:
        if os.stat(fnameO).st_size==0:
            return (False, 'Blabber output is empty', 0.0)

        finp = open (fnameO)
    except:
        return (False, 'Blabber output not found', 0.0)

    charge = 0
    for line in finp:
        if line.startswith ('M  CHG'):
            items = line.split()
            if int(items[2]):
                for c in range (4,len(items),2): charge+=int(items[c])
            break
    finp.close()

    suppl = SDMolSupplier(fnameO)
    molO = suppl.next()

    if clean:
        removefile (fnameI)
        removefile (fnameO)

    return (True, molO, charge)

def protonateFile(fnameI, pH= 7.4, mokaPath= os.environ.get('MOKA_PATH'), clean= True):
    """Adjusts the ionization state of the molecule.

        In this implementation, it uses blabber_sd from Molecular Discovery
        The result is a tuple containing:
           1) True/False: describes the success of the protonation for this compound
           2) (if True ) The name of the protonated molecules and its formal charge
              (if False) The error message
    """
    stderrf = open (os.devnull, 'w')
    stdoutf = open (os.devnull, 'w')

    pre, ext = os.path.splitext(fnameI)
    fnameO = pre + '.protonated' + ext

    mokaRun = mokaPath+'/blabber_sd'
    if mokaPath == None or  not mokaRun.is_file():
        return (False, 'Moka path not found.')
    call = [mokaRun, fnameI, 
            '-p',  str(pH),
            '-o',  fnameO]

    try:
        retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
    except:
        return (False, 'Blabber execution error')

    stdoutf.close()
    stderrf.close()

    if 'blabber110' in mokaPath: # in old blabber versions, error is reported as '0'
        if retcode == 0:
            return (False, 'Blabber 1.0 execution error')
    else:
        if retcode != 0:
            return (False, 'Blabber execution error')

    try:
        if os.stat(fnameO).st_size==0:
            return (False, 'Blabber output is empty')
    except:
        return (False, 'Blabber output not found')

    return (True, fnameO)

def convert3D(mol, corinaPath= os.environ.get('CORINA_PATH'), clean= True):
    """Converts the 2D structure of the molecule "moli" to 3D

        In this implementation, it uses CORINA from Molecular Networks
        The result is a tuple containing:
           1) suucTrue/False: describes the success of the 3D conversion for this compound
           2) (if True ) The name of the 3D molecule
              (if False) The error mesage
    """

    stderrf = open (os.devnull, 'w')
    stdoutf = open (os.devnull, 'w')

    name = getName(mol)
    fnameI = '.'.join([name, 'sdf'])
    if not fnameI.is_file():
        with open(fnameI, 'w') as f:
            writeSDF(mol, f)

    fnameO = '.'.join([name, '3D', 'sdf'])
    call = [corinaPath+'/corina',
            '-dwh','-dori',
            '-ttracefile=corina.trc',
            '-it=sdf', fnameI,
            '-ot=sdf', fnameO]

    try:
        retcode = subprocess.call(call, stdout=stdoutf, stderr=stderrf)
    except:
        return (False, 'Corina execution error')

    stdoutf.close()
    stderrf.close()

    if retcode != 0:
        return (False, 'Corina execution error')

    if not os.path.exists(fnameO):
        return (False, 'Corina output not found')

    suppl = SDMolSupplier(fnameO)
    molO = suppl.next()

    if clean:
        removefile(fnameI)
        removefile(fnameO)
        removefile('corina.trc')

    return (True, molO)