import urllib
from rdkit import Chem
from rdkit.Chem import AllChem,Draw,Descriptors
import os
import sys
import getopt
import re
from PIL import ImageTk, Image
import glob

def viewSDF (sdf,out,prop):

    if sdf:
        sdf_id = sdf
    else:
        sdf_id = 'file.sdf'

    if out:
        out_id = out
    else:
        out_id = 'sdf.html'

    if prop:
        prop_id = prop
    else:
        prop_id = 'name'

    # Read SDF
    suppl = Chem.SDMolSupplier(sdf_id)

    # Create directory to save images
    if not os.path.exists('images/'):
        os.makedirs('images/')
        
    # Get property headers
    prop_names=[prop_id]
    for i in suppl[0].GetPropNames():
        prop_names.append(i.lower())

    # Create a list containing the sdf data
    l=[]
    l.append(prop_names)
    for mol in suppl:
        if mol is None: continue
        
        #Save images
        name=mol.GetProp(prop_id)
        Draw.MolToFile(mol,'images/'+name+'.png')
       
        #Store properties
        l1=[]
        l1.append('<img src="images/'+name+'.png" height=200 width=200>')
        for i in prop_names[1:]:
            l1.append(mol.GetProp(mol.GetPropNames()[prop_names.index(i)-1]))

        l.append(l1)
        
    saveListAsHTML(l, out_id)
        
def saveListAsHTML(l,out):

    fo = open (out,'w+')
    fo.write('<table border="0.5">\n')

    for column in l:
        fo.write('<tr>',)    
        for i in column:
            fo.write('<td>'+i+'</td>')      
        fo.write('</tr>\n')
        
    fo.write('</table>')
    fo.close()

def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'viewSDF -f file.sdf -o sdf.html -p name' 

def main ():
    sdf = None
    out = None
    prop = None
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'f:o:p:', [])
    except getopt.GetoptError:
       usage()
       print "False, Arguments not recognized"
       sys.exit(1)

    if args:
       usage()
       print "False, Arguments not recognized"

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-f':
                sdf = arg
            elif opt in '-o':
                out = arg
            elif opt in '-p':
                prop = arg

    viewSDF(sdf,out,prop)

if __name__ == '__main__':    
    main()
