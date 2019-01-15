#!/home/ssericksen/anaconda2/bin/python2.7

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Scaffolds.MurckoScaffold as ms

try:
    raw_smi_file = sys.argv[1]
except:
    print('')
    print('usage:  get_murcko_scaffids.py   raw_smiles_file.smi')
    print('')
    print('note: smiles file should contain lines with smiles')
    print('       strings in field 1 and molecule name in field 2')
    print('       (space-delimited)')
    print('')
    print('')
    exit


with open(raw_smi_file, 'r') as fh:
   data = fh.readlines()

mol_dict = {}
ms_dict = {}

for line in data:
    
    # get the molid
    molid = line.split()[1].strip()

    # get its smiles string
    smi = line.split()[0].strip()

    # get RDKit canonical smiles
    can = Chem.MolToSmiles( Chem.MolFromSmiles(smi) )

    # get canonical murcko smiles from canonical smiles
    murcko = ms.MurckoScaffoldSmiles(can)
    murcko = Chem.MolToSmiles( Chem.MolFromSmiles(murcko) )

    # get generic murcko smiles
    gen_murcko = Chem.MolToSmiles( ms.MakeScaffoldGeneric( Chem.MolFromSmiles(murcko) ) )
    
    # for each molid key, add the smi, murcko mol, and murcko smiles
    mol_dict[molid] = [ smi, can, murcko, gen_murcko ]

    # bin the mols into the different murcko scaffolds observed
    if gen_murcko in ms_dict:
        ms_dict[gen_murcko].append( molid )
    else:
        ms_dict[gen_murcko] = [ molid ]

uniq_scaffs = ms_dict.keys()

print "molid,scaffid"
for molid in mol_dict:
    scaff = mol_dict[molid][3]
    scaffid = uniq_scaffs.index(scaff)
    print('{},{}').format( molid, str(scaffid) ) 



