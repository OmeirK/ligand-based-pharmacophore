# Extract coordinates of ligands in the affinity_json if they
# were downloaded

import os
import json
import shutil
import argparse
from rdkit import Chem

parser = argparse.ArgumentParser()

parser.add_argument('--affinity_json', '-j', help='The path to the .json file specifying ligands with known affinity')
parser.add_argument('--ligand_dir', '-ld', help='Path to the directory containing aligned template ligands')
parser.add_argument('--out_dir', '-od', help='Path to a directory containing output files')

args = parser.parse_args()

def main():
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.affinity_json) as f:
        aff_data = json.load(f)
    
    w = Chem.SDWriter(f'{args.out_dir}/ligs_with_affinity.sdf')
    for lig_mol in os.listdir(args.ligand_dir):
        lig_data = lig_mol.split('.')
        caseid = lig_data[0] + '.' + lig_data[1]
        ligid = lig_data[3]
        
        if caseid in aff_data:
            if ligid in aff_data[caseid]:
                print(lig_mol)
                print(caseid, ligid, aff_data[caseid][ligid])

                aff = aff_data[caseid][ligid]['aff']
                aff_unit = aff_data[caseid][ligid]['unit']
                aff_type = aff_data[caseid][ligid]['type']

                mol = Chem.MolFromMolFile(f'{args.ligand_dir}/{lig_mol}')
                mol.SetProp('affinity', str(aff))
                mol.SetProp('affinity_unit', aff_unit)
                mol.SetProp('affinity_type', aff_type)

                w.write(mol)
                shutil.copy(f'{args.ligand_dir}/{lig_mol}', f'{args.out_dir}/{lig_mol}')



if __name__=='__main__':
    main()
