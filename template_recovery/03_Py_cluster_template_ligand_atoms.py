import os
import glob
import tqdm
import argparse
import numpy as np
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures

parser = argparse.ArgumentParser()

parser.add_argument('--ligand_dir', '-ld', help='Directory containing .mol files of ligand atoms from template structures')

args = parser.parse_args()

# For each mol in a list, get the 3D coordinates
def get_mol_coords(mol):
    conf = mol.GetConformer()
    n_ats = mol.GetNumAtoms()
    coord_arr = np.zeros((n_ats,3))

    for at_id in range(0,n_ats):
        at_pos = conf.GetAtomPosition(at_id)

        coord_arr[at_id][0] = at_pos.x
        coord_arr[at_id][1] = at_pos.y
        coord_arr[at_id][2] = at_pos.z

        #print(at_id, coord_arr[at_id])
    
    
    return coord_arr

# Define chemical features that can be used for each molecule
def detect_pharmacophore_atoms(mols):
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
    pharm_data = {}
    i = 0
    for mol in tqdm.tqdm(mols):
        pharm_data[i] = {'acceptor' : {'idx': [], 'coords': []},
                         'donor' : {'idx': [], 'coords' : []},
                         'apolar' : {'idx': [], 'coords': []},
                         'aromatic' : {'idx': [], 'coords': []},
                         'halogen' : {'idx': [], 'coords': []}}
        mol_coords = get_mol_coords(mol)

        # Get acceptor, donor, and aromatic atoms
        feats = factory.GetFeaturesForMol(mol)
        for j, feat in enumerate(feats):
            feat_atoms = feat.GetAtomIds()
            feat_type = feat.GetFamily()
            
            #print(j, feat_type, feat_atoms)
            if (feat_type == 'Acceptor') or (feat_type == 'NegIonizable'):
                for at_id in feat_atoms:
                    at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
                    if at_id not in pharm_data[i]['acceptor']['idx']: 
                        pharm_data[i]['acceptor']['idx'].append(at_id)
                        pharm_data[i]['acceptor']['coords'].append(mol_coords[at_id])
            elif (feat_type == 'Donor') or (feat_type == 'PosIonizable'):
                for at_id in feat_atoms:
                    at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
                    if at_id not in pharm_data[i]['donor']['idx']: 
                        pharm_data[i]['donor']['idx'].append(at_id)
                        pharm_data[i]['donor']['coords'].append(mol_coords[at_id])
            elif feat_type == 'Aromatic':
                for at_id in feat_atoms:
                    if at_id not in pharm_data[i]['aromatic']['idx']: # Rings can have overlapping atoms
                        pharm_data[i]['aromatic']['idx'].append(at_id)
                        pharm_data[i]['aromatic']['coords'].append(mol_coords[at_id])
            elif (feat_type == 'Hydrophobe') or (feat_type == 'LumpedHydrophobe'):
                for at_id in feat_atoms:
                    if at_id not in pharm_data[i]['apolar']['idx']: 
                        pharm_data[i]['apolar']['idx'].append(at_id)
                        pharm_data[i]['apolar']['coords'].append(mol_coords[at_id])

            if at_elem in ['F', 'Cl', 'Br', 'I']: # Also check for halogen atoms
                pharm_data[i]['halogen']['idx'].append(at_id)
                pharm_data[i]['halogen']['coords'].append(mol_coords[at_id])
            
        print(i, pharm_data[i])
        i += 1

    return pharm_data

def main():
    lig_mols = []
    for mol_f in os.listdir(args.ligand_dir):
        print(mol_f)
        mol = Chem.MolFromMolFile(f'{args.ligand_dir}/{mol_f}')
        if mol is not None:
            lig_mols.append(mol)
    
    lig_pharm_data = detect_pharmacophore_atoms(lig_mols)



if __name__=='__main__':
    main()
