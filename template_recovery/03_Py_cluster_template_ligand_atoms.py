import os
import glob
import tqdm
import argparse
import numpy as np
from pymol import cmd, stored
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures

parser = argparse.ArgumentParser()

parser.add_argument('--ligand_dir', '-ld', help='Directory containing .mol files of ligand atoms from template structures')
parser.add_argument('--outdir', '-od', help='Directory where .xyz files of pharmacophore atoms will be stored (default = ./)', default='./')
parser.add_argument('--cluster_radius', '-r', help='The radius to use when clustering pharmacophore atoms, in Angstroms (default = 2.0)', type=float, default=2.0)
parser.add_argument('--min_size', '-m', help='The minumum number of atoms that must be retained in a cluster (default = 3)', type=int, default=3)

ATOM_TYPES = ['acceptor', 'donor', 'apolar', 'aromatic', 'halogen']

# Define the colors for each atom type to use in the visualization
ATYPE_COLORS = {'polar' : 'violetpurple',
                'donor' : 'lightblue',
                'acceptor' : 'deepsalmon',
                'halogen' : 'palegreen',
                'aromatic' : 'brightorange',
                'apolar' : 'paleyellow'
}

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
# Fit RDKit featues to match E-FTMap atom types
def detect_pharmacophore_atoms(mols):
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
    pharm_data = {}
    i = 0
    for mol in tqdm.tqdm(mols):
        pharm_data[i] = {'acceptor' : {'idx': [], 'coords': [], 'at_elem': []},
                         'donor' : {'idx': [], 'coords' : [], 'at_elem': []},
                         'apolar' : {'idx': [], 'coords': [], 'at_elem': []},
                         'aromatic' : {'idx': [], 'coords': [], 'at_elem': []},
                         'halogen' : {'idx': [], 'coords': [], 'at_elem': []}}
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
                        pharm_data[i]['acceptor']['at_elem'].append(at_elem)
                        pharm_data[i]['acceptor']['coords'].append(mol_coords[at_id])
            elif (feat_type == 'Donor') or (feat_type == 'PosIonizable'):
                for at_id in feat_atoms:
                    at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
                    if at_id not in pharm_data[i]['donor']['idx']: 
                        pharm_data[i]['donor']['idx'].append(at_id)
                        pharm_data[i]['donor']['at_elem'].append(at_elem)
                        pharm_data[i]['donor']['coords'].append(mol_coords[at_id])
            elif feat_type == 'Aromatic':
                for at_id in feat_atoms:
                    at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
                    if at_id not in pharm_data[i]['aromatic']['idx']: # Rings can have overlapping atoms
                        pharm_data[i]['aromatic']['idx'].append(at_id)
                        pharm_data[i]['aromatic']['at_elem'].append(at_elem)
                        pharm_data[i]['aromatic']['coords'].append(mol_coords[at_id])
            elif (feat_type == 'Hydrophobe') or (feat_type == 'LumpedHydrophobe'):
                for at_id in feat_atoms:
                    at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
                    if at_id not in pharm_data[i]['apolar']['idx']: 
                        pharm_data[i]['apolar']['idx'].append(at_id)
                        pharm_data[i]['apolar']['at_elem'].append(at_elem)
                        pharm_data[i]['apolar']['coords'].append(mol_coords[at_id])
        
        # Check for halogen atoms
        for at_id, atom in enumerate(mol.GetAtoms()):
            at_elem = mol.GetAtomWithIdx(at_id).GetSymbol()
            if at_elem in ['F', 'Cl', 'Br', 'I']: # Also check for halogen atoms
                pharm_data[i]['halogen']['idx'].append(at_id)
                pharm_data[i]['halogen']['at_elem'].append(at_elem)
                pharm_data[i]['halogen']['coords'].append(mol_coords[at_id])
            
        #print(i, pharm_data[i])
        i += 1

    return pharm_data

# Save pharmacophore atom coordinates
def save_pharmacophore_atoms(pharm_data, outdir):
    for atype in ATOM_TYPES:
        outlines = []
        for i in pharm_data:
            for j, elem in enumerate(pharm_data[i][atype]['at_elem']):
                coords = pharm_data[i][atype]['coords'][j]
                outlines.append(f'{elem} {coords[0]} {coords[1]} {coords[2]}\n')


        with open(f'{outdir}/pharm_{atype}.xyz', 'w') as fo:
            fo.write(f'{len(outlines)}\n')
            fo.write(f'{atype}\n')
            fo.write(''.join(outlines))

# Save clusters of pharmacophore atoms extacted
# from xyz files 
def create_xyz_clusters(outdir, radius=2.0, min_size=3):
    for xyz_f in glob.glob(f'{outdir}/pharm_*.xyz'):
        atype = os.path.basename(xyz_f).split('_')[1][:-4]

        os.makedirs(f'{outdir}/clusts_{atype}', exist_ok=True)

        cmd.reinitialize()
        cmd.load(xyz_f, f'pharm_{atype}')
        
        print(f'\tClustering atoms in {xyz_f}...')
        clust_lig_atoms(f'pharm_{atype}', radius, atype)
        
        for obj in cmd.get_object_list(f'clust.*.{atype}'):
            if cmd.count_atoms(obj) >= min_size:
                cmd.save(f'{outdir}/clusts_{atype}/{obj}.xyz', obj)

# Cluster all atoms within a given {pymol_obj} 
# Clusters are iteratively grown to include members within
# a radius defined by distance {radius}
def clust_lig_atoms(pymol_obj, radius, clust_prefix):
    stored.all_ats = []
    cmd.iterate(pymol_obj, 'stored.all_ats.append(index)')
    
    # For debugging purposes:
    tot_ats = cmd.count_atoms(pymol_obj)
    grouped_ats = 0

    skip_ats = [] # Atoms that have already been clustered
    out_clusts = []
    clust_idx = 0
    for i, at_idx in enumerate(stored.all_ats):
        if at_idx in skip_ats:
            continue
        
        cmd.select('tmp', f'{pymol_obj} and index {at_idx}')
        old_size = 1 # Single atom
        cmd.select('tmp', f'{pymol_obj} within {radius} of tmp') # Expand the cluster once
        curr_size = cmd.count_atoms('tmp')

        # Expand a cluster until it can no longer be expanded
        while curr_size != old_size:
            cmd.select('tmp',  f'{pymol_obj} within {radius} of tmp')
            curr_size2 = cmd.count_atoms('tmp')

            old_size = curr_size
            curr_size = curr_size2

        stored.clust = []
        cmd.iterate('tmp', 'stored.clust.append(index)')
        skip_ats += stored.clust

        cmd.create(f'grpclust.{clust_idx}.{curr_size}.{clust_prefix}', 'tmp')

        out_clusts.append((f'grpclust.{clust_idx}.{curr_size}.{clust_prefix}', curr_size))

        grouped_ats += curr_size

        clust_idx += 1

    #print(grouped_ats, tot_ats) # Check that all atoms were clustered
    if grouped_ats != tot_ats:
        print(f'WARNING: Not all atoms were clustered!! Something went wrong!')

    # Sort clusters by size, and rename them
    out_clusts.sort(key=lambda x: x[1], reverse=True)
    for i, vals in enumerate(out_clusts):
        clust = vals[0]
        size = vals[1]
        cmd.set_name(clust, f'clust.{i}.{size}.{clust_prefix}')
        
def main():
    os.makedirs(args.outdir, exist_ok=True)


    lig_mols = []
    for mol_f in os.listdir(args.ligand_dir):
        mol = Chem.MolFromMolFile(f'{args.ligand_dir}/{mol_f}')
        if mol is not None:
            lig_mols.append(mol)
    
    print(f'Detecting pharmacophore atoms...')
    lig_pharm_data = detect_pharmacophore_atoms(lig_mols)
    
    print(f'Saving pharmacophore cooridnates to: {args.outdir}')
    save_pharmacophore_atoms(lig_pharm_data, args.outdir)

    print(f'Saving clusters of pharmacophore atoms...')
    create_xyz_clusters(args.outdir, radius=args.cluster_radius, min_size=args.min_size)


if __name__=='__main__':
    main()
