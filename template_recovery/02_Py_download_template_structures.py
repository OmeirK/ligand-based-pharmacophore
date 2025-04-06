import os
import json
import tqdm
import shutil
import argparse
import requests
import subprocess
from pymol import cmd, stored

parser = argparse.ArgumentParser()

parser.add_argument('--template_json', '-i', help='The output from 01_Py_recover_templates_from_the_pdb.py')
parser.add_argument('--outdir', '-o', help='Path to a directory where files will be outputted. If it does not exist it will be created. (default = ./template_structures)', default='./template_structures')
parser.add_argument('--max_templates', '-m', help='The maximum number of template PDBs to download (default = 100)', default=100, type=int)

args = parser.parse_args()

def compile_download_data(template_data, cutoff):
    dl_data = {}
    for i, result in enumerate(template_data['result_set']):
        if i >= cutoff:
            break

        pdb_ch = result['identifier'].split('.')
        pdb = pdb_ch[0]
        chain = pdb_ch[1]
        
        score = result['score']
        
        if pdb not in dl_data:
            dl_data[pdb] = []

        dl_data[pdb].append(chain)

        # Save the top scoring template to be used as a reference
        # structure
        if i == 0:
            ref_pdb = pdb

            
    
    return dl_data, ref_pdb

# Download cif files
def download_cif_files(dl_data, outdir, split_outdir):
    print('Downloading tempaltes from RCSB...')
    for pdb in dl_data:
        pdb_coords = requests.get(f'https://files.rcsb.org/download/{pdb}.cif')
        
        dl_cif = f'{outdir}/{pdb}.cif.gz'
        with open(dl_cif, 'w') as fo:
            fo.write(pdb_coords.text)
        
        pdb_chains = dl_data[pdb]

        # Only save relevant chains
        cmd.reinitialize()
        cmd.load(dl_cif)
        cmd.remove('solvent')
        
        for ch in dl_data[pdb]:
            # Save the receptor and ligand
            cmd.select('pdb_ch', f'({pdb} and segi {ch}) or (byres hetatm within 6 of ({pdb} and segi {ch}))')
            cmd.save(f'{split_outdir}/{pdb}.{ch}.pdb', 'pdb_ch')


# Align all relevant chains to the reference PDB structure
# Collect information for bound ligands
def align_receptors(dl_data, ref_pdb, split_dir, out_dir):
    ref_pdb_ch = dl_data[ref_pdb][0]
    
    cmd.reinitialize()
    cmd.load(f'{split_dir}/{ref_pdb}.{ref_pdb_ch}.pdb', 'ref')
    cmd.save(f'{out_dir}/reference_receptor.pdb', f'ref')
    
    print(f'Aligning all receptors to {split_dir}/{ref_pdb}.{ref_pdb_ch}.pdb')
    # All pdbs are should be loaded in order of their score
    for pdb in dl_data:
        for ch in dl_data[pdb]:
            try:
                cmd.load(f'{split_dir}/{pdb}.{ch}.pdb', f'{pdb}.{ch}')
                cmd.align(f'{pdb}.{ch}', 'ref')
                cmd.save(f'{out_dir}/{pdb}.{ch}_aligned.pdb', f'{pdb}.{ch}')
                cmd.delete(f'{pdb}.{ch}')
            except:
                print(f'Error aligning {pdb}.{ch}! Skipping it')
                continue
        

# Extract ligand coordinates from each aligned model
# Protonate them with OpenBabel so that they can be used in RDKit
# Store results in out_dir
def extract_ligands(dl_data, aln_dir, out_dir):
    print('Attempting to assign bond orders to template ligands...')
    for pdb in tqdm.tqdm(dl_data):
        for ch in dl_data[pdb]:
            cmd.reinitialize()
            cmd.load(f'{aln_dir}/{pdb}.{ch}_aligned.pdb', f'{pdb}_{ch}')
            stored.ligands = []
            cmd.iterate(f'{pdb}_{ch} and hetatm', 'stored.ligands.append((resi,resn,segi))')
            lig_l = list(set(stored.ligands))

            for resi, resn, segi in lig_l:
                #print(f'Convert {pdb}.{ch}.{resn}.{resi}.{segi}.mol')
                cmd.save(f'tmp.mol', f'{pdb}_{ch} and resi {resi} and resn {resn} and segi {segi}')

                obabel_cmd = f'obabel -imol tmp.mol -omol -O{out_dir}/{pdb}.{ch}.{resi}.{resn}.{segi}.mol'
                subprocess.run(obabel_cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    os.remove('tmp.mol')

def main():
    with open(args.template_json) as f:
        template_data = json.load(f)

    dl_data, ref_pdb = compile_download_data(template_data, args.max_templates)

    # Create output directories
    os.makedirs(args.outdir, exist_ok=True)
    dl_dir = f'{args.outdir}/raw_cifs'
    os.makedirs(dl_dir, exist_ok=True)
    split_dir = f'{args.outdir}/split_receptors'
    os.makedirs(split_dir, exist_ok=True)
    
    aln_dir = f'{args.outdir}/aligned_receptors'
    os.makedirs(aln_dir, exist_ok=True)
    
    aln_lig_dir = f'{args.outdir}/aligned_ligands'
    os.makedirs(aln_lig_dir, exist_ok=True)

    download_cif_files(dl_data, dl_dir, split_dir)
    align_receptors(dl_data, ref_pdb, split_dir, aln_dir)
    extract_ligands(dl_data, aln_dir, aln_lig_dir)

    shutil.rmtree(split_dir)
    shutil.rmtree(dl_dir)

if __name__=='__main__':
    main()

