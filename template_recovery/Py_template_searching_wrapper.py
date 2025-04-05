'''
A wrapper script which runs the template searching pipeline
for a set of ligands that bind to a given protein_sequence
'''

import os
import argparse
import subprocess

BIN='/projectnb/docking/omeir/CASP_16/ligand_assessment/ligand_pipeline_dev/01_template_search'

parser = argparse.ArgumentParser()

parser.add_argument('--protein_sequence' , '-p', help='The sequence of the protein you want to dock a ligand to')
parser.add_argument('--ligand_list', '-l', help='A file containing the name and SMILES strings of ligands that you want to do. Each line should be formatted as LIGAND_NAME SMILES, where the ligand name and SMILES are separated by any white space')
parser.add_argument('--outdir', '-o', help='Directory where search results will be stored')

parser.add_argument('--seqid_cutoff', '-c', help='The sequence identity cutoff to use in the search. Value must be between 0.0 and 1.0 (default = 0.5)', default=0.5, type=float)
parser.add_argument('--max_templates', '-max_tmp', help='The maximum number of templates to download. This value is in place to prevent large sets of template structures from being downloaded. The top --max_templates PDB files with the best alignment scores will be downloaded (default = 100)', default=100, type=int)
parser.add_argument('--min_coverage', '-min_cov', help='The minimum MCS coverage a tempalte ligand must have to be displayed in template ensemble visualizations (default = 0.2)', default=0.2, type=float)

args = parser.parse_args()

def read_ligand_list(lig_list_file):
    ligand_data = []
    with open(lig_list_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l_data = line.strip().split()
                lig_name = l_data[0]
                lig_smi = l_data[1]

                ligand_data.append((lig_name, lig_smi))
    return ligand_data

def main():
    os.mkdir(args.outdir) # Throw an error if this directory exists

    cmd1 = f'python3 {BIN}/01_Py_recover_templates_from_the_pdb.py -s {args.protein_sequence} -ic {args.seqid_cutoff} -o {args.outdir}/pdb_templates.json'
    subprocess.run(cmd1.split())

    cmd2 = f'python3 {BIN}/02_Py_download_template_structures.py -i {args.outdir}/pdb_templates.json -o {args.outdir}/template_structures -m {args.max_templates}'
    subprocess.run(cmd2.split())

    ligand_data = read_ligand_list(args.ligand_list)

    for name, smiles in ligand_data:
        print(f'Finding best templates for {name} ({smiles})...')
        cmd3 = f'python3 {BIN}/03_Py_visualize_best_templates_v2.py -s {smiles} -d {args.outdir}/template_structures -m {args.min_coverage} -o {args.outdir}/{name}_templates'
        subprocess.run(cmd3.split())

if __name__=='__main__':
    main()
