'''
A wrapper script which runs the template searching pipeline
for a set of ligands that bind to a given protein_sequence
'''

import os
import json
import argparse
import subprocess

BIN='./template_recovery/'

parser = argparse.ArgumentParser()

parser.add_argument('--protein_sequence' , '-p', help='The sequence of the protein you want to dock a ligand to')
parser.add_argument('--outdir', '-o', help='Directory where search results will be stored')

parser.add_argument('--seqid_cutoff', '-c', help='The sequence identity cutoff to use in the search. Value must be between 0.0 and 1.0 (default = 0.0)', default=0.9, type=float)
parser.add_argument('--max_templates', '-max_tmp', help='The maximum number of templates to download. This value is in place to prevent large sets of template structures from being downloaded. The top --max_templates PDB files with the best alignment scores will be downloaded. Input -1 to download all recovered templates from the PDB (default = 100)', default=100, type=int)

args = parser.parse_args()

def main():
    os.mkdir(args.outdir) # Throw an error if this directory exists

    cmd1 = f'python3 {BIN}/01_Py_recover_templates_from_the_pdb.py -s {args.protein_sequence} -ic {args.seqid_cutoff} -o {args.outdir}/pdb_templates.json'
    subprocess.run(cmd1.split())

    # Download all tempaltes if max_template is not given
    if args.max_templates == -1:
        with open(f'{args.outdir}/pdb_templates.json') as fo:
            template_data = json.load(fo)
            n_templates = len(template_data)
            cmd2 = f'python3 {BIN}/02_Py_download_template_structures.py -i {args.outdir}/pdb_templates.json -o {args.outdir}/template_structures -m {n_templates}'
    else:
        cmd2 = f'python3 {BIN}/02_Py_download_template_structures.py -i {args.outdir}/pdb_templates.json -o {args.outdir}/template_structures -m {args.max_templates}'

    subprocess.run(cmd2.split())

if __name__=='__main__':
    main()
