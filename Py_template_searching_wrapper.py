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
parser.add_argument('--seqid_cutoff', '-c', help='The sequence identity cutoff to use in the search. Value must be between 0.0 and 1.0 (default = 0.9)', default=0.9, type=float)
parser.add_argument('--max_templates', '-max_tmp', help='The maximum number of templates to download. This value is in place to prevent large sets of template structures from being downloaded. The top --max_templates PDB files with the best alignment scores will be downloaded. Input -1 to download all recovered templates from the PDB (default = 100)', default=100, type=int)
parser.add_argument('--cluster_radius', '-r', help='The radius to use for clustering pharmacophore atoms in template ligands, in Angstroms (default = 1.0)', default=1.0, type=float)
parser.add_argument('--min_cluster_size', '-mc', help='The minimum number of atoms within a cluster of pharmacophore atoms (default = 3)', default=3, type=int)
parser.add_argument('--affinity_ligs_only', '-alo', help='Enable this flag to also extract a pharmacophore from ligands with experimental binding affinities reported in the RCSB PDB', action='store_true', default=False)

args = parser.parse_args()

def main():
    os.mkdir(args.outdir) # Throw an error if this directory exists

    cmd1 = f'python3 {BIN}/01_Py_recover_templates_from_the_pdb.py -s {args.protein_sequence} -ic {args.seqid_cutoff} -o {args.outdir}/pdb_templates.json'
    subprocess.run(cmd1.split())

    # Download all templates if max_template is not given
    if args.max_templates == -1:
        with open(f'{args.outdir}/pdb_templates.json') as fo:
            template_data = json.load(fo)
            n_templates = len(template_data)
            cmd2 = f'python3 {BIN}/02_Py_download_template_structures.py -i {args.outdir}/pdb_templates.json -o {args.outdir}/template_structures -m {n_templates}'
    else:
        cmd2 = f'python3 {BIN}/02_Py_download_template_structures.py -i {args.outdir}/pdb_templates.json -o {args.outdir}/template_structures -m {args.max_templates}'

    subprocess.run(cmd2.split())

    cmd3 = f'python3 {BIN}/03_Py_cluster_template_ligand_atoms.py -ld={args.outdir}/template_structures/aligned_ligands/ -od={args.outdir}/pharmacophore_extraction/ --min_size {args.min_cluster_size} --cluster_radius {args.cluster_radius}'

    subprocess.run(cmd3.split())

    cmd4 = f'python3 {BIN}/04_Py_visualize_pharmacophore.py -r={args.outdir}/template_structures/aligned_receptors/reference_receptor.pdb -pd={args.outdir}/pharmacophore_extraction/ -od={args.outdir}/'

    subprocess.run(cmd4.split())

    # If the user specifies --affinity_ligs_only, create an additional pharmacophore
    # from ligands with experimental affinity data
    if args.affinity_ligs_only == True:
        cmd_u1 = f'python3 {BIN}/util01_Py_get_affinity_ligs.py -j {args.outdir}/pdb_templates.json -o {args.outdir}/affinity_ligs.json'
        subprocess.run(cmd_u1.split())

        cmd_u2 = f'python3 {BIN}/util02_Py_compile_affinity_ligs.py -j {args.outdir}/affinity_ligs.json -ld {args.outdir}/template_structures/aligned_ligands/ -od {args.outdir}/template_structures/aligned_affinity_ligands/'
        subprocess.run(cmd_u2.split())
    
        cmd_u3 = f'python3 {BIN}/03_Py_cluster_template_ligand_atoms.py -ld={args.outdir}/template_structures/aligned_affinity_ligands/ -od={args.outdir}/affinity_ligs_pharmacophore_extraction/ --min_size {args.min_cluster_size} --cluster_radius {args.cluster_radius}'
        subprocess.run(cmd_u3.split())
        
        cmd_u4 = f'python3 {BIN}/04_Py_visualize_pharmacophore.py -r={args.outdir}/template_structures/aligned_receptors/reference_receptor.pdb -pd={args.outdir}/affinity_ligs_pharmacophore_extraction/ -od={args.outdir}/affinity_ligs_pharmacophore_extraction'
        subprocess.run(cmd_u4.split())

if __name__=='__main__':
    main()
