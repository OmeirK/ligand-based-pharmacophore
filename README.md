## Introduction

This script will extract template structures depositied in the PDB and generate a visualization of pharmacophore regions where protein-ligand interactions frequently occur. 

The user provides a protein sequence as the input, and can specify a sequence identity cutoff that is used to recover homolgous structures. By default, 100 structures with the best sequence identity are downloaded.

Ligand atoms are clustered on the basis of their pharmacophore atom types, classified using the RDKit definitions.
Ligand atoms are clustered using the approach described in the [E-FTMap publication](https://doi.org/10.1021/acs.jcim.3c01969). The output of this code is formatted to mirror the output of the E-FTMap webserver, however it only clusters ligand atoms from crystal structures deposited in the PDB.

## Installation
Install the dependencies for running the code with:
```
conda env create --file ligand_pharmacophore.yml
```

## Running the code
A wrapper script can be used to run the entire pipeline. To produce the example output, run:
```
python3 Py_template_searching_wrapper.py -p MLLLPLPLLLFLLCSRAEAGEIIGGTESKPHSRPYMAYLEIVTSNGPSKFCGGFLIRRNFVLTAAHCAGRSITVTLGAHNITEEEDTWQKLEVIKQFRHPKYNTSTLHHDIMLLKLKEKASLTLAVGTLPFPSQKNFVPPGRMCRVAGWGRTGVLKPGSDTLQEVKLRLMDPQACSHFRDFDHNLQLCVGNPRKTKSAFKGDSGGPLLCAGVAQGIVSYGRSDAKPPAVFTRISHYRPWINQILQAN -o=chymase_example -c=0.9 -max_tmp=20
```

## Interpreting the output
Within the output directory, pharmacophores can be visualized in pymol using the command:
```
pymol template_pharmacophore_visual.pse
```

The PDB file ```template_pharmacophore.pdb``` can be used to access coordinates of pharmacophore points in a single file. Coordinates for each pharmacophore cluster are written under different headers.

Pharmacophore clusters are ranked on the basis of their size, with larger clusters being assigned a better rank. Clusters are named using the following convention: ```clust.{rank}.{size}.{atom_type}```. A rank of 0 indicates that it is the largest cluster for a given atom type.

Within the visualization, acceptor atoms are shown in red, donor in blue, apolar in yellow, aromatic in organge, and halogen in green. If a cluster contains more atoms, it is represented as a more opaque surface.

Recovered template structures can be viewed in the ```output/template_structures``` directory. All template structures are aligned to a reference protein with the highest sequence identity to the provided sequence
