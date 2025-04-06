## Introduction

This script will extract template structures depositied in the PDB and generate a visualization of pharmacophore regions where protein-ligand interactions frequently occur. 

The user provides a protein sequence as the input, and can specify a sequence identity cutoff that is used to recover homolgous structures. By default, 100 structures with the best sequence identity are downloaded.

Ligand atoms are clustered on the basis of their pharmacophore atom types, classified using the RDKit definitions.
Ligand atoms are clustered using the approach described in the [E-FTMap publication](https://doi.org/10.1021/acs.jcim.3c01969).

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
pymol pharmacophore_visual.pse
```
Pharmacophore clusters are ranked on the basis of their size, with larger clusters being assigned a better rank. Clusters are named using the following convention: ``` clust.{rank}.{size}.{atom_type}```
