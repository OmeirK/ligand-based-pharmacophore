import os
import json
import glob
import argparse
from Bio import SeqIO
from rdkit import Chem

parser = argparse.ArgumentParser()

parser.add_argument('--rec_fa', '-r', help='fasta file with the receptor sequences')
parser.add_argument('--fragalysis_dir', '-fd', help='Path to the directory containing fragalysis results')
parser.add_argument('--msa_dir', '-msa', help='Directory containing MSA info, in the OF3 format. NOTE: Only tested for single chain proteins! (default = None)', default=None)
parser.add_argument('--outfile', '-o', help='Name of the output .json file containing of3 inputs')

args = parser.parse_args()

ALPHA='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def main():
    fasta_sequences = SeqIO.parse(open(args.rec_fa),'fasta')
    print(fasta_sequences)
    
    rec_seq_l = []
    A_i = 0
    for i, fa in enumerate(fasta_sequences):
        print(str(fa.seq))

        ch=ALPHA[A_i]
        
        rec_seq_l.append((ch, str(fa.seq)))

        A_i += 1

    print(rec_seq_l)

    # Read ligand info
    
    lig_data = {}
    case_dirs = glob.glob(f'{args.fragalysis_dir}/*/')

    for i, cd in enumerate(case_dirs):
        case = os.path.basename(cd[:-1])
        sdf_l = glob.glob(f'{cd}/{case}*_ligand.sdf')

        lig_data[case] = []

        for sdf in sdf_l:
            if '-alt' in os.path.basename(sdf):
                continue
            else:
                suppl = Chem.SDMolSupplier(sdf)
                m = suppl[0]
                print(sdf)
                print(m)
                print(m.GetPropsAsDict())
                smi = m.GetProp('smi')
                lig_data[case].append(smi)
            


    #for i, case in enumerate(df['Case_ID']):
    #    smi = df['Mol_SMILES'].iloc[i]
    #
    #    if case not in lig_data:
    #        lig_data[case] = []
    #
    #    if smi not in lig_data[case]:
    #        lig_data[case].append(smi)


    of3_inps = {'seeds': [1370180479], 'num_seeds': 1, 'queries': {}}
    for c in lig_data:
        of3_inps['queries'][c] = {'chains': []}
        for r_ch, r_seq in rec_seq_l:
            
            if args.msa_dir == None:
                of3_inps['queries'][c]['chains'].append({'molecule_type': 'protein',
                                                          'chain_ids': r_ch,
                                                          'sequence': r_seq
                                                         }
                                                   )

            else:
                of3_inps['queries'][c]['chains'].append({'molecule_type': 'protein',
                                                          'chain_ids': r_ch,
                                                          'sequence': r_seq,
                                                          'main_msa_file_paths': os.path.abspath(args.msa_dir)
                                                         }
                                                   )
        for i, smi in enumerate(lig_data[c]):
            lig_ch_i = A_i + i
            lig_ch = ALPHA[lig_ch_i]
            of3_inps['queries'][c]['chains'].append({'molecule_type': 'ligand',
                                                      'chain_ids': lig_ch,
                                                      'smiles': smi
                                                     }
                                                   )


        print(c, of3_inps['queries'][c])


    with open(args.outfile, 'w') as fo:
        json.dump(of3_inps, fo, indent=4)


        

            
if __name__=='__main__':
    main()
