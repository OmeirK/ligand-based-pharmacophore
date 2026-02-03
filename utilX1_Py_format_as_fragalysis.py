import os
import glob
import shutil
import argparse
import numpy as np
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import rdFMCS

parser = argparse.ArgumentParser()

parser.add_argument('--template_dir', '-td', help='Directory of aligned template structures, created by my ligand-pharmacophore code')
parser.add_argument('--outdir', '-od', help='Path to a directory where output files will be stored')

args = parser.parse_args()

# Check if molecules are enantiomers
# Code adapted from from:
#https://github.com/rdkit/rdkit/discussions/7169
def are_enantiomers(mol1, mol2):
#def are_enantiomers(smi1, smi2):
    #mol1 = Chem.MolFromSmiles(smi1)
    #assert mol1
    #mol2 = Chem.MolFromSmiles(smi2)
    #assert mol2
    can_smi1 = Chem.MolToSmiles(mol1)
    can_smi2 = Chem.MolToSmiles(mol2)
    return ("@" in can_smi1 and "@" in can_smi2
        and can_smi1.replace(
            "@@", "__DOUBLE_AT__"
        ).replace("@", "@@").replace(
            "__DOUBLE_AT__", "@"
        ) == can_smi2)

# Map atoms in the target ligand to atoms in the template
def cs_sym_mappings(target_mol, template_mol, cs_smarts): # accounts for target_mol symmetry
    cs_patt = Chem.MolFromSmarts(cs_smarts)
    target_cs_matches = target_mol.GetSubstructMatches(cs_patt, uniquify=False)
    template_cs_matches = template_mol.GetSubstructMatches(cs_patt, uniquify=False)

    # Debugging #
    #print(target_cs_matches)
    #print(template_cs_matches)
    #print(Chem.MolToSmiles(target_mol), Chem.MolToSmiles(template_mol), cs_smarts)
    #print(template_mol.HasSubstructMatch(cs_patt))
    
    mappings = set()
    for target_cs_match in target_cs_matches:
        for template_cs_match in template_cs_matches:
            mapping = tuple(sorted(zip(target_cs_match, template_cs_match), key=lambda x: x[1]))
            mappings.add(mapping)
    
    mol_sym_matches = target_mol.GetSubstructMatches(target_mol, uniquify=False)
    mappings_reduced = []
    while len(mappings) > 0:
        mapping = list(mappings.pop())
        mappings_reduced.append(mapping)
        redundant_mappings = {tuple((mol_sym_match[i], j) for i, j in mapping) for mol_sym_match in mol_sym_matches}
        mappings -= redundant_mappings
    return mappings_reduced

def is_identical(mol1, mol2):
    can_smi1 = Chem.MolToSmiles(mol1)
    can_smi2 = Chem.MolToSmiles(mol2)

    if can_smi1 == can_smi2:
        return True
    else:
        return False

# Check if molecules physically overlap
def is_overlap(mol1, mol2):
    is_overlap = False
    
    # Check for MCS
    res=rdFMCS.FindMCS([mol2, mol1])
    sym_mappings = cs_sym_mappings(mol1, mol2, res.smartsString)

    mol1_at = mol1.GetConformer()
    mol2_at = mol2.GetConformer()

    for mi, mapping in enumerate(sym_mappings):
        mol1_mcs_pos = []
        mol2_mcs_pos = []
        for i,j in mapping:
            mol1_at = mol1.GetConformer().GetAtomPosition(j)
            mol1_c = [mol1_at.x, mol1_at.y, mol1_at.z]
            mol1_mcs_pos.append(mol1_c)

            mol1_at = mol1.GetConformer().GetAtomPosition(i)
            mol1_c = [mol1_at.x, mol1_at.y, mol1_at.z]
            mol1_mcs_pos.append(mol1_c)
        
        mol1_mcs_pos = np.array(mol1_mcs_pos)
        mol1_mcs_pos = np.array(mol1_mcs_pos)
        
        rmsd = np.sqrt(((mol1_mcs_pos - mol1_mcs_pos)**2).sum(-1).mean())

        if rmsd < 1.0:
            is_overlap = True
            print('\t', 'RMSD', rmsd, is_overlap)
    
    return is_overlap

# If all else fails, just remove hydrogens with pymol. Then
# try and load the updated molecule
def fix_with_pymol(molf):
    cmd.reinitialize()
    cmd.load(molf)
    cmd.remove('elem H')
    cmd.save(molf)

def fix_invalid_mol(molf, tmpdir='tmp_mols/'):
    print(molf)
    os.makedirs(tmpdir, exist_ok=True)
    
    problems = Chem.DetectChemistryProblems(Chem.MolFromMolFile(molf, sanitize=False))
    if problems[0].GetType() == 'AtomValenceException':
        # Parse the error message
        atom_id_l = []
        elem_l = []
        charge_val_l = []
        err_at_l = []
        # Get the problematic atom info
        for p in problems:
            #print(p, p.Message())
            message = p.Message()
            m_data = message.split()
            atom_id = int(m_data[5])
            elem = m_data[6].strip(',')
            charge_val = int(m_data[7].strip(','))

            err_at_l.append((atom_id, elem, charge_val))
        
        # Read the atom info and edit it
        with open(molf) as f:
            mol_lines = f.readlines()
        
        new_mol_lines = []
        curr_atid = 0
        for l in mol_lines:
            l_data = l.split()
            # Check atom lines
            if len(l_data) == 16:
                is_err = False
                for atom_id, elem, charge_val in err_at_l:
                    if curr_atid == atom_id:
                        if l_data[3] == elem:
                            is_err = True
                
                if is_err:
                    #print(elem, l_data[3], l_data[5], charge_val)
                    newl = l[:38] + '0' + l[39:]
                    #print(curr_atid, l.strip())
                    #print(curr_atid, newl.strip(), '< New')
                    new_mol_lines.append(newl)
                else:
                    new_mol_lines.append(l)

                curr_atid += 1

            # Alter the formal charge to be neutral
            elif l.startswith('M  CHG'):
                newl = f'{l_data[0]}  {l_data[1]}  {l_data[2]}  {l_data[3]}  0\n'
                newl = []
                l_data = l.strip().split()
                for i, dat in enumerate(l_data):
                    #print(i, i%2)
                    if i < 2:
                        newl.append(dat)
                    elif (i % 2 == 0) and i >= 3:
                        print('Charge 0')
                        newl.append('0')
                    else:
                        newl.append(dat)

                
                newl.append('\n')
                
                newl = '  '.join(newl)
                #print(l.strip())
                #print(newl.strip(), '< New')
                
                new_mol_lines.append(newl)
            else:
                new_mol_lines.append(l)


        tmp_molf = f'{tmpdir}/{os.path.basename(molf)}'
        #print(tmp_molf)
        with open(tmp_molf, 'w') as fo:
            fo.write(''.join(new_mol_lines))
        
        mol = Chem.MolFromMolFile(tmp_molf)

        if mol is None:
            fix_with_pymol(tmp_molf)
            print('\tPyMOL fix attempt...')
            mol = Chem.MolFromMolFile(tmp_molf)
            
            if mol is None:
                print(f'\tSanitization failed!')
                return None, molf
            else:
                print(f'\tFormal charges fixed!')
                return mol, molf
        else:
            print(f'\tFormal charges fixed!')
            return mol, molf

    else:
        print(f'\tProblem type {problems[0].GetType()} not accounted for in the code :(')
        return None, molf

# Remove ligands with a size below the threshold
def filter_ligs(lig_mol_l, case_id, threshold=4):
    failed_mols = []
    valid_mols = []
    hac_l = []

    enan_data = {}
    alt_data = {}

    for mol_path in lig_mol_l:
        molf = mol_path
        mol = Chem.MolFromMolFile(molf)

        # Attempt a fix
        if mol is None:
            mol, new_molf = fix_invalid_mol(molf)
            molf = new_molf

        if mol is not None:
            #continue #Debug
            ha_count = mol.GetNumHeavyAtoms()
            smi = Chem.MolToSmiles(mol)
            #print(molf, ha_count, smi)
            mol.SetProp('smi', smi)
            mol.SetProp('path', molf)
            
            # Check for potential enantiomers or alt conformations
            for i, hac in enumerate(hac_l):
                if ha_count == hac:
                    is_enan = are_enantiomers(mol, valid_mols[i])
                    is_iden = is_identical(mol, valid_mols[i])
                    print('\t', is_enan, is_iden)

                    if (is_enan) or (is_iden):
                        is_ov = is_overlap(mol, valid_mols[i])

                    # Store info for alternate enantiomers
                    if is_enan and is_ov:
                        enan_molf = valid_mols[i].GetProp('path')
                        
                        if molf not in enan_data:
                            enan_data[molf] = []
                        if enan_molf not in enan_data:
                            enan_data[enan_molf] = []

                        enan_data[molf].append(enan_molf)
                        enan_data[enan_molf].append(molf)
                    
                    # Store info for alternate conformers
                    if is_iden and is_ov:
                        alt_molf = valid_mols[i].GetProp('path')
                        
                        if molf not in alt_data:
                            alt_data[molf] = []
                        if alt_molf not in alt_data:
                            alt_data[alt_molf] = []

                        alt_data[molf].append(alt_molf)
                        alt_data[alt_molf].append(molf)
                        

            valid_mols.append(mol)
            hac_l.append(ha_count)

                
        else:
            print(molf, mol, 'FAILED')
            failed_mols.append(os.path.abspath(molf))
        
    return valid_mols, enan_data, alt_data, failed_mols

def copy_files(copy_rec, copy_mols, pdbid, enan_idx, outdir, exclude_mols=[], alt_mols=[]):
    print(pdbid)
    dirname = f'{pdbid}.e{enan_idx}'

    outpath = f'{outdir}/{dirname}'
    os.makedirs(outpath, exist_ok=True)

    recf_name = os.path.basename(copy_rec)
    new_recf_name = f'{dirname}{recf_name.strip(pdbid)[:-4]}.pdb'

    shutil.copy(copy_rec, f'{outpath}/{new_recf_name}')
    print(f'\t{recf_name} ->', new_recf_name)
    
    alt_idx = 0
    for m in copy_mols:
        molf = m.GetProp('path')

        if molf in exclude_mols:
            continue

        molf_name = os.path.basename(molf)
        if molf in alt_mols:
            print(f'ALTMOL!')
            #newname = f'{dirname}-{molf_name.strip(pdbid)[:-4]}-alt{alt_idx}_ligand.mol'
            newname = f'{dirname}-{molf_name[len(pdbid)+1:]}-alt{alt_idx}_ligand.mol'
            alt_idx += 1
        else:
            #newname = f'{dirname}-{molf_name.strip(pdbid)[:-4]}_ligand.mol'
            newname = f'{dirname}-{molf_name[len(pdbid)+1:]}_ligand.mol'

        shutil.copy(molf, f'{outpath}/{newname}')
        print(f'\t{molf_name} ->', newname)
    pass

def main():
    os.makedirs(args.outdir, exist_ok=True)

    rec_pdbs = glob.glob(f'{args.template_dir}/aligned_receptors/*_aligned.pdb')
    
    all_failed = []
    for pdb_f in rec_pdbs:
        pdbid = os.path.basename(pdb_f).split('_aligned')[0]
        print(pdbid)

        lig_mols = glob.glob(f'{args.template_dir}/aligned_ligands/{pdbid}*.mol')
        
        # No bound ligands found for the pdb.ch
        if len(lig_mols) == 0:
            continue
        
        valid_mols, enan_data, alt_data, failed_mols = filter_ligs(lig_mols, pdbid)
        
        all_failed += failed_mols
        print(len(valid_mols), len(failed_mols))

        
        if (len(enan_data) == 0) and (len(alt_data) == 0):
            #continue #Debug
            copy_files(pdb_f, valid_mols, f'{pdbid}', 0, args.outdir, exclude_mols=[], alt_mols=[])
        elif len(alt_data) > 0:
            #print('alt_data:', alt_data)
            exclude_alt_confs = []
            for alt_molf in alt_data:
                if alt_molf in exclude_alt_confs:
                    continue

                alt_confs = alt_data[alt_molf]

                copy_files(pdb_f, valid_mols, f'{pdbid}', 0, args.outdir, exclude_mols=[], alt_mols=alt_confs)
                exclude_alt_confs += alt_confs
        else:
            #continue #Debug
            for i, enan in enumerate(enan_data):
                copy_files(pdb_f, valid_mols, f'{pdbid}', i, args.outdir, exclude_mols=enan_data[enan], alt_mols=[])
        
    with open(f'{args.outdir}/failed_ligand_processing.txt', 'w') as fo:
        fo.write('\n'.join(all_failed))


if __name__=='__main__':
    main()
