import json
import requests
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--template_json', '-j', help='The path to the .json file with the recovered templates')
parser.add_argument('--output_json', '-o', help='The path to a .json file where output files will be stored (default = affinity_ligs.json)', default='affinity_ligs.json')

args = parser.parse_args()

AFF_PRIORITY = {
'IC50': 0,
'Ki': 1,
'Kd': 2
}

def main():
    print(f'Finding templates in {args.template_json} that bind ligands with known affinity...')
    with open(args.template_json) as f:
        template_data = json.load(f)
    
    # Save ligands that have an affinity deposited in the PDB
    saved_data = {}
    for result in template_data["result_set"]:
        case = result['identifier']
        case_d = case.split('.')
        pdb = case_d[0]
        segi = case_d[1]
        data = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
        results = data.json()
        
        try:
            affinity_data = results['rcsb_binding_affinity']
        except:
            continue
        
        for data in affinity_data:
            ligand_id = data['comp_id']
            aff = data['value']
            aff_type = data['type']
            aff_unit = data['unit']

            # Only save ligands with IC50, Ki, or Kd
            if aff_type not in AFF_PRIORITY:
                continue

            #print(case, pdb, ligand_id, aff, aff_unit, aff_type)

            if case not in saved_data:
                saved_data[case] = {}

            if ligand_id not in saved_data[case]:
                saved_data[case][ligand_id] = {'aff': aff, 'type': aff_type, 'unit': aff_unit}
            else:
                # Check if affinity should be updated
                # Prioritize Kd > Ki > IC50
                priority = AFF_PRIORITY[aff_type]

                if priority > AFF_PRIORITY[saved_data[case][ligand_id]['type']]:
                    #print(f'\tUpdating affinity -> {aff} {aff_unit} {aff_type}')
                    saved_data[case][ligand_id] = {'aff': aff, 'type': aff_type, 'unit': aff_unit}
                
                if (aff_type == saved_data[case][ligand_id]['type']) and (aff < saved_data[case][ligand_id]['aff']):
                    #print(f'\tUpdating affinity -> {aff} {aff_unit} {aff_type}')
                    saved_data[case][ligand_id] = {'aff': aff, 'type': aff_type, 'unit': aff_unit}
                    
    with open(args.output_json, 'w') as fo:
        json.dump(saved_data, fo, indent=4)
    


if __name__=='__main__':
    main()
