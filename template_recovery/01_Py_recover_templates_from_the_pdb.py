import json
import requests
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--sequence', '-s', help='The sequence of the protein used for the template search')
parser.add_argument('--identity_cutoff', '-ic', help='The sequence identity cutoff to use in the search. Value must be between 0.0 and 1.0 (default = 0.5)', default=0.5, type=float)
parser.add_argument('--evalue_cutoff', '-ec', help='The evalue cutoff to use in the search (default = 0.1)', default=0.1, type=float)
parser.add_argument('--outfile', '-o', help='The name of a .json file with the recovered templates (default = pdb_templates.json)', default='pdb_templates.json')

args = parser.parse_args()

# Build a query json for the RCSB API search
# Find all PDBs with sequence ID above identity_cutoff, with
# non polymer ligands bound
query_data = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "sequence_type": "protein",
                    "value": args.sequence,
                    "identity_cutoff": args.identity_cutoff,
                    "evalue_cutoff": args.evalue_cutoff
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
                    "operator": "exists"
                }
            }
        ]
    },
    
    "request_options": {
        "return_all_hits": True
    },
    "return_type": "polymer_instance"
}

my_search = json.dumps(query_data)

data = requests.get(f"https://search.rcsb.org/rcsbsearch/v2/query?json={my_search}")
results = data.json()

print(f'{results["total_count"]} hits were found in the search!')
print(f'Writing output to {args.outfile}...')

with open(args.outfile, 'w') as fo:
    json.dump(results, fo, indent=4)

