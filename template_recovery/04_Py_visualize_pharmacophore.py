'''
Create a .pse file with pharmacophore atoms loaded onto a reference receptor
'''

import os
import glob
import argparse
from pymol import cmd

parser = argparse.ArgumentParser()
parser.add_argument('--ref_receptor', '-r', help='A .pdb file of a reference receptor which will be used for the visualization')
parser.add_argument('--pharm_dir', '-pd', help='Path to the pharmacophore_extraction/ directory with clustered pharmacophore atoms')
parser.add_argument('--out_pse', '-o', help='Name of the output .pse file (default = visual.pse)', default='visual.pse')

args = parser.parse_args()

# Dictionary that will color each pharm type on a spectrum from their main color -> grey
COLOR_DICT = {'donor' : {0: [0, 0, 255], 1: [86, 86, 255], 2: [138, 138, 255], 3: [200, 200, 255], 4: [220, 220, 220]},
              'acceptor' : {0: [255, 0, 0], 1: [255, 86, 86], 2: [255, 138, 138], 3: [255, 200, 200], 4: [220, 220, 220]},
              'apolar' : {0: [255, 255, 0], 1: [255, 255, 86], 2: [255, 255, 138], 3: [255, 255, 226], 4: [220, 220, 220]},
              'aromatic': {0 : [255, 138, 0], 1: [255, 178, 86], 2: [255, 200, 138], 3: [255, 229, 198], 4: [220, 220, 220]},
              'halogen' : {0: [0, 255, 0], 1: [86, 255, 86], 2: [138, 255, 138], 3: [200, 255, 200], 4: [220, 220, 200]}}

# Make clusters with less atoms more transparent
def density_clust_by_size(pharm_type, clust_n='clust', obj_l=[]):
   if len(COLOR_DICT[pharm_type]) == 0:
      return 

   newcolor = COLOR_DICT[pharm_type][0]
   color_n = f'{pharm_type}_0'

   if obj_l == []:
      obj_l = cmd.get_object_list(f'{clust_n}*.{pharm_type}')

   max_rank = -1
   for obj in obj_l:
      rank = int(obj.split('.')[2])
      if rank > max_rank:
         max_rank = rank

   for obj in obj_l:
      size = int(obj.split('.')[2])
      rank_pct = size/max_rank

      transparency = 1.0 - (0.6*rank_pct)
      if transparency > 0.8:
        transparency = 0.8

      #print(obj, rank_pct, transparency)

      cmd.set_color(color_n, newcolor)
      cmd.color(color_n, obj)
      cmd.set('transparency', transparency, obj)

def color_clust_by_rank(pharm_type, clust_n='ligclust', obj_l=[]):
   if len(COLOR_DICT[pharm_type]) == 0:
      return 

   if obj_l == []:
      obj_l = cmd.get_object_list(f'{clust_n}*.{pharm_type}')

   max_rank = -1
   for obj in obj_l:
      rank = int(obj.split('.')[2])
      if rank > max_rank:
         max_rank = rank
   
   print(f'\tColoring {pharm_type}')
   #print(max_rank)
   #print(COLOR_DICT[pharm_type])
   
   for obj in obj_l:
      #print(f'Coloring {obj}')
      size = int(obj.split('.')[2])
      rank_pct = size/max_rank

      if size == 1:
         #print('color 4', obj, rank_pct)
         color_n = f'{pharm_type}_4'
         newcolor = COLOR_DICT[pharm_type][4]
      elif rank_pct >= 0.9:
         #print('color 0', obj, rank_pct)
         color_n = f'{pharm_type}_0'
         newcolor = COLOR_DICT[pharm_type][0]
      elif (rank_pct >= 0.5) and (rank_pct < 0.9):
         #print('color 1', obj, rank_pct)
         color_n = f'{pharm_type}_1'
         newcolor = COLOR_DICT[pharm_type][1]
      elif (rank_pct >= 0.3) and (rank_pct < 0.5):
         #print('color 2', obj, rank_pct)
         color_n = f'{pharm_type}_2'
         newcolor = COLOR_DICT[pharm_type][2]
      elif (rank_pct >= 0.1) and (rank_pct < 0.3):
         #print('color 3', obj, rank_pct)
         color_n = f'{pharm_type}_3'
         newcolor = COLOR_DICT[pharm_type][3]
      elif ((rank_pct >= 0.0) and (rank_pct < 0.2)):
         #print('color 4', obj, rank_pct)
         color_n = f'{pharm_type}_4'
         newcolor = COLOR_DICT[pharm_type][4]
      
      cmd.set_color(color_n, newcolor)
      cmd.color(color_n, obj)

def main():
    # Get cluster files
    cmd.reinitialize()
    cmd.load(args.ref_receptor, 'receptor')

    clust_dirs = glob.glob(f'{args.pharm_dir}/clusts_*/')
    print(clust_dirs)
    
    # Load cluster atoms for each atom type, ordered by rank
    atype_l = []
    for cd in clust_dirs:
        atype = os.path.basename(cd[:-1]).split('_')[1]
        atype_l.append(atype)
        print(cd, atype)
        
        xyz_inf = []
        for xyz_f in os.listdir(cd):
            xyz_data = xyz_f.split('.')
            rank = int(xyz_data[1])
            
            xyz_inf.append((xyz_f, rank))
        
        xyz_inf.sort(key=lambda x: x[1])
        
        # Load files in order
        for xyz_f, rank in xyz_inf:
            cmd.load(f'{cd}/{xyz_f}')

    # Set visualization settings
    cmd.disable('*') # Speeds the code up (anecdotally)
    cmd.bg_color('white')
    cmd.set('solvent_radius', 0.2, 'clust.*')
    cmd.set('surface_quality', 1, 'clust.*')
    cmd.set('transparency', 0.5, 'clust.*')
    cmd.hide('everything', f'clust.*')
    
    # Create a gaussian surface for each cluster
    # This will alter surface representations
    cmd.enable('clust.*')

    cmd.alter(f'clust*', 'b=50')
    cmd.alter(f'clust*', 'q=1')
    cmd.set('gaussian_resolution', 1)


    for atype in atype_l:
        color_clust_by_rank(atype, 'clust')
        
        surf_l = []
        for obj in cmd.get_object_list(f'clust.*.{atype}'):
            cmd.set_name(obj, f'ats_{obj}')
            cmd.map_new(name=f'mapats{obj}', type='gaussian', grid=0.5, selection=f'ats_{obj}')
            cmd.isosurface(f'{obj}', f'mapats{obj}')
            surf_l.append(f'{obj}')
            

        density_clust_by_size(atype, 'clust', obj_l=surf_l)
    
    cmd.group('lig_ats', 'ats_*')
    cmd.group('maps', 'mapats*')

    # Create scenes for visualization
    cmd.enable()
    cmd.enable('receptor')
    cmd.orient()
    cmd.disable('mapats*')
    cmd.hide('everything', 'mapats*')

    cmd.disable('ats_*')
    cmd.hide('everything', 'ats*')

    cmd.color('grey70', 'receptor and elem C')
    cmd.show('surface', f'clust.*')

    # Scenes for each atom type
    for atype in atype_l:
        cmd.disable('clust.*')
        cmd.enable(f'clust.*.{atype}')
        cmd.scene(f'{atype}', 'store')
    
    cmd.disable('clust.*')

    # Top 3, Top5, and all pharmacophore features
    for i in range(0, 3):
        cmd.enable(f'clust.{i}.*')

    cmd.scene('Pharm_Top3', 'store')
    
    for i in range(3, 5):
        cmd.enable(f'clust.{i}.*')
    
    cmd.scene('Pharm_Top5', 'store')

    cmd.enable(f'clust.*')
    cmd.scene('Pharm_All', 'store')
    
    cmd.scene('Pharm_Top3')
    cmd.save(args.out_pse)

if __name__=='__main__':
    main()
