
import numpy as np
from biopandas.pdb import PandasPdb
import sys
import pandas as pd
import Bio.PDB
import json
import os
pdb_name = sys.argv[1]
pdb_path = sys.argv[2]
surface_path = sys.argv[3]
center_path = sys.argv[4]
full_path = sys.argv[5]




aa_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
           'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
           'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


def getKeysByValue(dict_elements, value_find):
    keys_list = list()
    items_list = dict_elements.items()
    for item  in items_list:
        if item[1] == value_find:
            keys_list.append(item[0])
    return  keys_list[0]

def getCenterPoint(center, pdb_atoms):
    as_points = []
    y_sum = 0
    x_sum = 0
    z_sum = 0

    for j in center['res']:
        as_point = int(j)
        as_coords = pdb_atoms[pdb_atoms['residue_number'] == as_point]
        try:
            as_coords = as_coords[as_coords["atom_name"] == "CA"]
            print('try')
            print(as_coords)
            nbr = int(as_coords["atom_number"]) - 1
        except:
            try:
                as_coords = as_coords[as_coords["atom_name"] == "CA"]
                as_coords = as_coords[as_coords["alt_loc"] == "A"]
                print(as_coords)
                nbr = int(as_coords["atom_number"]) - 1
            except:
               print('Bad input structure')

        as_points.append(as_coords)
        x_sum = x_sum + as_coords["x_coord"][nbr]
        y_sum = y_sum + as_coords["y_coord"][nbr]
        z_sum = z_sum + as_coords["z_coord"][nbr]

    center_coords = np.array((x_sum / len(center['res']), y_sum / len(center['res']), z_sum / len(center['res'])), dtype=float)
    return center_coords

def getPathToSurface(center, surface_nbr, sur_try):
    old_center = center
    surface = sur_try['res'][surface_nbr]
    path_name = sur_try['AA'][surface_nbr]
    first_coords = pdb_atoms[pdb_atoms['residue_number'] == surface]

    try:
        first_coords = first_coords[first_coords["atom_name"] == "CA"]
        nbr = int(first_coords["atom_number"])- 1
    except:
        first_coords = first_coords[first_coords["atom_name"] == "CA"]
        first_coords = first_coords[first_coords["alt_loc"] == "A"]
        nbr = int(first_coords["atom_number"]) - 1
        print(first_coords)
    surf_coords = np.array((first_coords["x_coord"][nbr], first_coords["y_coord"][nbr], first_coords["z_coord"][nbr]),
                           dtype=float)
    distance = np.linalg.norm(center - surf_coords)

    x_dif = (center[0] - first_coords["x_coord"])[nbr] / (distance * 2)
    y_dif = (center[1] - first_coords["y_coord"])[nbr] / (distance * 2)
    z_dif = (center[2] - first_coords["z_coord"])[nbr] / (distance * 2)

    curent_path = {}
    aa_first = sur_try['AA'][surface_nbr]
    curent_path[int(surface)] = getKeysByValue(aa_dict, aa_first)
    for i in range(0, int(distance * 2)):
        point = np.array((center[0] - x_dif, center[1] - y_dif, center[2] - z_dif), dtype=float)
        close_atoms = ns.search(point, 2)
        try:
            for j in range(0, len(close_atoms) + 1):
                curent_path[close_atoms[j].parent.id[1]] = str(close_atoms[j].parent.resname)
        except:
            pass
        center = point
    center = old_center
    f_name = str(path_name) + str(int(surface)) + ".txt"
    pymol_text ='sele resid '
    new_dir = os.path.join(pdb_path, 'detected_paths_'+str(pdb_name))
    try:
        os.makedirs(new_dir)
    except OSError:
        pass
    path = os.path.join(new_dir, f_name)
    with open(path, 'w') as file:
        for i in curent_path.keys():
            file.write(str(i) + ' ' + str(curent_path[i]) + '\n')
            pymol_text = pymol_text + str(i) + '+'
    return pymol_text[0:-1]




# pdb_name = '5uro'
# pdb_path = 'C:/Users/hp/Desktop/balcony_surface/'
# surface_path = 'C:/Users/hp/Desktop/balcony_surface/5uro_surface.txt'
# center_path = 'C:/Users/hp/Desktop/select/AS/5uro_active_site.txt'

all_path = pdb_path+ pdb_name + '.pdb'
pdb_try = PandasPdb().read_pdb(full_path)
sur_try = pd.read_csv(str(surface_path), delimiter=r"\s+")
sur_try = sur_try.dropna()
center = pd.read_fwf(str(center_path))
pdb_atoms = pdb_try.df["ATOM"]
parser = Bio.PDB.PDBParser(QUIET=True)
structures = parser.get_structure(str(pdb_name), full_path)
structure = structures[0]
atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
ns = Bio.PDB.NeighborSearch(atoms)

new_center = getCenterPoint(center, pdb_atoms)

df_all = pd.DataFrame(columns=['aa_surface', 'pymol sele'])

print(sur_try.columns)
for n in range(0,len(sur_try['res'])):
    line_pymol = getPathToSurface(new_center, n, sur_try)
    df_all.loc[n] = [sur_try['res'][n]] + [line_pymol]

df_all.to_csv('pymol_selection' + str(pdb_name) + '.csv')