#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 11:07:15 2019

@author: Hemanta
Here I Compare two approaches of getting neighbors, one by using get_bonded structure
and theother by using CrystalNN()
To get the conncection properties by using bonded structure just need to uncomment the commented part.
This works for compact oxides set very well. I did the test for 200 compounds.
In this script I am using same approach to get the neighbors and the site_geometry
This is the final script for the generation of big table.
The modification is using new criteria to get the motif type
"""
import time
start = time.time()
import os
import pandas as pd
import re, ast
import numpy as np
from pymatgen import Structure
from pymatgen.analysis.local_env import CrystalNN
from robocrys.condense.site import SiteAnalyzer
from PrincipalComponent import PrincipalComponent
from collections import Counter
from robocrys.condense.fingerprint import get_site_fingerprints
Principal_Component_Matrix = PrincipalComponent()
from pymatgen.io.vasp import Poscar

#read the list of the data from your folder
os.chdir("/Users/Hemanta/Desktop/20_test_comp")
#f = open("file_data.dat", "r")
f = open("file_test_motif", "r")
not_worked = open("not_worked_global_table.txt", "w")
for cline, lines in enumerate(f):
    try:
        header = lines.split()[0]
        print(cline, header)
        
                
        if ".cif" in header:
            structure_from_cif = Structure.from_file(header)
        else:
            poscar = Poscar.from_file(str(header))
            structure_from_cif = poscar.structure

        latvec_a, latvec_b, latvec_c = structure_from_cif.lattice.matrix
    
    
        anion = structure_from_cif.species[-1]
        elements = [structure_from_cif.sites[i].species_string for i in range(structure_from_cif.num_sites)]
        element_remove = [unwanted_motif for unwanted_motif in elements if unwanted_motif.startswith(str(anion))]
        final_elem_list = elements[0: -len(element_remove)]
        
        y = CrystalNN()
        

        nn_info_data = [y.get_nn_info(structure_from_cif, i) for i in range(len(final_elem_list))]
        
        nn_info_elem_list = [[coord_site['site'].species_string for coord_site in nn_info_data[i]] for i in range(len(final_elem_list))]
        
        nn_info_neighbor_list = [[str(structure_from_cif.species[i]) +str(i), str([coord_site['site'] for coord_site in nn_info_data[i]])] for i in range(len(final_elem_list))]
        
          
        nn_info_header = [[str(structure_from_cif.species[i]), str([coord_site['site'] for coord_site in nn_info_data[i]])] for i in range(len(final_elem_list))]
        
    
        final_key_list1 = [str(i[0]) + ''.join('%s%d' % t for t in Counter(j).items()) for i, j in zip(nn_info_header, nn_info_elem_list)] 
    
        final_key_list=  [str(i[0]) + ''.join('%s%d' % t for t in Counter(j).items()) for i, j in zip(nn_info_neighbor_list, nn_info_elem_list)] 
        
        nn_info_surrounding_data = [[ast.literal_eval(i) for i in re.findall(r'\([ ,.\d-]+\)', j[1])] for j in nn_info_neighbor_list]
        
        data = dict(zip(final_key_list, nn_info_surrounding_data))
        
            
#        final_key_list=  [str(i[0]) + ''.join('%s%d' % t for t in Counter(j).items()) for i, j in zip(nn_info_neighbor_list, nn_info_elem_list)] 
        
        
        matrix = np.zeros(shape = (len(final_key_list) , len(final_key_list)))
         
            
        for i in range(len(final_key_list)):
            array1 = np.array(data[final_key_list[i]])
            for j in range(len(final_key_list)):
                for nx in range(-1, 2):
                    for ny in range(-1, 2):
                        for nz in range(-1, 2):
                            n_common = 0
                            array2_shifted = np.array([np.array(i) + np.round(nx * np.array(latvec_a),4) + np.round(ny * np.array(latvec_b),4) + 
                                                   np.round(nz * np.array(latvec_c),4) for i in data[final_key_list[j]]])
        
                            for neighbor_1 in array1:
                                for neighbor_2 in array2_shifted:
                                    if np.allclose(neighbor_1, neighbor_2, rtol=1e-02, atol=1e-04, equal_nan=False):
                                        n_common += 1
        
                            if nx == ny == nz == 0 and i == j:
                                pass
                            elif n_common != 0:
                                matrix[i][j] = n_common
                                
                     
        pd.set_option("display.max_rows", len(final_key_list))
        pd.set_option("display.max_columns", len(final_key_list))
        pd.set_option("display.width", 200)
        df = pd.DataFrame(matrix, index = final_key_list, columns = final_key_list)
        print(df)

     

        def get_motif_type(structure_from_cif):
            new_list = []
            
            for j in range(len(final_elem_list)):
                site_fingerprint = list(get_site_fingerprints(structure_from_cif)[j].items())
                parameter1 = max(site_fingerprint,
                               key=lambda x: x[1] if "wt" not in x[0] else 0)
            
                parameter2 = max(site_fingerprint,
                               key=lambda x: x[1] if "wt" in x[0] else 0)
            
                y1 = parameter2[0].split()[-1]
                
                if parameter2[1] > parameter1[1]:
                
                    parameter_3 = max(site_fingerprint,
                                   key=lambda x: x[1] if str(y1) in x[0].split()[-1] and "wt" not in x[0]  else 0)
                
                    geometry = " ".join(parameter_3[0].split()[:-1])
                    new_list.append(geometry)
                else:
                    geometry = " ".join(parameter1[0].split()[:-1])
                    new_list.append(geometry)
            return new_list         
        
# =============================================================================
#         
#         def get_motif_type(structure_from_cif):
#             new_list = []
#             y1 = y.get_bonded_structure(structure_from_cif)
#             for i in range(len(final_key_list1)):
#                 motif_type = SiteAnalyzer(y1).get_site_geometry(i)
#                 new_list.append(motif_type['type'])
#             return new_list
# =============================================================================
    
        motif_types = ['_'.join(map(str, i)) for i in zip(final_key_list1, get_motif_type(structure_from_cif))]
        
        
        print(motif_types)
        visited_list = []
        for ix, iy in np.ndindex(matrix.shape):
            
             if (not int(matrix[ix][iy]) is 0) and str(final_key_list1[ix])+ "-X-" + str(final_key_list1[iy]) not in visited_list:
                 visited_list.append(str(final_key_list1[ix])+ "-X-" + str(final_key_list1[iy]))
    
                 Principal_Component_Matrix.AddANode(final_key_list1[ix], matrix[ix][iy], 
                    motif_types[iy], motif_types[ix])
    
        type_df = pd.DataFrame(Principal_Component_Matrix.Graph).T
        writer = pd.ExcelWriter('scompound.xlsx',
                                   engine='xlsxwriter',
                                    options={'strings_to_urls': False})
        type_df.to_excel(writer,  sheet_name = 'Sheet1' , na_rep = 'Nan')
        writer.save()
        print("test")
    except:
        not_worked.write(str(header)+"\n")

end = time.time() 
print("Time taken", end - start)

