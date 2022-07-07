from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator, center_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure, Molecule
import glob,os
import numpy as np
from mp_api import MPRester
USER_API_KEY="JBkm9BAbwqy33jXXfIsRhz8C0R2v9gkI"
from surface_cleave_adsorption import gen_slab_from_mp_id

def is_slab_symmetrical(slab):
    atom_list = slab.species
    coords_old = slab.frac_coords
    slab_layer = []
    for i in range(len(coords_old)):
        if coords_old[i][2] not in slab_layer:
            slab_layer.append(coords_old[i][2])
    slab_layer.sort()
    number_of_slab_layers = len(slab_layer)
    for j in range(int(number_of_slab_layers / 2)):
        slab_layer_i = []
        slab_layer_j = []
        for k in range(len(coords_old)):
            if coords_old[k][2] == slab_layer[j]:
                slab_layer_i.append(atom_list[k])
            elif coords_old[k][2] == slab_layer[-1 - j]:
                slab_layer_j.append(atom_list[k])
            else:
                continue
        slab_layer_i.sort()
        slab_layer_j.sort()
        if slab_layer_i == slab_layer_j:
            continue
        else:
            return False
    return True

mp_id = "mp-2642"
with MPRester(USER_API_KEY) as m:
    doc = m.summary.get_data_by_id(mp_id)
    structure = doc.structure
    #print(structure)
    #print(doc.formula_pretty)
structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
slabs = generate_all_slabs(structure, max_index=1, min_slab_size=8.0, min_vacuum_size=15.0, center_slab=True,
                           symmetrize=True)
slab_count = 1
alloy_name = mp_id + '_' + doc.formula_pretty
for slab in slabs:
    mi_string = "".join([str(i) for i in slab.miller_index])
    if is_slab_symmetrical(slab):
        fname1 = alloy_name + '_' + mi_string + '_u' + '%i.cif' % (slab_count)
        CifWriter(slab).write_file(fname1)
    else:
        Lattice_new = Lattice.from_parameters(a=slab.lattice.a, b=slab.lattice.b, c=slab.lattice.c,
                                          alpha=180-slab.lattice.alpha, beta=180-slab.lattice.beta, gamma=slab.lattice.gamma)
        flip_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
        atom_list = slab.species
        coords_old = slab.frac_coords
        coords_new = np.dot(coords_old,flip_array)
        slab_new = Structure(Lattice_new, atom_list, coords_new)
        print(slab_new)
        fname1 = alloy_name + '_' + mi_string + '_u' + '%i.cif'%(slab_count)
        fname2 = alloy_name + '_' + mi_string + '_d' + '%i.cif'%(slab_count)
        CifWriter(slab).write_file(fname1)
        CifWriter(slab_new).write_file(fname2)
print("Done!")