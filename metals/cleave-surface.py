import os.path

from pymatgen.core.structure import Structure, Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#import nglview as nv
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
#from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure, Molecule
import glob
def my_write_gjf(fname, atoms):
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom in atoms:
        f.write('%s %f %f %f\n'%(atom.symbol, atom.position[0], atom.position[1], atom.position[2]))
    cell = atoms.cell
    for i in cell:
        f.write('Tv ')
        for j in i:
            f.write('%f '%j)
        f.write('\n')
    f.write('\n')

cif_lists = ["mp-1007692_HfPd.cif"]
print(cif_lists)
for alloy_species in cif_lists:
    alloy_name = alloy_species.split('.')[0]
    struct = Structure.from_file(alloy_species)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    slabs = generate_all_slabs(struct, max_index=1, min_slab_size=8.0, min_vacuum_size=15.0, center_slab=True, symmetrize=False)
    slab_count = 1
    #print(os.path.exists(alloy_name))
    if not os.path.exists(alloy_name):
        os.mkdir(alloy_name)
    for slab in slabs:
        #slab = slab.get_orthogonal_c_slab()
        mi_string = "".join([str(i) for i in slab.miller_index])
        atoms = AseAtomsAdaptor.get_atoms(slab)
        if len(atoms) < 50:
            fname1 = alloy_name + '_' + mi_string + '_' + '%i.cif'%(slab_count)
            fname2 = alloy_name + '_' + mi_string + '_' + '%i.gjf' % (slab_count)
            w = CifWriter(slab)
            w.write_file(fname1)
            my_write_gjf(fname2, atoms)
            slab_count += 1
            os.system('move ' + fname1 + ' ' + alloy_name)
            os.system('move ' + fname2 + ' ' + alloy_name)
        else:
            continue
print('Done')