from pymatgen.core.structure import Structure, Molecule
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator, center_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure, Molecule
import glob,os

def my_write_gjf(fname, atoms, adsorbate):
    atom_z = []
    for atom in atoms:
        atom_z.append(atom.position[2])
    for i in range(len(adsorbate)):
        atom_z.remove(atom_z[-1-i])
    atom_z.sort()
    print(len(atoms),len(atom_z))
    z_value_fix = 0.5*atom_z[0] + 0.5*atom_z[-1]
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom in atoms:
        if atom.position[2] > z_value_fix:
            f.write('%s %s %f %f %f\n'%(atom.symbol, str(0), atom.position[0], atom.position[1], atom.position[2]))
        else:
            f.write('%s %s %f %f %f\n'%(atom.symbol, str(-1), atom.position[0], atom.position[1], atom.position[2]))
    cell = atoms.cell
    for i in cell:
        f.write('Tv ')
        for j in i:
            f.write('%f '%j)
        f.write('\n')
    f.write('\n')
    f.close()

H = Molecule("H", [[0, 0, 0]])
O = Molecule("O", [[0, 0, 0]])
OH = Molecule("OH", [[0, 0, 0], [-0.793, 0.384, 0.422]])
OOH = Molecule("OOH", [[0, 0, 0], [-1.067, -0.403, 0.796], [-0.696, -0.272, 1.706]])
cif_lists = glob.glob('*.cif')

cif_lists = glob.glob('*.cif')
print(cif_lists)
for alloy_species in cif_lists:
    alloy_name = alloy_species.split('.')[0]
    struct = Structure.from_file(alloy_species)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    slabs = generate_all_slabs(struct, max_index=2, min_slab_size=7.0, min_vacuum_size=15.0, center_slab=True)
    slab_count = 1
    os.mkdir(alloy_name)
    for slab in slabs:
        #slab = slab.get_orthogonal_c_slab()
        mi_string = "".join([str(i) for i in slab.miller_index])
        atoms = AseAtomsAdaptor.get_atoms(slab)
        if len(atoms) < 50:
            fname = alloy_name + '_' + mi_string + '_' + '%i.gjf'%(slab_count)
            my_write_gjf(fname, atoms)
            slab_count += 1
            os.system('move ' + fname + ' ' + alloy_name)
        else:
            continue

for slab_cif in cif_lists:
    slab_name = slab_cif.split('.')[0]
    slab1 = Structure.from_file(slab_cif)
    slab = center_slab(slab1)
    sg = SpacegroupAnalyzer(slab, symprec=0.01, angle_tolerance=5.0)
    cs = sg.get_crystal_system()
    print(cs)
    asf = AdsorbateSiteFinder(slab)
    ads_sites = asf.find_adsorption_sites()
    adsorbate = O
    print(adsorbate)
    print(len(adsorbate))
    ads_structs = asf.generate_adsorption_structures(adsorbate, min_lw=6.0)
    num = len(ads_structs)
    print('%i adducts in total.'%num)
    for i in range(num):
        atoms = AseAtomsAdaptor.get_atoms(ads_structs[i])
        fname = slab_name + '_%i.gjf'%(i)
        my_write_gjf(fname, atoms, adsorbate)
print('Done2')