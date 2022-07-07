from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, Structure, Lattice, ReconstructionGenerator, center_slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar,Incar,Kpoints
from pymatgen.core.structure import Structure, Molecule
import glob,os
import numpy as np
from mp_api import MPRester
USER_API_KEY="JBkm9BAbwqy33jXXfIsRhz8C0R2v9gkI"
def my_write_gjf(fname, atoms, adsorbate, ads_struc):
    atom_z = []
    for site_i in range(len(ads_struc)-len(adsorbate)):
        atom_z.append(ads_struc[site_i].c)
    #for i in range(len(adsorbate)):
    #    atom_z.remove(atom_z[-1-i])
    atom_z.sort()
    #print(len(atoms),len(atom_z))
    z_value_fix = 0.5*atom_z[0] + 0.5*atom_z[-1]
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom_i in range(len(atoms)):
        atom = atoms[atom_i]
        if ads_struc[atom_i].c > z_value_fix:
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
def my_write_gjf2(fname, atoms, adsorbate, ads_struc):
    atom_z = []
    for site_i in range(len(ads_struc)-len(adsorbate)):
        atom_z.append(ads_struc[site_i].c)
    #for i in range(len(adsorbate)):
    #    atom_z.remove(atom_z[-1-i])
    atom_z.sort()
    #print(len(atoms),len(atom_z))
    z_value_fix = 0.5*atom_z[0] + 0.5*atom_z[-1]
    f = open(fname,'w')
    f.write('#\n\nOK\n\n0 1\n')
    for atom_i in range(len(atoms)-len(adsorbate)):
        atom = atoms[atom_i]
        if ads_struc[atom_i].c > z_value_fix:
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

def is_slab_symmetrical(slab):
    atom_list = slab.species
    coords_old = slab.frac_coords
    mi_string = "".join([str(i) for i in slab.miller_index])
    slab_layer = []
    print(mi_string)
    print(atom_list)
    print(coords_old)
    for i in range(len(coords_old)):
        print("%.8f" % (coords_old[i][2]))
        if "%.8f" % (coords_old[i][2]) not in slab_layer:
            print(slab_layer,"%.8f" % (coords_old[i][2]))
            slab_layer.append("%.8f" % (coords_old[i][2]))
    slab_layer.sort()
    number_of_slab_layers = len(slab_layer)
    print(number_of_slab_layers)
    for j in range(int(number_of_slab_layers / 2)):
        slab_layer_i = []
        slab_layer_j = []
        for k in range(len(coords_old)):
            print(coords_old[k])
            if "%.8f" % (coords_old[k][2]) == slab_layer[j]:
                slab_layer_i.append(atom_list[k])
            elif "%.8f" % (coords_old[k][2]) == slab_layer[-1 - j]:
                slab_layer_j.append(atom_list[k])
            else:
                continue
        slab_layer_i.sort()
        slab_layer_j.sort()
        print(slab_layer_j)
        print(slab_layer_i)
        if slab_layer_i == slab_layer_j:
            continue
        else:
            return False
    return True
def slab_flipping(slab):
    Lattice_new = Lattice.from_parameters(a=slab.lattice.a, b=slab.lattice.b, c=slab.lattice.c,
                                          alpha=180 - slab.lattice.alpha, beta=180 - slab.lattice.beta,
                                          gamma=slab.lattice.gamma)
    flip_array = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    atom_list = slab.species
    coords_old = slab.frac_coords
    coords_new = np.dot(coords_old, flip_array)
    slab_new = Structure(Lattice_new, atom_list, coords_new)
    return slab_new
def gen_slab_from_mp_id(mp_id,adsorbate):
    with MPRester(USER_API_KEY) as m:
        doc = m.summary.get_data_by_id(mp_id)
        structure = doc.structure
        print(structure)
        print(doc.formula_pretty)
    structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    slabs = generate_all_slabs(structure, max_index=1, min_slab_size=8.0, min_vacuum_size=15.0, center_slab=True,
                               symmetrize=False)
    slab_count = 1
    alloy_name = mp_id + '_' + doc.formula_pretty
    if not os.path.exists(alloy_name):
        os.mkdir(alloy_name)
    for slab in slabs:
        mi_string = "".join([str(i) for i in slab.miller_index])
        slab = center_slab(slab)
        slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
        sg = SpacegroupAnalyzer(slab, symprec=0.01, angle_tolerance=5.0)
        cs = sg.get_crystal_system()
        asf = AdsorbateSiteFinder(slab)
        ads_sites = asf.find_adsorption_sites()
        #adsorbate = O
        ads_structs = asf.generate_adsorption_structures(adsorbate, min_lw=6.0)
        num = len(ads_structs)
        print('%i adducts in total.' % num)
        # print(ads_structs[1])
        # for a in ads_structs[1]:
        #    print(a)
        # print(ads_structs[1][0])
        for i in range(num):
            atoms = AseAtomsAdaptor.get_atoms(ads_structs[i])
            # print(atoms.positions)
            fname = slab_name + '_%i.gjf' % (i)
            my_write_gjf(fname, atoms, adsorbate, ads_structs[i])
            os.system('move ' + fname + ' ' + alloy_name)
        my_write_gjf2(slab_name + '.gjf', atoms, adsorbate, ads_structs[i])
        os.system('move ' + slab_name + '.gjf ' + alloy_name)
        slab_count += 1
    return

H = Molecule("H", [[0, 0, 0]])
O = Molecule("O", [[0, 0, 0]])
OH = Molecule("OH", [[0, 0, 0], [-0.793, 0.384, 0.422]])
OOH = Molecule("OOH", [[0, 0, 0], [-1.067, -0.403, 0.796], [-0.696, -0.272, 1.706]])

def gen_slabs_from_cif(cif_file):
    alloy_name = cif_file.split('.')[0]
    struct = Structure.from_file(cif_file)
    struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
    slabs = generate_all_slabs(struct, max_index=1, min_slab_size=8.0, min_vacuum_size=15.0, center_slab=False, symmetrize=False)
    if not os.path.exists(alloy_name):
        os.mkdir(alloy_name)
    return slabs, alloy_name
def gen_adsorbate_on_slab(slab,adsorbate,alloy_name,slab_count):
    mi_string = "".join([str(i) for i in slab.miller_index])
    slab = center_slab(slab)
    slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
    sg = SpacegroupAnalyzer(slab, symprec=0.01, angle_tolerance=5.0)
    cs = sg.get_crystal_system()
    asf = AdsorbateSiteFinder(slab)
    ads_sites = asf.find_adsorption_sites()
    ads_structs = asf.generate_adsorption_structures(adsorbate, min_lw=6.0)
    num = len(ads_structs)
    print('%i adducts in total.'%num)
    for i in range(num):
        atoms = AseAtomsAdaptor.get_atoms(ads_structs[i])
        fname = slab_name + '_%i.gjf'%(i)
        my_write_gjf(fname, atoms, adsorbate, ads_structs[i])
        os.system('move ' + fname + ' ' + alloy_name)
        my_write_gjf2(slab_name+'.gjf', atoms, adsorbate, ads_structs[i])
        os.system('move ' + slab_name + '.gjf ' + alloy_name)
        slab_count += 1
    return
def slab_io_from_cifs(cif_file):
    slabs, alloy_name = gen_slabs_from_cif(cif_file)
    mi_string_old = 'old'
    slab_count = 1
    for slab in slabs:
        mi_string = "".join([str(i) for i in slab.miller_index])
        print(mi_string)
        if mi_string == mi_string_old:
            slab_count += 1
            slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
        else:
            mi_string_old = mi_string
            slab_count = 1
            slab_name = alloy_name + '_' + mi_string + '_' + '%i' % (slab_count)
        if is_slab_symmetrical(slab):
            center_slab(slab)
            fname1 = slab_name+"u.cif"
            CifWriter(slab).write_file(fname1)
            atoms = AseAtomsAdaptor.get_atoms(slab)
            adsorbate = []
            fname2 = slab_name + "u.gjf"
            my_write_gjf(fname2, atoms, adsorbate, slab)
            os.system('move ' + fname1 + ' ' + alloy_name)
            os.system('move ' + fname2 + ' ' + alloy_name)
        else:
            slab_new = slab_flipping(slab)
            center_slab(slab_new)
            center_slab(slab)
            fname1a = slab_name + "u.cif"
            fname1b = slab_name + "d.cif"
            CifWriter(slab).write_file(fname1a)
            CifWriter(slab_new).write_file(fname1b)
            atoms_u = AseAtomsAdaptor.get_atoms(slab)
            atoms_d = AseAtomsAdaptor.get_atoms(slab_new)
            adsorbate = []
            fname2a = slab_name + "u.gjf"
            fname2b = slab_name + "d.gjf"
            my_write_gjf(fname2a, atoms_u, adsorbate, slab)
            my_write_gjf(fname2b, atoms_d, adsorbate, slab)
            os.system('move ' + fname1a + ' ' + alloy_name)
            os.system('move ' + fname1b + ' ' + alloy_name)
            os.system('move ' + fname2a + ' ' + alloy_name)
            os.system('move ' + fname2b + ' ' + alloy_name)
    return
if __name__ == '__main__':
    print("slab")
    params = {"SYSTEM": "VASP-OPT","ENCUT": 500,"EDIFF":1e-5, \
              "SIGMA":0.2,"ICHARG":2,"ENCUT":450,"ISTART":0,\
              "ISMEAR":1,"ALGO":"Fast","PREC": "Normal", \
              "EDIFF":1E-3,"ISPIN":1,"LREAL":"Auto","NSW":500,\
              "IBRION":2,"EDIFFG":0.2,"ISIF":2,"IVDW":12,"NPAR":12}
    print(Incar(params))
    filename = 'INCAR'
    Incar(params).write_file(filename)
    Mag = True
    if Mag:
        params["ISPIN"]=2
        Incar(params).write_file(filename+"_1")
    Relax_lattice = True
    if Relax_lattice:
        params["ISIF"]=3
        Incar(params).write_file(filename+"_2")