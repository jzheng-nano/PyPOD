#from pymatgen.ext.matproj import MPRester
import os.path,os

from pymatgen.io.cif import CifWriter
from mp_api import MPRester
USER_API_KEY="JBkm9BAbwqy33jXXfIsRhz8C0R2v9gkI"
#elements_list = ["Al", "Si", "P", "S",\
            #"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",\
            #"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",\
            #"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"]
elements_list2 = ["Al", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "Se",\
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"]
elements_list3 = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", \
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au"]
elements_list1 = ["Ru", "Os", "Co", "Ni"]# "Cu", "Rh", "Pd", "Pt", "Ag", "Au", "Ir"]
elements_list = elements_list3
with MPRester(USER_API_KEY) as m:
    results_a = []
    results_b = []
    for element_i in range(len(elements_list)):
        for element_j in range(element_i+1,len(elements_list)):
            #for element3 in elements_list:
            #if element1 != element2:  # and element1 != element3 and element2 != element3:
            results = m.summary.search(chemsys=elements_list[element_i]+"-"+elements_list[element_j],
                                       #nelements_max=2,
                                       #is_magnetic=False,
                                       fields=["formula_pretty",
                                                "material_id",
                                                "is_stable",
                                                "theoretical",
                                                "nelements",
                                                "is_magnetic",
                                                "nsites",
                                                "structure"])
            for i in range(len(results)):
                result_0 = results[i]
                #print(result_0.material_id,result_0.nelements)
                #print(result_0.structure)
                if result_0.nelements==2:
                    if result_0.is_stable==True or result_0.theoretical==False:  # and results[i].get('spacegroup.symbol')=='Fm-3m':
                        results_a.append(results[i])
    for element_i in range(len(elements_list)):
        for element_j in range(element_i+1,len(elements_list)):
            for element_k in range(element_j + 1, len(elements_list)):
            #for element3 in elements_list:
            #if element1 != element2:  # and element1 != element3 and element2 != element3:
                results = m.summary.search(chemsys=elements_list[element_i] + "-" + elements_list[element_j] + "-" + elements_list[element_k],
                                           fields=["formula_pretty",
                                                    "material_id",
                                                    "is_stable",
                                                    "theoretical",
                                                    "nelements",
                                                    "is_magnetic",
                                                    "nsites",
                                                    "structure"])
                for i in range(len(results)):
                    result_0 = results[i]
                    #print(result_0)
                    #
                    if result_0.nelements==3:
                        if result_0.is_stable==True or result_0.theoretical==False:  # and results[i].get('spacegroup.symbol')=='Fm-3m':
                            results_b.append(results[i])
    print(len(results_a))
    #for j in range(10):
    #    print(results[j].material_id, results[j].formula_pretty,results[j].is_magnetic)
    fn = open("BinaryAlloy.txt", "w")
    fn1 = open("TernaryAlloy.txt", "w")
    if not os.path.exists("BinaryAlloy"):
        os.mkdir("BinaryAlloy")
    if not os.path.exists("TernaryAlloy"):
        os.mkdir("TernaryAlloy")
    for i in range(len(results_a)):
        #print(results_a[i].material_id, results_a[i].formula_pretty,results_a[i].is_magnetic,results_a[i].nsites)
        fn.write(str(results_a[i].material_id) + '   ')
        fn.write(results_a[i].formula_pretty + '   ')
        fn.write(str(results_a[i].is_magnetic) + '    ')
        fn.write(str(results_a[i].theoretical) + '    ')
        fn.write(str(results_a[i].is_stable) + '    ')
        fn.write(str(results_a[i].nsites) + '\n')
        fname = results_a[i].material_id + "_" + results_a[i].formula_pretty + ".cif"
        CifWriter(results_a[i].structure).write_file(fname)
        print(fname)
        os.system('move ' + fname + ' BinaryAlloy')
        print('move ' + fname + ' BinaryAlloy')
    for i in range(len(results_b)):
        #print(results_b[i].material_id, results_b[i].formula_pretty,results_b[i].is_magnetic,results_b[i].nsites)
        fn1.write(str(results_b[i].material_id) + '   ')
        fn1.write(results_b[i].formula_pretty + '   ')
        fn1.write(str(results_b[i].is_magnetic) + '    ')
        fn1.write(str(results_b[i].theoretical) + '    ')
        fn1.write(str(results_b[i].is_stable) + '    ')
        fn1.write(str(results_b[i].nsites) + '\n')
        fname = results_b[i].material_id + "_" + results_b[i].formula_pretty + ".cif"
        CifWriter(results_b[i].structure).write_file(fname)
        os.system('move ' + fname + ' TernaryAlloy')
    fn.close()
    fn1.close()



