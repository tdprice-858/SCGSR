from ase.io import read, write
from ase.visualize import view
from ase import Atoms
from acat.adsorption_sites import SlabAdsorptionSites
import sys, os
sys.path.insert(0, "/Users/tdprice/Documents/py/code/pynta_main/pynta")
from pynta.main import Pynta
from ase.geometry import get_distances
import numpy as np
from acat.utilities import get_mic
from copy import deepcopy
import ase
from xtb.ase.calculator import XTB
#from deepmd.calculator import DP
from ase.calculators import calculator
from sella import Sella, Constraints


def generate_unique_placements(slab,cas):
    nslab = len(slab)
    middle = sum(slab.cell)/2.0

    unique_single_sites = get_unique_sites(cas,about=middle)

    unique_site_pairs = dict() # (site1site,site1morph,),(site2site,site2morph),xydist,zdist
    for unique_site in unique_single_sites:
        uni_site_fingerprint = (unique_site["site"],unique_site["morphology"])
        for site in cas.get_sites():
            site_fingerprint = (site["site"],site["morphology"])
            bd,d = get_distances([unique_site["position"]], [site["position"]], cell=slab.cell, pbc=(True,True,False))
            xydist = np.linalg.norm(bd[0][0][:1])
            zdist = bd[0][0][2]

            fingerprint = (uni_site_fingerprint,site_fingerprint,round(xydist,3),round(zdist,3))

            if fingerprint in unique_site_pairs.keys():
                current_sites = unique_site_pairs[fingerprint]
                current_dist = np.linalg.norm(sum([s["position"][:1] for s in current_sites])/2-middle[:1])
                possible_dist = np.linalg.norm((unique_site["position"][:1]+site["position"][:1])/2-middle[:1])
                if possible_dist < current_dist:
                    unique_site_pairs[fingerprint] = [unique_site,site]
            else:
                unique_site_pairs[fingerprint] = [unique_site,site]

    unique_site_pairs_lists = list(unique_site_pairs.values())
    unique_site_lists = [[unique_site] for unique_site in unique_single_sites]

    single_site_bond_params_lists = []
    for unique_site_list in unique_site_lists:
        pos = deepcopy(unique_site_list[0]["position"])
        single_site_bond_params_lists.append([{"site_pos": pos,"ind": None, "k": 100.0, "deq": 0.0}])

    double_site_bond_params_lists = []
    for unique_site_pair_list in unique_site_pairs_lists:
        bond_params_list = []
        for site in unique_site_pair_list:
            pos = deepcopy(site["position"])
            bond_params_list.append({"site_pos": pos,"ind": None, "k": 100.0, "deq": 0.0})
        double_site_bond_params_lists.append(bond_params_list)

    return unique_site_lists,unique_site_pairs_lists,single_site_bond_params_lists,double_site_bond_params_lists

def get_unique_sites(cas, unique_composition=False,
                     unique_subsurf=False,
                     return_signatures=False,
                     return_site_indices=False,
                     about=None, site_list=None):
    """Function mostly copied from the ACAT software
    Get all symmetry-inequivalent adsorption sites (one
    site for each type).

    Parameters
    ----------
    unique_composition : bool, default False
        Take site composition into consideration when
        checking uniqueness.

    unique_subsurf : bool, default False
        Take subsurface element into consideration when
        checking uniqueness.

    return_signatures : bool, default False
        Whether to return the unique signatures of the
        sites instead.

    return_site_indices: bool, default False
        Whether to return the indices of each unique
        site (in the site list).

    about: numpy.array, default None
        If specified, returns unique sites closest to
        this reference position.

    """

    if site_list is None:
        sl = cas.site_list
    else:
        sl = site_list
    key_list = ['site', 'morphology']
    if unique_composition:
        if not cas.composition_effect:
            raise ValueError('the site list does not include '
                             + 'information of composition')
        key_list.append('composition')
        if unique_subsurf:
            key_list.append('subsurf_element')
    else:
        if unique_subsurf:
            raise ValueError('to include the subsurface element, ' +
                             'unique_composition also need to be set to True')
    if return_signatures:
        sklist = sorted([[s[k] for k in key_list] for s in sl])
        return sorted(list(sklist for sklist, _ in groupby(sklist)))
    else:
        seen_tuple = []
        uni_sites = []
        if about is not None:
            sl = sorted(sl, key=lambda x: get_mic(x['position'],
                                                  about, cas.cell, return_squared_distance=True))
        for i, s in enumerate(sl):
            sig = tuple(s[k] for k in key_list)
            if sig not in seen_tuple:
                seen_tuple.append(sig)
                if return_site_indices:
                    s = i
                uni_sites.append(s)

        return uni_sites

QE_slab = read('/Users/tdprice/slabs/QE_slab.xyz')
vasp_slab = read('/Users/tdprice/slabs/slab.xyz')
new_vasp_slab = read('/Users/tdprice/slabs/vasp_slab_2pbc.xyz')
print(f"QE slab pbc = {QE_slab.pbc}")
print(f"vasp slab pbc = {vasp_slab.pbc}")
print(f"2 pbc vasp slab = {new_vasp_slab.pbc}")
atoms = read('/Users/tdprice/Desktop/Desktop/02_pt-mgo-ethylene/01_bare/opt_PBE_400_221/vasprun.xml', index=':')
#view(atoms)
write('/Users/tdprice/Desktop/Desktop/02_pt-mgo-ethylene/01_bare/opt_PBE_400_221/test.xyz', atoms)
print(os.getcwd())
print(os.path.join(os.getcwd(),'j','xtbharm.traj'))

'''
cas_vasp_3pbc = SlabAdsorptionSites(vasp_slab, "fcc111" ,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal='Ag')
cas_vasp_2pbc = SlabAdsorptionSites(new_vasp_slab, "fcc111" ,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal='Ag')
cas_QE_slab = SlabAdsorptionSites(QE_slab, "fcc111" ,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal='Ag')

unique_site_lists,\
unique_site_pairs_lists,\
single_site_bond_params_lists,\
double_site_bond_params_lists = generate_unique_placements(vasp_slab,cas_vasp_3pbc)

print(f"vasp_3pbc \n unique_site_lists ={unique_site_lists} \n single_site_bond_params_lists = {single_site_bond_params_lists}\
\n double_site_bond_params_lists = {double_site_bond_params_lists}")
print(len(unique_site_lists))
print(len(unique_site_pairs_lists))
print(len(single_site_bond_params_lists))
print(len(double_site_bond_params_lists))

unique_site_lists,\
unique_site_pairs_lists,\
single_site_bond_params_lists,\
double_site_bond_params_lists = generate_unique_placements(new_vasp_slab,cas_vasp_2pbc)


print(f"vasp_2pbc \n unique_site_lists ={unique_site_lists} \n single_site_bond_params_lists = {single_site_bond_params_lists}\
\n double_site_bond_params_lists = {double_site_bond_params_lists}")
print(len(unique_site_lists))
print(len(unique_site_pairs_lists))
print(len(single_site_bond_params_lists))
print(len(double_site_bond_params_lists))

unique_site_lists,\
unique_site_pairs_lists,\
single_site_bond_params_lists,\
double_site_bond_params_lists = generate_unique_placements(QE_slab,cas_QE_slab)

print(f"QE 2 pbc \n unique_site_lists ={unique_site_lists} \n single_site_bond_params_lists = {single_site_bond_params_lists}\
\n double_site_bond_params_lists = {double_site_bond_params_lists}")
print(len(unique_site_lists))
print(len(unique_site_pairs_lists))
print(len(single_site_bond_params_lists))
print(len(double_site_bond_params_lists))


class HarmonicallyForcedXTB():
    def get_energy_forces(self):
        energy = 0.0
        forces = np.zeros(self.atoms.positions.shape)
        if hasattr(self.parameters, "atom_bond_potentials"):
            for atom_bond_potential in self.parameters.atom_bond_potentials:
                E, F = get_energy_forces_atom_bond(self.atoms, **atom_bond_potential)
                energy += E
                forces += F

        if hasattr(self.parameters, "site_bond_potentials"):
            for site_bond_potential in self.parameters.site_bond_potentials:
                E, F = get_energy_forces_site_bond(self.atoms, **site_bond_potential)
                energy += E
                forces += F

        return energy[0][0], forces

    def calculate(self, atoms=None, properties=None, system_changes=calculator.all_changes):
        parent_class.calculate(self, atoms=atoms, properties=properties, system_changes=system_changes)
        energy, forces = self.get_energy_forces()
        self.results["energy"] += energy
        self.results["free_energy"] += energy
        self.results["forces"] += forces
parent_class= XTB
hfxtb = HarmonicallyForcedXTB(parent_class, method="GFN1-xTB")
#vasp_slab.set_constraint(out_constraints)
vasp_slab.calc = hfxtb

print(vasp_slab)
opt = Sella(vasp_slab,trajectory="xtbharm.traj",order=0)
#opt.run(fmax=0.02,steps=1)'''