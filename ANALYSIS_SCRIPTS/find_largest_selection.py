import MDAnalysis as mda
import shared_data as sd
from utility_functions import get_timing


@get_timing
def find_largest_selection():
    n_atoms_list = []
    for selection_name in sd.SELECTION_DICT:
        pdb = f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb"
        universe = mda.Universe(pdb)
        selection = universe.select_atoms("all")
        n_atoms_list.append((selection_name, selection.n_atoms))
        print(f"{selection_name}: N atoms = {selection.n_atoms}")
    print("Largest selection:")
    largest_selection = sorted(n_atoms_list, key=lambda x: x[1], reverse=True)[0]
    print(f"{largest_selection[0]}: {largest_selection[1]} atoms")

find_largest_selection()
