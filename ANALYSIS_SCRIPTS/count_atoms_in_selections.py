import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf


pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
universe = mda.Universe(pdb)
selection_all = universe.select_atoms("all")
atom_counter = 0
for selection_name, selection_keyword in sd.SELECTION_DICT.items():
    pdb_selection = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
    pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
    universe = mda.Universe(pdb)
    selection_keyword = uf.get_selection_keyword_indices(pdb_selection, selection_name, selection_keyword)
    if selection_keyword is not None:
        selection = universe.select_atoms(selection_keyword)
        atom_counter += selection.n_atoms
        print(f"Added {selection.n_atoms} from {selection_name}")
print(f"Total atoms in selections: {atom_counter}")
print(f"Total atoms: {selection_all.n_atoms}")
