import MDAnalysis as mda
import shared_data as sd


n_atoms_list = []
for selection_name in sd.SELECTION_DICT:
    pdb = f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb"
    universe = mda.Universe(pdb)
    selection = universe.select_atoms("all")
    n_atoms_list.append((selection_name, selection.n_atoms))
    print(f"{selection_name}: N atoms = {selection.n_atoms}")
print("Largest selection:")
print(sorted(n_atoms_list, key=lambda x: x[1], reverse=True)[0])
