import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf


# pdb_selection = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
# pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
# for selection_name, selection_keyword in sd.SELECTION_DICT.items():
#     selection_keyword_indices = uf.get_selection_keyword_indices(pdb_selection, selection_name, selection_keyword)
#     for system_name, system_path in sd.SYSTEMS.items():
#         for trj_index in range(1, sd.N_TRJ+1):
#             xtc = f"{system_path}/T{trj_index}/aligned_common_atoms.xtc"
#             universe = mda.Universe(pdb, xtc)
#             selection = universe.select_atoms(selection_keyword_indices)
#             selection.write(f"{system_path}/{selection_name}.pdb")
#             selection.write(f"{system_path}/T{trj_index}/{selection_name}.xtc", frames="all")
#             print(f"wrote {system_path}/T{trj_index}/{selection_name}.xtc, n_atoms={selection.n_atoms}")

pdb_selection = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
selection_keyword_indices = uf.get_selection_keyword_indices(pdb_selection, "sphere_layer1", sd.SELECTION_DICT["sphere_layer1"])
system_path = sd.SYSTEMS['LNC']
xtc = f"{system_path}/T1/aligned_common_atoms.xtc"
universe = mda.Universe(pdb, xtc)
selection = universe.select_atoms(selection_keyword_indices)
selection.write(f"{system_path}/sphere_layer1.pdb")
selection.write(f"{system_path}/T1/sphere_layer1.xtc", frames="all")
print(f"wrote {system_path}/T1/sphere_layer1.xtc, n_atoms={selection.n_atoms}")