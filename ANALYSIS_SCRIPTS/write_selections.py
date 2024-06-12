import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf
import os
from multiprocessing import Pool


def write_selections(trj_index):
    pdb_selection = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
    pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
    for system_path in sd.SYSTEMS.values():
        for selection_name, selection_keyword in sd.SELECTION_DICT.items():
            selection_keyword_indices = uf.get_selection_keyword_indices(pdb_selection, selection_name, selection_keyword)
            if selection_keyword_indices is not None:
                n_atoms = len(selection_keyword_indices.split(" ")[1:])
                if n_atoms < 50:
                    print(f"selection {selection_name} is too small, n_atoms={n_atoms}")
                else:
                    xtc = f"{system_path}/T{trj_index}/aligned_common_atoms.xtc"
                    universe = mda.Universe(pdb, xtc)
                    selection = universe.select_atoms(selection_keyword_indices)
                    selection.write(f"{system_path}/{selection_name}.pdb")
                    selection.write(f"{system_path}/T{trj_index}/{selection_name}.xtc", frames="all")
                    print(f"wrote {system_path}/T{trj_index}/{selection_name}.xtc, n_atoms={selection.n_atoms}")
            else:
                print(f"No atoms in selection for {system_path}/T{trj_index}/{selection_name}.xtc")


@uf.get_timing
def write_selections_parallel():
    pool = Pool()
    pool.map(write_selections, trj_indices)


trj_indices = range(1, sd.N_TRJ+1)
write_selections_parallel()
