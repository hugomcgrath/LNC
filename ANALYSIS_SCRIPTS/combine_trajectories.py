import MDAnalysis as mda
import utility_functions as uf
import shared_data as sd
from multiprocessing import Pool
import os


def combine_trajectories(selection_name):
    pdb = f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb"
    if os.path.exists(pdb):
        for partition_name, combination_indices_array in combination_indices_arrays.items():
            for cv_index in range(sd.N_CROSS_VALIDATION_COMBINATIONS):
                for pair_name, pair_systems in sd.PAIRS.items():
                    xtc_list = []
                    for system_path in pair_systems.values():
                        for trj_index in combination_indices_array[cv_index]:
                            xtc = f"{system_path}/T{trj_index}/{selection_name}.xtc"
                            xtc_list.append(xtc)
                    universe = mda.Universe(pdb, xtc_list)
                    selection = universe.select_atoms("all")
                    selection.write(f"{sd.BASE_DIR}/{pair_name}/{selection_name}.pdb")
                    selection.write(f"{sd.BASE_DIR}/{pair_name}/{selection_name}_{partition_name}_{cv_index}.xtc", frames="all")
                    print(f"wrote {sd.BASE_DIR}/{pair_name}/{selection_name}_{partition_name}_{cv_index}.xtc")


@uf.get_timing
def combine_trajectories_parallel():
    pool = Pool()
    pool.map(combine_trajectories, sd.SELECTION_DICT)


combination_indices_arrays = uf.get_combination_indices()
combine_trajectories_parallel()
