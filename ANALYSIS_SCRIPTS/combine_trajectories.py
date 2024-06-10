import MDAnalysis as mda
import utility_functions as uf
import shared_data as sd
import numpy as np
from MDAnalysis.coordinates.memory import MemoryReader


# if os.path.exists(f"{sd.BASE_DIR}/cross_validation_indices_train.npy"):
#     combination_indices_arrays = {
#         "train": np.load(f"{sd.BASE_DIR}/cross_validation_indices_train.npy"),
#         "validation": np.load(f"{sd.BASE_DIR}/cross_validation_indices_validation.npy"),
#         "test": np.load(f"{sd.BASE_DIR}/cross_validation_indices_test.npy"),
#     }
# else:
#     combination_indices_arrays = uf.get_combination_indices()
#     for partition_name, combination_indices_array in combination_indices_arrays.items():
#         np.save(f"{sd.BASE_DIR}/cross_validation_indices_{partition_name}.npy", combination_indices_array)

combination_indices_arrays = uf.get_combination_indices()
for selection_name in sd.SELECTION_DICT:
    pdb = f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb"
    for partition_name, combination_indices_array in combination_indices_arrays.items():
        for cv_index in range(sd.N_CROSS_VALIDATION_COMBINATIONS):
            for pair_name, pair_systems in sd.PAIRS.items():
                xtc_list = []
                for system_name, system_path in pair_systems.items():
                    for trj_index in combination_indices_array[cv_index]:
                        xtc = f"{system_path}/T{trj_index}/{selection_name}.xtc"
                        xtc_list.append(xtc)
                universe = mda.Universe(pdb, xtc_list)
                selection = universe.select_atoms("all")
                if sd.NORMALIZE_BY_AVERAGE:
                    X = np.zeros((len(universe.trajectory), selection.n_atoms, 3))
                    for frame_index, _ in enumerate(universe.trajectory):
                        x = selection.positions
                        average_norm = np.mean(np.linalg.norm(x, axis=1))
                        X[frame_index, :, :] = x / average_norm
                    universe = mda.Universe(pdb, X, format=MemoryReader, order="fac")
                    selection = universe.select_atoms("all")
                selection.write(f"{sd.BASE_DIR}/{pair_name}/{selection_name}.pdb")
                selection.write(f"{sd.BASE_DIR}/{pair_name}/{selection_name}_{partition_name}_{cv_index}.xtc", frames="all")
                print(f"wrote {sd.BASE_DIR}/{pair_name}/{selection_name}_{partition_name}_{cv_index}.xtc")