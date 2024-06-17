import shared_data as sd
import MDAnalysis as mda
from multiprocessing import Pool
import utility_functions as uf


# write trajectories with common atoms
def write_common_trajectories(trj_index):
    for system_path in sd.SYSTEMS.values():
        pdb_common = f"{system_path}/common_atoms.pdb"
        universe = mda.Universe(f"{system_path}/system.pdb", 
                                f"{system_path}/T{trj_index}/traj_pruned.xtc")
        selection = universe.select_atoms(
            (
                f"({sd.HEAVY_ATOMS_SELECTION_KEYWORD} and protein)"
                f" or {sd.NUCLEIC_PRUNED_SELECTION_KEYWORD}"
                " and not segid 9"
            )
        )
        selection.write(pdb_common)
        selection.write(f"{system_path}/T{trj_index}/common_atoms.xtc",
                        frames=universe.trajectory[:sd.TRJ_LEN])
        print(f"wrote {system_path}/T{trj_index}/common_atoms.xtc")


@uf.get_timing
def write_common_trajectories_parallel():
    pool = Pool()
    pool.map(write_common_trajectories, trj_indices)


trj_indices = range(1, sd.N_TRJ+1)
write_common_trajectories_parallel()
