import shared_data as sd
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import utility_functions as uf


@uf.get_timing
def align_trajectories():
    for trj_index in range(1, sd.N_TRJ+1):
        for system_path in sd.SYSTEMS.values():
            mobile = mda.Universe(f"{system_path}/common_atoms.pdb",
                                f"{system_path}/T{trj_index}/common_atoms.xtc")
            # only one common reference
            reference = mda.Universe(f"{sd.SYSTEMS['LNC']}/common_atoms.pdb")
            # align to common CA and P
            alignment_selection_string = f"({sd.PREALIGNMENT_SELECTION})"
            print(f"aligning {system_path}/T{trj_index}/common_atoms.xtc")
            AlignTraj(mobile, reference, alignment_selection_string, 
                    filename=f"{system_path}/T{trj_index}/aligned_common_atoms.xtc").run()
            print(f"aligned, writing {system_path}/T{trj_index}/aligned_common_atoms.xtc")


align_trajectories()
