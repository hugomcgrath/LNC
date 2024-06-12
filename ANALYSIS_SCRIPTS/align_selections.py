import shared_data as sd
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import subprocess
from multiprocessing import Pool
import utility_functions as uf
import os


def align_selections(trj_index):
    for system_path in sd.SYSTEMS.values():
        for selection_name in sd.SELECTION_DICT:
            pdb_mobile = f"{system_path}/{selection_name}.pdb"
            xtc_mobile = f"{system_path}/T{trj_index}/{selection_name}.xtc"
            if os.path.exists(pdb_mobile) and os.path.exists(xtc_mobile):
                mobile = mda.Universe(pdb_mobile, xtc_mobile)
                # only one common reference
                reference = mda.Universe(f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb")
                # align to atoms of the selection
                alignment_selection_string = f"({sd.ALIGNMENT_SELECTION})"
                AlignTraj(mobile, reference, alignment_selection_string, 
                        filename=f"{system_path}/T{trj_index}/{selection_name}_aligned.xtc").run()
                subprocess.run(f"mv {system_path}/T{trj_index}/{selection_name}_aligned.xtc {system_path}/T{trj_index}/{selection_name}.xtc", 
                            shell=True)
                print(f"aligned, writing {system_path}/T{trj_index}/{selection_name}.xtc")


@uf.get_timing
def align_selections_parallel():
    pool = Pool()
    pool.map(align_selections, trj_indices)


trj_indices = range(1, sd.N_TRJ+1)
align_selections_parallel()
