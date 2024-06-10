import shared_data as sd
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import subprocess


for system_name, system_path in sd.SYSTEMS.items():
    for trj_index in range(1, sd.N_TRJ+1):
        for selection_name in sd.SELECTION_DICT:
            mobile = mda.Universe(f"{system_path}/{selection_name}.pdb",
                                  f"{system_path}/T{trj_index}/{selection_name}.xtc")
            # only one common reference
            reference = mda.Universe(f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb")
            # align to atoms of the selection
            alignment_selection_string = f"({sd.ALIGNMENT_SELECTION})"
            print(f"aligning {system_path}/T{trj_index}/{selection_name}.xtc")
            AlignTraj(mobile, reference, alignment_selection_string, 
                      filename=f"{system_path}/T{trj_index}/{selection_name}_aligned.xtc").run()
            subprocess.run(f"mv {system_path}/T{trj_index}/{selection_name}_aligned.xtc {system_path}/T{trj_index}/{selection_name}.xtc", 
                           shell=True)
            print(f"aligned, writing {system_path}/T{trj_index}/{selection_name}.xtc")
