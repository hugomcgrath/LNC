import shared_data as sd
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj

# find intersections
selections = {}
for system_name, system_path in sd.SYSTEMS.items():
    pdb = f"{system_path}/system.pdb"
    universe = mda.Universe(pdb)
    selection = universe.select_atoms("all")
    selections[system_name] = set([(resid, resname, segid, name) for resid, resname, segid, name in zip(selection.resids, 
                                                                                                        selection.resnames, 
                                                                                                        selection.segids, 
                                                                                                        selection.names)])
intersection = list(selections['LNC'].intersection(selections['NONE']))

for system_name, system_path in sd.SYSTEMS.items():
    with open(f"{system_path}/system.pdb", "r") as file_old, open(f"{system_path}/common_atoms.pdb", "w") as file_new:
        lines_old = file_old.readlines()
        for line_old in lines_old:
            if "ATOM" not in line_old:
                file_new.write(line_old)
            else:
                resid = int(line_old[22:26].strip())
                resname = line_old[17:20].strip()
                segid = line_old[21].strip()
                name = line_old[12:16].strip()
                identifier = (resid, resname, segid, name)
                if identifier in intersection:
                    file_new.write(line_old)
    print(f"wrote {system_path}/common_atoms.pdb")

# write trajectories with common atoms
for system_name, system_path in sd.SYSTEMS.items():
    pdb_common_path = f"{system_path}/common_atoms.pdb"
    selection_string = "index"
    with open(pdb_common_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            if "ATOM" in line:
                index = line[6:11].strip()
                selection_string += f" {index}"
    for trj_index in range(1, sd.N_TRJ+1):
        universe = mda.Universe(f"{system_path}/system.pdb", 
                                f"{system_path}/T{trj_index}/traj_pruned.xtc")
        selection = universe.select_atoms(selection_string)
        print(f"writing {system_path}/T{trj_index}/common_atoms.xtc")
        selection.write(f"{system_path}/T{trj_index}/common_atoms.xtc",
                        frames="all")
        print(f"wrote {system_path}/T{trj_index}/common_atoms.xtc")