import MDAnalysis as mda
import shared_data as sd


for system_name, system_path in sd.SYSTEMS.items():
    pdb = f"{system_path}/system.pdb"
    for i in range(1, sd.N_TRJ+1):
        xtc = f"{system_path}/T{i}/traj_comp.xtc"
        try:
            u = mda.Universe(pdb, xtc)
            print(f"{xtc}: {len(u.trajectory)} frames")
        except:
            print(f"{xtc} missing")