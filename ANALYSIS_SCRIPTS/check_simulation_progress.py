import shared_data as sd
import os
import subprocess
import matplotlib.pyplot as plt
import MDAnalysis as mda


subprocess.run([f"bash ~/LNC/ANALYSIS_SCRIPTS/copy_trj.sh"], shell=True)
plt.rcParams.update({'font.size': 18})
fig, axs = plt.subplots(1, 2, sharey=True)
fig.subplots_adjust(wspace=0)
fig.set_size_inches(18, 6)
for i, (system_name, system_path) in enumerate(sd.SYSTEMS.items()):
    trj_size = []
    x = [_ for _ in range(1, sd.N_TRJ+1)]
    pdb = f"{system_path}/system.pdb"
    for j in x:
        xtc = f"{system_path}/T{j}/traj_comp.xtc"
        # if os.path.exists(xtc):
        #     result = subprocess.run([f"du -h {xtc}"], shell=True, stdout=subprocess.PIPE, text=True)
        #     size_as_str = result.stdout.split("\t")[0]
        #     if "G" in size_as_str:
        #         trj_size.append(float(size_as_str.removesuffix("G")))
        #     elif "M" in size_as_str:
        #         trj_size.append(float(size_as_str.removesuffix("M")) / 1000)
        # else:
        #     trj_size.append(0)
        u = mda.Universe(pdb, xtc)
        print(f"{xtc}: {len(u.trajectory)}")
        trj_size.append(len(u.trajectory))
    axs[i].bar(x, trj_size)
    axs[i].set_xticks([_ for _ in range(0, sd.N_TRJ+1, 2)])
    axs[i].set_xlabel("Trajectory index")
    axs[i].set_title(system_name)
# axs[0].set_ylabel("Trajectory file size (GB)")
axs[0].set_ylabel("# of frames")
plt.show()