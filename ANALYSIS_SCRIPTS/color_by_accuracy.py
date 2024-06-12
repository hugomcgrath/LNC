import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf
from glob import glob


def extract_results_file_config_id(result_path):
    return int(result_path.split("_")[-1].split(".")[0])


base_accuracy = 0.5
# RESULTS = f"{sd.BASE_DIR}/RESULTS_2024-6-12_9-30-20"
RESULTS = sd.RESULTS
pdb = f"{sd.SYSTEMS['LNC']}/common_atoms.pdb"
universe = mda.Universe(pdb)
pdb_selection_keywords = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
with open(f"{sd.SYSTEMS['LNC']}/config_indices_for_visualization.txt", "w") as file:
    result_paths = sorted(glob(f"{RESULTS}/*passthrough*validation*"), key=extract_results_file_config_id)
    for result_path in result_paths:
        cf = uf.get_config_file_from_file_path(result_path)
        accuracy_mean, _ = uf.get_accuracy_mean_and_sem(result_path)
        selection_name = cf.sel_name
        selection_keyword = sd.SELECTION_DICT[selection_name]
        selection_keyword_indices = uf.get_selection_keyword_indices(
            pdb_selection_keywords, 
            selection_name, 
            selection_keyword,
        )
        selection = universe.select_atoms(selection_keyword_indices)
        print(f"Setting tempfactors and occupancies for selection: {selection_name}")
        file.write(f"{selection_name}: {cf.config_index}\n")
        selection.tempfactors = [accuracy_mean for _ in range(selection.n_atoms)]
        selection.occupancies = [cf.config_index for _ in range(selection.n_atoms)]
selection_not_classified = universe.select_atoms("prop tempfactor == 0")
print(f"Setting tempfactors and occupancies for parts not used in classification")
selection_not_classified.tempfactors = [base_accuracy for _ in range(selection_not_classified.n_atoms)]
selection_not_classified.occupancies = [999 for _ in range(selection_not_classified.n_atoms)]
selection_all = universe.select_atoms("all")
print("Writing whole selection with tempfactors set")
selection_all.write(f"{sd.SYSTEMS['LNC']}/sphere_passthrough_linear_colored_by_accuracy.pdb")
