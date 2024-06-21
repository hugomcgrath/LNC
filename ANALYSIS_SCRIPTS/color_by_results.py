import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf
from glob import glob


def extract_results_file_config_id(result_path):
    return int(result_path.split("_")[-1].split(".")[0])


def color_ribosome_by_results(result_type, selection_type, dim_reduction_name, model_name, feature_type, partition_name):
    name = f"{selection_type}_{dim_reduction_name}_{model_name}_{feature_type}_{partition_name}_{result_type}"
    RESULTS = sd.RESULTS
    pdb = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
    universe = mda.Universe(pdb)
    with open(f"{sd.SYSTEMS['LNC']}/{name}_config_indices.txt", "w") as file:
        result_paths = sorted(glob(f"{RESULTS}/*{dim_reduction_name}_{model_name}_{feature_type}_{partition_name}*"), key=extract_results_file_config_id)
        for result_path in result_paths:
            cf = uf.get_config_file_from_file_path(result_path)
            accuracy_mean, accuracy_sem = uf.get_accuracy_mean_and_sem(result_path)
            selection_name = cf.sel_name
            selection_keyword = sd.SELECTION_DICT[selection_name]
            selection_keyword_indices = uf.get_selection_keyword_indices(
                pdb, 
                selection_name, 
                selection_keyword,
            )
            selection = universe.select_atoms(selection_keyword_indices)
            print(f"Setting tempfactors and occupancies for selection: {selection_name}")
            file.write(f"{selection_name}: {cf.config_index}\n")
            if result_type == "accuracy":
                result = accuracy_mean
                base_value = 0
            elif result_type == "sem":
                result = accuracy_sem
                base_value = -1
            elif result_type == "sd":
                result = accuracy_sem * sd.N_CROSS_VALIDATION_COMBINATIONS ** (1 / 2)
                base_value = -1
            selection.tempfactors = [result for _ in range(selection.n_atoms)]
            selection.occupancies = [cf.config_index for _ in range(selection.n_atoms)]
    selection_not_classified = universe.select_atoms("prop tempfactor == 0")
    print(f"Setting tempfactors and occupancies for parts not used in classification")
    selection_not_classified.tempfactors = [base_value for _ in range(selection_not_classified.n_atoms)]
    selection_not_classified.occupancies = [999 for _ in range(selection_not_classified.n_atoms)]
    selection_all = universe.select_atoms("all")
    print("Writing whole selection with tempfactors set")
    selection_all.write(f"{sd.SYSTEMS['LNC']}/{name}.pdb")


result_type = "accuracy"
selection_type = "grid"
dim_reduction_name = "passthrough"
model_name = "linear"
feature_type = "cartesian"
partition_name = "validation"

color_ribosome_by_results(
    result_type, 
    selection_type, 
    dim_reduction_name, 
    model_name, 
    feature_type, 
    partition_name
)
