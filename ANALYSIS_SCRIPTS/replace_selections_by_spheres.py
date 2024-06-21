import shared_data as sd
import MDAnalysis as mda
import utility_functions as uf
from glob import glob


def extract_results_file_config_id(result_path):
    return int(result_path.split("_")[-1].split(".")[0])


def replace_selections_by_spheres_and_color_by_results(result_type, selection_type, dim_reduction_name, model_name, feature_type, partition_name):
    name = f"{selection_type}_{dim_reduction_name}_{model_name}_{feature_type}_{partition_name}_{result_type}_spheres"
    RESULTS = sd.RESULTS
    pdb = f"{sd.SYSTEMS['LNC']}/added_antibiotic.pdb"
    universe = mda.Universe(pdb)
    result_paths = sorted(glob(f"{RESULTS}/*{dim_reduction_name}_{model_name}_{feature_type}_{partition_name}*"), key=extract_results_file_config_id)
    pdb_spheres = f"{sd.SYSTEMS['LNC']}/{name}.pdb"
    with open(pdb_spheres, "w") as file_out:
        for result_path in result_paths:
            cf = uf.get_config_file_from_file_path(result_path)
            accuracy_mean, _ = uf.get_accuracy_mean_and_sem(result_path)
            tempfactor = f"{accuracy_mean:.2f}".rjust(6)
            selection_name = cf.sel_name
            with open(f"{sd.SYSTEMS['LNC']}/{selection_name}.pdb") as file_in:
                lines = file_in.readlines()
                line = lines[2]
                beginning = line[:30]
                middle = line[55:61]
                end = line[67:]
                # print(line[30:54])
            selection_keyword = sd.SELECTION_DICT[selection_name]
            selection_keyword_indices = uf.get_selection_keyword_indices(
                pdb, 
                selection_name, 
                selection_keyword,
            )
            selection = universe.select_atoms(selection_keyword_indices)
            sphere_coordinates = selection.center_of_geometry()
            x = f"{sphere_coordinates[0]:.3f}".rjust(8)
            y = f"{sphere_coordinates[1]:.3f}".rjust(8)
            z = f"{sphere_coordinates[2]:.3f}".rjust(8)
            coords = x + y + z
            line_out = beginning + coords + middle + tempfactor + end
            file_out.write(line_out)

result_type = "accuracy"
selection_type = "grid"
dim_reduction_name = "passthrough"
model_name = "linear"
feature_type = "cartesian"
partition_name = "validation"

replace_selections_by_spheres_and_color_by_results(
    result_type, 
    selection_type, 
    dim_reduction_name, 
    model_name, 
    feature_type, 
    partition_name
)
