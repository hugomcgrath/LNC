import MDAnalysis as mda
import numpy as np
import shared_data as sd
import ml_models as mm
import utility_functions as uf
import sys
import timeit
from config_file_utils import ConfigFile
from datetime import timedelta
import tracemalloc


def classify_selections(config_file_path, config_file_name):

    cf = ConfigFile(f"{config_file_path}/{config_file_name}")

    pair_dir = f"{sd.BASE_DIR}/{cf.pair_name}"
    pdb = f"{pair_dir}/{cf.sel_name}.pdb"

    y = uf.get_y(cf)
    y_predicted_array = {"train": np.zeros((sd.N_CROSS_VALIDATION_COMBINATIONS, len(y["train"]))),
                     cf.partition_name: np.zeros((sd.N_CROSS_VALIDATION_COMBINATIONS, len(y[cf.partition_name])))}

    for cv_index in range(sd.N_CROSS_VALIDATION_COMBINATIONS):
        uni = {
            "train": mda.Universe(pdb, f"{pair_dir}/{cf.sel_name}_train_{cv_index}.xtc"),
            cf.partition_name: mda.Universe(pdb, f"{pair_dir}/{cf.sel_name}_{cf.partition_name}_{cv_index}.xtc")
        }
        X = {
            "train": uf.get_X(uni["train"], cf.config_dict), 
            cf.partition_name: uf.get_X(uni[cf.partition_name], cf.config_dict)
        }
        y_predicted = {}
        y_predicted["train"], y_predicted[cf.partition_name] = mm.ml_pipeline(
            X["train"],
            X[cf.partition_name],
            y["train"],
            cf.config_dict
        )
        y_predicted_array["train"][cv_index, :] = y_predicted["train"].reshape(len(y["train"]))
        y_predicted_array[cf.partition_name][cv_index, :] = y_predicted[cf.partition_name].reshape(len(y[cf.partition_name]))

    result_file_name = {"train": f"{sd.RESULTS}/y_{cf.pair_name}_{cf.sel_name}_{cf.dim_reduction_name}_{cf.model_name}_{cf.feature_type}_train_{cf.config_index}.npy",
                        cf.partition_name: f"{sd.RESULTS}/y_{cf.pair_name}_{cf.sel_name}_{cf.dim_reduction_name}_{cf.model_name}_{cf.feature_type}_{cf.partition_name}_{cf.config_index}.npy"}
    print(f"Writing {result_file_name['train']}")
    np.save(result_file_name["train"], y_predicted_array["train"])
    print(f"Writing {result_file_name[cf.partition_name]}")
    np.save(result_file_name[cf.partition_name], y_predicted_array[cf.partition_name])

if __name__ == "__main__":
    start_time = timeit.default_timer()
    
    tracemalloc.start()
    config_file_path = sys.argv[1]
    config_file_name = sys.argv[2]
    classify_selections(config_file_path, config_file_name)
    _, peak = tracemalloc.get_traced_memory()
    print(f"Peak memory usage: {peak / 10**9:.2f}GB")
    tracemalloc.stop()

    total_time = timeit.default_timer() - start_time
    print(f"Time {timedelta(seconds=total_time)}")
