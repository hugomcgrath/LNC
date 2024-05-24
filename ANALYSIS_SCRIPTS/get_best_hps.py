import shared_data as sd
from glob import glob
import yaml
import utility_functions as uf
from config_file_utils import ConfigFile
import pandas as pd

data = []
config_indices = []
for y_predicted_array_name in glob(f"{sd.RESULTS}/*validation*"):
    cf = uf.get_config_file_from_file_path(y_predicted_array_name)
    if cf.partition_name == "validation":
        y = uf.get_y(cf)
        accuracy_mean, _ = uf.get_accuracy_mean_and_sem(y_predicted_array_name)
        data_point = {}
        data_point["pair_name"] = cf.pair_name
        data_point["sel_name"] = cf.sel_name
        data_point["dim_reduction_name"] = cf.dim_reduction_name
        data_point["model_name"] = cf.model_name
        data_point["feature_type"] = cf.feature_type
        data_point["accuracy_mean"] = accuracy_mean
        data.append(data_point)
        config_indices.append(cf.config_index)
df = pd.DataFrame(data, index=config_indices)
config_indices_best = df.groupby(["pair_name", "sel_name"]).idxmax()["accuracy_mean"]

for i, config_index_best in enumerate(config_indices_best):
    config_file_path = f"{sd.CONFIG_FILES_DIR}/config_{config_index_best}.yaml"
    config_file_path_best = f"{sd.CONFIG_FILES_BEST_DIR}/config_{i}.yaml"
    cf = ConfigFile(config_file_path)
    cf.config_dict["config_index"] = config_index_best
    cf.config_dict["partition_name"] = "test"
    with open(config_file_path_best, "w") as config_file_best:
        yaml.dump(cf.config_dict, config_file_best)
