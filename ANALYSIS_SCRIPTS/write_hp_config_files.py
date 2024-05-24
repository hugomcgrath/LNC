import yaml
import itertools
import os
import shared_data as sd
from datetime import datetime


def backup():
    ct = datetime.now()
    time_identifier = f"{ct.year}-{ct.month}-{ct.day}_{ct.hour}-{ct.minute}-{ct.second}"
    os.system(f"cp -r {sd.CONFIG_FILES_DIR} {sd.CONFIG_FILES_DIR}_{time_identifier}")
    os.system(f"cp -r {sd.CONFIG_FILES_BEST_DIR} {sd.CONFIG_FILES_BEST_DIR}_{time_identifier}")
    os.system(f"cp -r {sd.RESULTS} {sd.RESULTS}_{time_identifier}")


def generate_configs(step_name, hps_dict):
    configs = []
    for name, hp_dict in hps_dict.items():
        hp_values_list = []
        for hp_name, hp_values in hp_dict.items():
            hp_values_list.append(hp_values)
        hp_values_combinations = itertools.product(*hp_values_list)
        for combination in hp_values_combinations:
            config_dict = {name: {}}
            for hp_name, hp_value in zip(hp_dict, combination):
                config_dict[name][hp_name] = hp_value
            configs.append({step_name: config_dict})
    return configs

user_input = input("Backup results and config files? ([y]/n): ")
if user_input != "n":
    backup()
    print("Backed up")
else:
    print("Not backing up")

user_input = input("Rewrite old results and config files? ([y]/n): ")
if user_input != "n":
    os.system(f"rm {sd.CONFIG_FILES_DIR}/* {sd.CONFIG_FILES_BEST_DIR}/* {sd.RESULTS}/*")
    print("Rewriting")
else:
    print("Not rewriting")
config_index_offset = len(os.listdir(sd.CONFIG_FILES_DIR))

configs_pairs = [{"pair_name": pair_name} for pair_name in sd.PAIRS]
configs_selections = [{"selection_name": selection_name} for selection_name in sd.SELECTION_DICT]
configs_dim_reduction = generate_configs("dim_reduction", sd.DIM_REDUCTION_HPS)
configs_models = generate_configs("model", sd.MODEL_HPS)
config_features = generate_configs("feature_type", sd.FEATURE_HPS)

config_combinations = itertools.product(*[
    configs_pairs, 
    configs_selections, 
    configs_dim_reduction, 
    configs_models, 
    config_features
])
for i, config_combination in enumerate(config_combinations):
    pair_name, selection_name, dim_reduction, model, feature_type = config_combination
    with open(f"{sd.CONFIG_FILES_DIR}/config_{i+config_index_offset}.yaml", "w") as config_file:
        config_file.write("---\n")
        yaml.dump({"config_index": i+config_index_offset}, config_file)
        yaml.dump(pair_name, config_file)
        yaml.dump(selection_name, config_file)
        yaml.dump(dim_reduction, config_file)
        yaml.dump(model, config_file)
        yaml.dump(feature_type, config_file)
        yaml.dump({"partition_name": "validation"}, config_file)
print("Don't forget to modify run file!")
print(f"#$ -t {1+config_index_offset}-{len(os.listdir(sd.CONFIG_FILES_DIR))}")