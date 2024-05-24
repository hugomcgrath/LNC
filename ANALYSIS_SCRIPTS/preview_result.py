import sys
import numpy as np
import shared_data as sd
from config_file_utils import ConfigFile
import utility_functions as uf

result_path = sys.argv[1]
partition_name = result_path.split("_")[-2]
config_id = result_path.split("/")[-1].split("_")[-1].strip(".npy")
y_predicted = np.load(result_path)
cf = ConfigFile(f"{sd.CONFIG_FILES_DIR}/config_{config_id}.yaml")
y = uf.get_y(cf)
accuracy_mean, accuracy_sem = uf.get_accuracy_mean_and_sem(y_predicted, y, partition_name, 5)

print(f"result array size: {y_predicted.shape}")
print("result array:")
print(y_predicted)
print("cross validation accuracy:")
print(f"{accuracy_mean}\u00B1{accuracy_sem}")
print(np.sum(y_predicted, axis=1))
