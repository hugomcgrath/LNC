import shared_data as sd
import yaml

class ConfigFile:

    def __init__(self, config_file_path):
        self.config_file_path = config_file_path
        self.unpack_config_file()

    def unpack_config_file(self):
        with open(self.config_file_path, "r") as config_file:
            self.config_dict = yaml.safe_load(config_file)
        self.pair_name = self.config_dict["pair_name"]
        self.sel_name = self.config_dict["selection_name"]
        self.dim_reduction_name = list(self.config_dict["dim_reduction"].items())[0][0]
        self.model_name = list(self.config_dict["model"].items())[0][0]
        self.feature_type = list(self.config_dict["feature_type"].items())[0][0]
        self.partition_name = self.config_dict["partition_name"]
        self.config_index = self.config_dict["config_index"]

