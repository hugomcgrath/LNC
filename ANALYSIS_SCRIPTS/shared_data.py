N_TRJ = 20
# pruned individual trajectory length
TRJ_LEN = ...
N_CROSS_VALIDATION_COMBINATIONS = 10

BASE_DIR = "/home/hmcgrat/LNC"
RESULTS = f"{BASE_DIR}/RESULTS"
CONFIG_FILES_DIR = f"{BASE_DIR}/HP_CONFIG_FILES"
CONFIG_FILES_BEST_DIR = f"{BASE_DIR}/HP_CONFIG_FILES_BEST"

SYSTEMS = {
    "LNC": f"{BASE_DIR}/LNC",
    "NONE": f"{BASE_DIR}/NONE",
}
PAIRS = {
    "LNC_NONE": {"LNC": SYSTEMS["LNC"],
                 "NONE": SYSTEMS["NONE"]},
}

ALIGNMENT_SELECTION = "name CA P"
HEAVY_ATOMS_SELECTION_KEYWORD = "type C O N S P"
NUCLEIC_PRUNED_SELECTION_KEYWORD = "(nucleic and (nucleicbackbone or name C1' C4' N1 N3 C5 C8))"
REFERENCE_POINT_FOR_SPHERE = "segid 9"
SELECTION_DICT = {
    ...
}

DIM_REDUCTION_HPS = {
    "passthrough": {},
    "pca": {"n_components": [0.05 + i * 0.05 for i in range(19)]},
    # tica doesn't seem to affect accuracy much
    # "tica": {"var_cutoff": [0.95, 0.9, 0.85, 0.8, 0.75], "lag": [1, 2, 5, 10, 25, 50]},
}
MODEL_HPS = {
    "linear": {},
    # "svc_linear": {},
    # "random_forest": {"n_estimators": [50, 100, 250]},
}
FEATURE_HPS = {
    "cartesian": {},
    # "cartesian_averaged": {"window_width": [10], "overlapping": [False]},
    # "rmsf": {"window_width": [10], "overlapping": [False]},
}
PARTITION_PERCENTAGES = {
    "train": 0.6,
    "validation": 0.2,
    "test": 0.2,
}
PARTITION_N_TRJS = {
    partition_name: int(partition_percentage * N_TRJ) for (partition_name, partition_percentage) in PARTITION_PERCENTAGES.items()
}
