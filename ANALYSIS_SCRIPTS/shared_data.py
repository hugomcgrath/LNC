from itertools import product


N_TRJ = 20
# pruned individual trajectory length
TRJ_LEN = 700
N_CROSS_VALIDATION_COMBINATIONS = 10
# setting to 1 means no time windowing
N_TIME_WINDOWS = 1
# if set to True, half of the labels are set to incorrect values
# this should ideally result in accuracy close to 0.5
WRONG_LABELS = True

BASE_DIR = "/home/hmcgrat/LNC"
RESULTS = f"{BASE_DIR}/RESULTS"
LOGS = f"{BASE_DIR}/LOGS"
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

PREALIGNMENT_SELECTION = "name CA P"
ALIGNMENT_SELECTION = "name CA P"
HEAVY_ATOMS_SELECTION_KEYWORD = "type C O N S P"
NUCLEIC_PRUNED_SELECTION_KEYWORD = "(nucleic and (nucleicbackbone or name C1' C4' N1 N3 C5 C8))"
REFERENCE_POINT = "segid 9"

# # selections around grid points
# GRID_SPACING = 30
# GRID_POINT_RADIUS = 15
# GRID_RANGE_X = range(-3*GRID_SPACING, 4*GRID_SPACING, GRID_SPACING)
# GRID_RANGE_Y = range(-3*GRID_SPACING, 4*GRID_SPACING, GRID_SPACING)
# GRID_RANGE_Z = range(-2*GRID_SPACING, 3*GRID_SPACING, GRID_SPACING)
# GRID_POINTS = product(GRID_RANGE_X, GRID_RANGE_Y, GRID_RANGE_Z)
# SELECTION_DICT = {
#     f"grid_{i}:{j}:{k}": "all" for i, j, k in GRID_POINTS
# }

# spherical layers
R = [15, 25, 31.50, 36.06, 39.69, 42.75, 45.43, 47.82, 50]
SELECTION_DICT = {
    f"sphere_layer{i+1}": f"byres sphlayer {R[i]} {R[i+1]} {REFERENCE_POINT}" for i, _ in enumerate(R[:-1])
}

DIM_REDUCTION_HPS = {
    "passthrough": {},
    # "pca": {"n_components": [0.05 + i * 0.05 for i in range(19)]},
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
    # "bat": {},
}
PARTITION_PERCENTAGES = {
    "train": 0.6,
    "validation": 0.2,
    "test": 0.2,
}
PARTITION_N_TRJS = {
    partition_name: int(partition_percentage * N_TRJ) for (partition_name, partition_percentage) in PARTITION_PERCENTAGES.items()
}
