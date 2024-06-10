import numpy as np
import shared_data as sd
from config_file_utils import ConfigFile
from scipy.stats import sem
import MDAnalysis as mda
from MDAnalysis.analysis.bat import BAT


def extract_trj(aligned_u, selection_keyword="all"):
    """
    Extract the trajectory from the universe object into a numpy array
    and calculate the mean structure of this trajectory.

    Parameters:

        aligned_u (object): Universe object created from the
        aligned trajectories.

    Returns:

        X (array (n_frames, 3 * n_atoms)): Trajectory extracted into
        an array.
    """
    trj = aligned_u.trajectory
    atomgroup = aligned_u.select_atoms(selection_keyword)
    n_frames = len(trj)
    n_atoms = len(atomgroup)
    # initialize empty array
    X = np.zeros((n_frames, 3 * n_atoms))
    # fill array, for every frame insert structure coordinates
    for frame, _ in enumerate(trj):
        if sd.NORMALIZE_BY_AVERAGE:
            x = atomgroup.positions
            average_norm = np.mean(np.linalg.norm(x, axis=1))
            X[frame, :] = (x / average_norm).ravel()
        else:
            X[frame, :] = atomgroup.positions.ravel()
    return X


def get_rmsf_features(X, window_width, overlapping):
    """
    Assign to each coordinate in each frame fluctuations, similar to RMSF,
    calculated from the surrounding frames (something like a moving average).
    """
    if overlapping:
        rmsf_features = np.zeros(X.shape)
        for i in range(X.shape[0]):
            window = X[min(i-window_width, 0):max(i+window_width,
                                                  X.shape[0]), :]
            rmsf_features[i] = np.std(window, axis=0)
    else:
        rmsf_features = np.zeros((X.shape[0]//window_width, X.shape[1]))
        for i in range(rmsf_features.shape[0]):
            window = X[i*window_width:(i+1)*window_width, :]
            rmsf_features[i] = np.std(window, axis=0)
    return rmsf_features


def get_cartesian_averaged_features(X, window_width, overlapping):
    if overlapping:
        cartesian_averaged_features = np.zeros(X.shape)
        for i in range(X.shape[0]):
            window = X[min(i-window_width, 0):max(i+window_width,
                                                  X.shape[0]), :]
            cartesian_averaged_features[i] = np.mean(window, axis=0)
    else:
        cartesian_averaged_features = np.zeros((X.shape[0]//window_width,
                                                X.shape[1]))
        for i in range(cartesian_averaged_features.shape[0]):
            window = X[i*window_width:(i+1)*window_width, :]
            cartesian_averaged_features[i] = np.mean(window, axis=0)
    return cartesian_averaged_features


def get_X(aligned_u, config_dict):
    feature_type, feature_hps = list(config_dict["feature_type"].items())[0]
    X = extract_trj(aligned_u)
    if feature_type == "cartesian":
        pass
    elif feature_type == "cartesian_averaged":
        X = get_cartesian_averaged_features(X, **feature_hps)
    elif feature_type == "rmsf":
        X = get_rmsf_features(X, **feature_hps)
    # missing bonds, doesn't work
    # elif feature_type == "bat":
    #     bat_object = BAT(aligned_u)
    #     bat_object.run()
    #     X = bat_object.results.bat
    return X


def get_y(config_file):
    feature_type, feature_hps = list(config_file.config_dict["feature_type"].items())[0]
    y = {}
    for partition_name, partition_n_trj in sd.PARTITION_N_TRJS.items():
        y[partition_name] = np.zeros(
            partition_n_trj * sd.TRJ_LEN * len(sd.PAIRS[config_file.pair_name])
        )
        if sd.WRONG_LABELS:
            y[partition_name][np.random.choice(range(partition_n_trj*sd.TRJ_LEN), size=partition_n_trj*sd.TRJ_LEN//2, replace=False)] = 1
            y[partition_name][np.random.choice(
                range(partition_n_trj*sd.TRJ_LEN, partition_n_trj*sd.TRJ_LEN*len(sd.PAIRS[config_file.pair_name])), 
                size=partition_n_trj*sd.TRJ_LEN//2, 
                replace=False
            )] = 1
        else:
            y[partition_name][:partition_n_trj*sd.TRJ_LEN] = 1

        if (feature_type == "cartesian") or (feature_type == "bat"):
            pass
        elif feature_type == "rmsf" or feature_type == "cartesian_averaged":
            if not feature_hps["overlapping"]:
                y[partition_name] = y[partition_name][::feature_hps["window_width"]]
        elif feature_type == "combined":
            if feature_hps["overlapping"]:
                y[partition_name] = np.concatenate((y[partition_name], y[partition_name]), axis=0)
            else:
                y = np.concatenate(
                    (y[partition_name], y[partition_name][::feature_hps["window_width"]]),
                    axis=0
                )
    return y


def get_combination_indices():
    combination_indices = {}
    for partition_name, partition_n_trj in sd.PARTITION_N_TRJS.items():
        combination_indices[partition_name] = np.zeros(
            (sd.N_CROSS_VALIDATION_COMBINATIONS, partition_n_trj), 
            dtype=int
        ) 
    for i in range(sd.N_CROSS_VALIDATION_COMBINATIONS):
        trj_indices = [_ for _ in range(1, sd.N_TRJ+1)]
        for partition_name, partition_n_trj in sd.PARTITION_N_TRJS.items():
            combination_indices[partition_name][i] = np.random.choice(trj_indices, 
                                                                      partition_n_trj, 
                                                                      replace=False)
            trj_indices = np.setdiff1d(trj_indices, combination_indices[partition_name][i])
    return combination_indices


def get_config_file_from_file_path(file_path, type_of_file="result"):
    if type_of_file == "result":
        config_id = file_path.split("_")[-1].strip(".npy")
    elif type_of_file == "log":
        config_id = int(file_path.split(".")[-1]) - 1
    cf = ConfigFile(f"{sd.CONFIG_FILES_DIR}/config_{config_id}.yaml")
    return cf


def get_accuracy_mean_and_sem(result_file_path):
    cf = get_config_file_from_file_path(result_file_path)
    y_predicted_array = np.load(result_file_path)
    y = get_y(cf)
    boolean_array = np.logical_not(np.logical_xor(y_predicted_array, y[cf.partition_name]))
    accuracy_array = np.sum(boolean_array, axis=1) / boolean_array.shape[1]
    accuracy_mean = np.mean(accuracy_array)
    accuracy_sem = sem(accuracy_array)
    return accuracy_mean, accuracy_sem


# make sure not to accidentaly run twice or relabel something correctly labeled
def relabel_chains(pdb_in, pdb_out, chain_map):
    with open(pdb_in, "r") as pdb_file_in, open(pdb_out, "w") as pdb_file_out:
        pdb_lines_in = pdb_file_in.readlines()
        for pdb_line_in in pdb_lines_in:
            if len(pdb_line_in) == 79 and pdb_line_in[21] in chain_map:
                pdb_line_in = list(pdb_line_in)
                pdb_line_in[21] = chain_map[pdb_line_in[21]]
                pdb_line_out = "".join(pdb_line_in)
                pdb_file_out.write(f"{pdb_line_out}")
            else:
                pdb_file_out.write(f"{pdb_line_in}")


def get_selection_keyword_indices(pdb, selection_name, selection_keyword):
    universe = mda.Universe(pdb)
    reference_point = universe.select_atoms(sd.REFERENCE_POINT).center_of_geometry()
    selection_keyword += f" and (({sd.HEAVY_ATOMS_SELECTION_KEYWORD}) or {sd.NUCLEIC_PRUNED_SELECTION_KEYWORD} and (not {sd.REFERENCE_POINT}))"
    segment_selection_dict = {
        "e1": f" and (prop x > {reference_point[0]} and prop y > {reference_point[1]} and prop z > {reference_point[2]})",
        "e2": f" and (prop x < {reference_point[0]} and prop y > {reference_point[1]} and prop z > {reference_point[2]})",
        "e3": f" and (prop x > {reference_point[0]} and prop y < {reference_point[1]} and prop z > {reference_point[2]})",
        "e4": f" and (prop x > {reference_point[0]} and prop y > {reference_point[1]} and prop z < {reference_point[2]})",
        "e5": f" and (prop x < {reference_point[0]} and prop y < {reference_point[1]} and prop z > {reference_point[2]})",
        "e6": f" and (prop x < {reference_point[0]} and prop y > {reference_point[1]} and prop z < {reference_point[2]})",
        "e7": f" and (prop x > {reference_point[0]} and prop y < {reference_point[1]} and prop z < {reference_point[2]})",
        "e8": f" and (prop x < {reference_point[0]} and prop y < {reference_point[1]} and prop z < {reference_point[2]})",
    }
    # if "grid" in selection_name:
    #     grid_point = np.array([float(grid_point_coordinate) for grid_point_coordinate in selection_name.split("_")[-1].split(":")])
    #     grid_point += reference_point
    #     selection_keyword += f" and point {grid_point[0]} {grid_point[1]} {grid_point[2]} {sd.GRID_POINT_RADIUS}"
    for segment_name, segment_selection_keyword in segment_selection_dict.items():
        if segment_name in selection_name:
            selection_keyword += segment_selection_keyword

    selection = universe.select_atoms(selection_keyword, updating=False)
    print(selection.center_of_geometry())
    selection_keyword_indices = "index"
    for index in selection.indices:
        selection_keyword_indices += f" {index}"
    return selection_keyword_indices


def get_time_windowed_data(data, time_window_index):
    time_windowed_data = {}
    for partition_name in data:
        n_subarrays = sd.N_TIME_WINDOWS * sd.PARTITION_N_TRJS[partition_name] * len(sd.PAIRS["LNC_NONE"])
        subarray_list = np.split(data[partition_name], n_subarrays, axis=0)
        time_windowed_data[partition_name] = np.concatenate(subarray_list[time_window_index::sd.N_TIME_WINDOWS], axis=0)
    return time_windowed_data
