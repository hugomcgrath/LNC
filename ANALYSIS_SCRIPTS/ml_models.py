import shared_data as sd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from pyemma.coordinates.transform import TICA


def ml_pipeline(X_train, X_test, y_train, config_dict):

    dim_reduction_name, dim_reduction_hps = list(config_dict["dim_reduction"].items())[0]
    if dim_reduction_name == "passthrough":
        dim_reduction = FunctionTransformer(**dim_reduction_hps)
    elif dim_reduction_name == "pca":
        dim_reduction = PCA(**dim_reduction_hps)
    elif dim_reduction_name == "tica":
        dim_reduction = TICA(**dim_reduction_hps)

    model_name, model_hps = list(config_dict["model"].items())[0]
    if model_name == "linear":
        model = LinearRegression(**model_hps)
    elif model_name == "svc_linear":
        model = LinearSVC(**model_hps)
    elif model_name == "random_forest":
        model = RandomForestClassifier(**model_hps)

    steps=[
        ("scaler", StandardScaler()),
        ("dim_reduction", dim_reduction),
        ("model", model),
    ]
    pipeline = Pipeline(steps)
    pipeline.fit(X_train, y_train)
    y_predicted_train = pipeline.predict(X_train)
    y_predicted_test = pipeline.predict(X_test)

    if model_name == "linear":
        y_predicted_train = apply_threshold(y_predicted_train)
        y_predicted_test = apply_threshold(y_predicted_test)

    return y_predicted_train, y_predicted_test


def apply_threshold(y, threshold=0.5):
    return np.array(y >= threshold, dtype=int)
