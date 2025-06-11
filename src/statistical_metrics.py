"""Statistical metrics to compare physical metrics on truth and predictions."""

import numpy as np


def compute_mae(predictions: dict, targets: dict) -> float:
    """
    Compute Mean Absolute Error (MAE) between predicted and true metric values.

    Parameters
    ----------
    predictions : dict
        Dictionary of predicted metric values.
    targets : dict
        Dictionary of true metric values.

    Returns
    -------
    float
        The Mean Absolute Error.
    """
    abs_errors = [
        abs(float(predictions[key]) - float(targets[key]))
        for key in predictions
        if key in targets
    ]
    return np.mean(abs_errors)


def compute_rmse(predictions: dict, targets: dict) -> float:
    """
    Compute Root Mean Squared Error (RMSE) between predicted and true metric values.

    Parameters
    ----------
    predictions : dict
        Dictionary of predicted metric values.
    targets : dict
        Dictionary of true metric values.

    Returns
    -------
    float
        The Root Mean Squared Error.
    """
    squared_errors = [
        (float(predictions[key]) - float(targets[key])) ** 2
        for key in predictions
        if key in targets
    ]
    return np.sqrt(np.mean(squared_errors))
