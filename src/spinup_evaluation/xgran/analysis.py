"""
Perform resampling analysis prior to data processing.

This module contains all the analysis functions for determining what variables
and metrics can run at different temporal granularities in climate data processing.

Functions:
- analyze_variable_availability: Core function to analyze variable availability
- analyze_metric_requirements: Analyze what metrics can run at what granularities
- get_runnable_metrics_at_granularity: Get metrics that can run at specific granularity
- get_optimal_granularities_for_metrics: Find finest granularity for each metric
- show_availability_summary: Comprehensive availability report
- get_maximum_granularity_with_all: Find granularity where most metrics can run
"""

import os

GRANULARITY_ORDER = ["10d", "1m", "3m", "1y"]


def get_available_granularities(variable_file_map):
    """Get all available granularities from the file map."""
    all_grans = set()
    for name, entries in variable_file_map.items():
        for entry in entries:
            # breakpoint()
            if name == "mesh_mask":
                continue
            all_grans.add(entry["granularity"])
    return sorted(all_grans, key=GRANULARITY_ORDER.index)


def analyse_variable_availability(variable_file_map, required_vars=None):
    """
    Core function: analyze what variables are available at what granularities.

    Handles both direct files and resampling possibilities

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        required_vars: List of variables to analyze (if None, analyzes all)

    Returns
    -------
        Dictionary containing:
        - variable_availability: Dict of var -> list of achievable granularities
        - available_granularities: List of all granularities in the system
        - direct_files: Dict of var -> list of granularities with direct files
    """
    available_grans = get_available_granularities(variable_file_map)

    # If no specific variables requested, analyze all
    if required_vars is None:
        all_vars = set()
        for _ in variable_file_map.values():
            all_vars.update(variable_file_map.keys())
        required_vars = all_vars

    variable_availability = {}
    direct_files = {}

    for var in required_vars:
        variable_availability[var] = []

        # Check direct file availability
        direct_available = []
        for entry in variable_file_map.get(var, []):
            if os.path.exists(entry["file"]):
                direct_available.append(entry["granularity"])

        direct_files[var] = direct_available
        print(f"{var}: direct files at {direct_available}")

        # For each target granularity, check if achievable
        for target_gran in available_grans:
            # Direct file available?
            if target_gran in direct_available:
                variable_availability[var].append(target_gran)
                print(f"{var}: achievable at {target_gran} (direct)")
                continue

            # Can we resample from a finer granularity?
            target_rank = GRANULARITY_ORDER.index(target_gran)
            can_resample = any(
                GRANULARITY_ORDER.index(direct_gran) < target_rank
                for direct_gran in direct_available
            )

            if can_resample:
                variable_availability[var].append(target_gran)
                print(f"{var}: achievable at {target_gran} (via resampling)")
            else:
                print(f"{var}: NOT achievable at {target_gran}")

    return {
        "variable_availability": variable_availability,
        "available_granularities": available_grans,
        "direct_files": direct_files,
    }


def analyze_metric_requirements(variable_file_map, metric_requirements):
    """
    Analyze what metrics can run at what granularities.

    Built on top of analyze_variable_availability

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        metric_requirements: Dictionary mapping metric names to required variables

    Returns
    -------
        Dictionary containing all variable availability info plus:
        - runnable_metrics: Dict of granularity -> list of runnable metric names
    """
    # Get all variables needed by metrics
    all_vars = set()
    for vars_list in metric_requirements.values():
        all_vars.update(vars_list)

    # Analyze variable availability
    analysis = analyse_variable_availability(variable_file_map, all_vars)

    # Determine runnable metrics per granularity
    runnable_metrics = {}
    for gran in analysis["available_granularities"]:
        runnable_metrics[gran] = []

        for metric_name, required_vars in metric_requirements.items():
            if all(
                gran in analysis["variable_availability"].get(var, [])
                for var in required_vars
            ):
                runnable_metrics[gran].append(metric_name)

    # Add metrics analysis to the result
    analysis["runnable_metrics"] = runnable_metrics

    # Build a mapping: granularity -> list of variables available at that granularity
    granularity_to_variables = {
        gran: [] for gran in analysis["available_granularities"]
    }
    for var, grans in analysis["variable_availability"].items():
        for gran in grans:
            granularity_to_variables[gran].append(var)

    analysis["granularity_to_variables"] = granularity_to_variables

    # invert the direct files so we index by granularity and return the
    # available variables
    direct_from_file_by_granularity = {
        gran: [] for gran in analysis["available_granularities"]
    }
    for var, direct_grans in analysis["direct_files"].items():
        for gran in direct_grans:
            direct_from_file_by_granularity[gran].append(var)

    analysis["direct_from_file_by_granularity"] = direct_from_file_by_granularity

    return analysis


def get_runnable_metrics_at_granularity(
    variable_file_map, metric_requirements, target_granularity
):
    """
    Get what metrics can run at a specific granularity.

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        metric_requirements: Dictionary mapping metric names to required variables
        target_granularity: Granularity to check (e.g., "1m", "3m")

    Returns
    -------
        List of metric names that can run at the target granularity
    """
    analysis = analyze_metric_requirements(variable_file_map, metric_requirements)
    return analysis["runnable_metrics"].get(target_granularity, [])


# Find the maximum granularity where each metric can be computed
def get_optimal_granularities_for_metrics(variable_file_map, metric_requirements):
    """
    Get the finest granularity where each metric can run.

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        metric_requirements: Dictionary mapping metric names to required variables

    Returns
    -------
        Dictionary mapping metric names to their finest achievable granularity
    """
    analysis = analyze_metric_requirements(variable_file_map, metric_requirements)

    optimal_grans = {}
    for metric_name in metric_requirements:
        # Find finest granularity where this metric can run
        for gran in GRANULARITY_ORDER:  # Finest to coarsest
            if metric_name in analysis["runnable_metrics"].get(gran, []):
                optimal_grans[metric_name] = gran
                break

    return optimal_grans


# Find the maximum granularity where all metrics can be computed
def get_maximum_granularity_with_all(
    variable_file_map, metric_requirements, GRANULARITY_ORDER=None
):
    """
    Find the granularity where the most metrics can run.

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        metric_requirements: Dictionary mapping metric names to required variables
        GRANULARITY_ORDER: Optional custom granularity order

    Returns
    -------
        String representing the best granularity for running the most metrics
    """
    if GRANULARITY_ORDER is None:
        GRANULARITY_ORDER = ["10d", "1m", "3m", "1y"]

    analysis = analyze_metric_requirements(variable_file_map, metric_requirements)

    # Find granularity with most runnable metrics
    best_gran = None
    max_metrics = 0

    for gran, metrics in analysis["runnable_metrics"].items():
        if len(metrics) > max_metrics:
            max_metrics = len(metrics)
            best_gran = gran

    print(
        (
            f"\nBest granularity: {best_gran} "
            f"(can run {max_metrics} metrics: "
            f"{analysis['runnable_metrics'][best_gran]})"
        )
    )
    return best_gran


def show_availability_summary(variable_file_map, metric_requirements):
    """
    Comprehensive summary of what's possible.

    Args:
        variable_file_map: Dictionary mapping variables to file entries
        metric_requirements: Dictionary mapping metric names to required variables

    Returns
    -------
        Full analysis dictionary for further processing
    """
    analysis = analyze_metric_requirements(variable_file_map, metric_requirements)

    print("\n=== AVAILABILITY SUMMARY ===")

    # Variable availability
    print("\nVariable availability:")
    for var, grans in analysis["variable_availability"].items():
        direct = analysis["direct_files"][var]
        resampled = [g for g in grans if g not in direct]
        print(f"  {var}:")
        print(f"    Direct: {direct}")
        if resampled:
            print(f"    Resampled: {resampled}")

    # Metric possibilities
    print("\nMetric possibilities:")
    for gran, metrics in analysis["runnable_metrics"].items():
        if metrics:
            print(f"  At {gran}: {metrics}")

    # Optimal granularities
    optimal = get_optimal_granularities_for_metrics(
        variable_file_map, metric_requirements
    )
    print("\nOptimal granularities (finest possible):")
    for metric, gran in optimal.items():
        print(f"  {metric}: {gran}")

    return analysis
