"""Parameter grids used by the two ArchLink clustering stages."""

FULL_INITIAL_GRID = {
    "resolution_parameters": [1, 5, 10, 30, 60, 90, 110, 130, 150, 200],
    "bandwidths": [0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
    "partgraph_ratios": [60, 100, 80],
    "max_edges": [60, 80, 100],
}

FAST_INITIAL_GRID = {
    "resolution_parameters": [30, 90],
    "bandwidths": [0.1, 0.2],
    "partgraph_ratios": [80],
    "max_edges": [80, 100],
}

FULL_FINAL_GRID = {
    "resolution_parameters": [1, 5, 10, 30, 60, 90, 110, 130, 150, 200],
    "max_edges": [60, 70, 80, 85, 100],
}

FAST_FINAL_GRID = {
    "resolution_parameters": [30, 90],
    "max_edges": [80, 100],
}


def normalize_clustering_mode(mode: str) -> str:
    """Validate and normalize the clustering mode."""
    normalized = str(mode or "full").strip().lower()
    if normalized not in {"full", "fast"}:
        raise ValueError(
            f"Unknown clustering mode '{mode}'. Expected 'full' or 'fast'."
        )
    return normalized


def get_initial_clustering_grid(mode: str) -> dict:
    """Return the parameter grid for the first Leiden scan."""
    normalized = normalize_clustering_mode(mode)
    source = FAST_INITIAL_GRID if normalized == "fast" else FULL_INITIAL_GRID
    return {key: list(values) for key, values in source.items()}


def get_final_clustering_grid(mode: str) -> dict:
    """Return the parameter grid for the post-RF Leiden scan."""
    normalized = normalize_clustering_mode(mode)
    source = FAST_FINAL_GRID if normalized == "fast" else FULL_FINAL_GRID
    return {key: list(values) for key, values in source.items()}
