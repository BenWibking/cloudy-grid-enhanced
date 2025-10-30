#!/usr/bin/env python3
"""
Convert Cloudy grid outputs into HDF5 cooling tables using a native *.in file.

This script understands the native Cloudy grid syntax directly, so a
CIAOLoop-generated *.run file is no longer required. Helper routines are
adapted from the `cloudy_cooling_tools/cloudy_grids` package.
"""

from __future__ import annotations

import argparse
import copy
import math
import operator
import re
import sys
from functools import reduce
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import h5py
import numpy as np

NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def get_grid_indices(dims: List[int], index: int) -> List[int]:
    """Return the multidimensional indices corresponding to a flat index."""
    if not dims:
        return []
    indices: List[int] = []
    for dim in reversed(dims):
        indices.append(index % dim)
        index //= dim
    return list(reversed(indices))


def loadTemps(mapFile: str, _gridDimension: List[int]) -> List[float]:
    """Read the temperature column in a single ASCII cooling map."""
    temps: List[float] = []
    with open(mapFile, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            temps.append(float(line.split()[0]))
    return temps


def loadMap(
    mapFile: str,
    gridDimension: List[int],
    indices: List[int],
    gridData: List[np.ndarray],
    grid_temp: List[float],
) -> None:
    """Populate the grid data arrays with values from a map file."""
    temp: List[float] = []
    heat: List[float] = []
    cool: List[float] = []
    mmw: List[float] = []

    with open(mapFile, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            temp.append(float(parts[0]))
            heat.append(float(parts[1]))
            cool.append(float(parts[2]))
            mmw.append(float(parts[3]))

    def new_empty_grid(my_dims: List[int]) -> np.ndarray:
        empty = np.empty(shape=my_dims, dtype=float)
        empty.fill(np.nan)
        return empty

    if len(gridData) == 0:
        my_dims = copy.deepcopy(gridDimension)
        my_dims.append(len(grid_temp))
        gridData.append(np.array(grid_temp, dtype=float))
        gridData.append(new_empty_grid(my_dims))
        gridData.append(new_empty_grid(my_dims))
        gridData.append(new_empty_grid(my_dims))

    grid_idx = np.searchsorted(grid_temp, temp)

    heat_arr = np.asarray(heat, dtype=float)
    cool_arr = np.asarray(cool, dtype=float)
    mmw_arr = np.asarray(mmw, dtype=float)

    gridData[1][tuple(indices)][grid_idx] = heat_arr
    gridData[2][tuple(indices)][grid_idx] = cool_arr
    gridData[3][tuple(indices)][grid_idx] = mmw_arr


def _is_comment(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return True
    return stripped.startswith("#") or stripped.lower().startswith("c ")


def _extract_parameter_name(tokens: Sequence[str]) -> str:
    name_tokens: List[str] = []
    for token in tokens:
        if token in {"=", ","}:
            break
        if NUMBER_RE.fullmatch(token):
            break
        name_tokens.append(token)
    if not name_tokens:
        name_tokens = list(tokens)
    return " ".join(name_tokens).strip()


def _parse_grid_values(line: str) -> List[float]:
    matches = NUMBER_RE.findall(line)
    if len(matches) < 3:
        raise ValueError(f"Unable to parse grid specification: '{line}'")
    start, stop, step = (float(matches[i]) for i in range(3))
    if math.isclose(step, 0.0):
        raise ValueError(f"Grid step cannot be zero (line: '{line}')")
    count = int(round((stop - start) / step)) + 1
    values = [start + i * step for i in range(count)]
    # Guard against round-off so the final point matches the requested upper bound.
    values[-1] = stop
    return values


def _parse_in_file(path: Path) -> Tuple[List[str], List[List[float]], List[float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    parameter_names: List[str] = []
    parameter_values: List[List[float]] = []
    temperature_grid: List[float] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if _is_comment(line):
            i += 1
            continue

        lower = line.lower()
        if " vary" in lower:
            tokens = line.split()
            vary_idx = [j for j, token in enumerate(tokens) if token.lower() == "vary"]
            if not vary_idx:
                i += 1
                continue
            command_tokens = tokens[:vary_idx[0]]
            parameter_name = _extract_parameter_name(command_tokens)

            j = i + 1
            grid_line = None
            while j < len(lines):
                candidate = lines[j].strip()
                if _is_comment(candidate):
                    j += 1
                    continue
                if candidate.lower().startswith("grid"):
                    grid_line = candidate
                    break
                if " vary" in candidate.lower():
                    break
                j += 1

            if grid_line is None:
                raise RuntimeError(f"No grid specification found after '{line}' in {path}")

            values = _parse_grid_values(grid_line)
            if parameter_name.lower().startswith("constant temperature"):
                temperature_grid = values
            else:
                parameter_names.append(parameter_name)
                parameter_values.append(values)

            i = j + 1
            continue

        i += 1

    return parameter_names, parameter_values, temperature_grid


def _infer_prefix(lines: Iterable[str], fallback: str) -> str:
    save_patterns = [
        (re.compile(r'save\s+cooling\s+"([^"]+)"', re.IGNORECASE), "_cooling"),
        (re.compile(r'save\s+heating\s+"([^"]+)"', re.IGNORECASE), "_heating"),
        (re.compile(r'save\s+grid\s+"([^"]+)"', re.IGNORECASE), "_grid"),
    ]

    for pattern, suffix in save_patterns:
        match = next((pattern.search(line) for line in lines if pattern.search(line)), None)
        if match:
            filename = Path(match.group(1)).name
            stem = Path(filename).stem
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
            return stem

    return fallback


def _find_map_file(data_dir: Path, prefix: str, index: int) -> Path:
    candidates = [
        data_dir / f"{prefix}_run{index}.dat",
        data_dir / f"{prefix}_run{index}.txt",
        data_dir / f"{prefix}_run{index:04d}.dat",
        data_dir / f"{prefix}_run{index:04d}.txt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate map file for index {index} with prefix '{prefix}'")


def convert_cooling_from_in(
    in_file: Path,
    output_file: Path,
    data_dir: Path | None = None,
    prefix_override: str | None = None,
) -> None:
    parameter_names, parameter_values, temperature_grid = _parse_in_file(in_file)
    lines = in_file.read_text(encoding="utf-8").splitlines()
    prefix = prefix_override or _infer_prefix(lines, in_file.stem)

    base_dir = data_dir or in_file.parent
    base_dir = base_dir.resolve()

    if parameter_values:
        grid_dimension = [len(values) for values in parameter_values]
        total_runs = reduce(operator.mul, grid_dimension, 1)
    else:
        grid_dimension = []
        total_runs = 1

    map_files = [str(_find_map_file(base_dir, prefix, idx + 1)) for idx in range(total_runs)]

    temp_values: List[float] = []
    for map_file in map_files:
        temp_values.extend(loadTemps(map_file, grid_dimension))
    grid_temp = sorted(set(temp_values))

    grid_data: List[np.ndarray] = []
    for idx, map_file in enumerate(map_files):
        indices = get_grid_indices(list(grid_dimension), idx)
        loadMap(map_file, list(grid_dimension), indices, grid_data, grid_temp)

    if temperature_grid:
        grid_temp_array = np.array(grid_temp)
        log_temp = np.log10(grid_temp_array)
        if not (
            np.allclose(grid_temp_array, np.array(temperature_grid), rtol=1e-5, atol=1e-8)
            or np.allclose(log_temp, np.array(temperature_grid), rtol=1e-5, atol=1e-8)
        ):
            print(
                "Warning: temperature samples in ASCII data do not match the grid defined "
                "in the .in file.",
                file=sys.stderr,
            )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_file, "w") as output:
        names = ["Temperature", "Heating", "Cooling", "MMW"]
        for i, data in enumerate(grid_data):
            dataset = output.create_dataset(names[i], data=data, dtype=">f8")
            dataset.attrs["Dimension"] = np.array(data.shape, dtype=">i8")
            dataset.attrs["Rank"] = np.array(len(data.shape), dtype=">i8")

        for idx, values in enumerate(parameter_values):
            array = np.array(values, dtype=float)
            dataset = output.create_dataset(f"Parameter{idx + 1}", data=array, dtype=">f8")
            dataset.attrs["Dimension"] = np.array(array.shape, dtype=">i8")
            dataset.attrs["Name"] = parameter_names[idx]


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Convert a Cloudy *.in grid run into an HDF5 cooling table."
    )
    parser.add_argument("in_file", type=Path, help="Cloudy grid input file (e.g., my_grid.in)")
    parser.add_argument("output", type=Path, help="Destination HDF5 file")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="Directory containing <prefix>_run*.dat files (defaults to the .in file directory).",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Override for the ASCII file prefix (normally inferred from save commands).",
    )
    args = parser.parse_args(argv)

    convert_cooling_from_in(args.in_file, args.output, args.data_dir, args.prefix)


if __name__ == "__main__":
    main()
