#!/usr/bin/env python3
"""
Display progress for an ongoing Cloudy grid run by inspecting per-point
files (e.g. `grid########_isrf_ism_heating.txt`). Only grid points with
completed heating (or cooling) files are counted as finished, so the
progress reflects actual completed simulations.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import time
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

GRID_CMD_RE = re.compile(
    r"grid\s+range\s+from\s+([+-]?\d+(?:\.\d+)?)\s+to\s+([+-]?\d+(?:\.\d+)?)\s+step\s+([+-]?\d+(?:\.\d+)?)",
    re.IGNORECASE,
)
HDEN_RE = re.compile(r"HDEN\s*=\s*([+-]?\d+(?:\.\d+)?)", re.IGNORECASE)
TEMP_RE = re.compile(r"CONSTANT\s+TEMP\s*([+-]?\d+(?:\.\d+)?)", re.IGNORECASE)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Show Cloudy grid progress")
    parser.add_argument(
        "--prefix",
        default="isrf_ism",
        help="Prefix used in per-point filenames (default: %(default)s)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Cloudy input file to read grid definitions from (optional)",
    )
    parser.add_argument(
        "--detail",
        action="store_true",
        help="List unique log n_H and log T values completed",
    )
    return parser.parse_args()


def parse_grid_spec_from_text(text: str) -> List[Tuple[float, float, float]]:
    return [tuple(map(float, match.groups())) for match in GRID_CMD_RE.finditer(text)]


def infer_grid_specs(input_path: Optional[Path]) -> List[Tuple[float, float, float]]:
    if input_path and input_path.exists():
        return parse_grid_spec_from_text(input_path.read_text(encoding="utf-8"))

    per_point_inputs = sorted(Path(".").glob("grid?????????_*.in"))
    if per_point_inputs:
        return parse_grid_spec_from_text(per_point_inputs[0].read_text(encoding="utf-8"))
    return []


def count_points(start: float, end: float, step: float) -> int:
    if step == 0:
        raise ValueError("Grid step cannot be zero")
    n_steps = (end - start) / step
    n_steps = int(round(n_steps))
    return abs(n_steps) + 1


def collect_completed_indices(prefix: str) -> List[Tuple[int, float]]:
    indices = []
    for path in Path(".").glob(f"grid?????????_{prefix}_heating.txt"):
        match = re.match(r"^grid(\d{9})_", path.name)
        if match:
            mtime = path.stat().st_mtime
            indices.append((int(match.group(1)), mtime))
    if not indices:
        for path in Path(".").glob(f"grid?????????_{prefix}_cooling.txt"):
            match = re.match(r"^grid(\d{9})_", path.name)
            if match:
                mtime = path.stat().st_mtime
                indices.append((int(match.group(1)), mtime))
    indices.sort(key=lambda x: x[0])
    return indices


def collect_completed_values(indices: Iterable[int]) -> Tuple[List[float], List[float]]:
    hden_values: List[float] = []
    temp_values: List[float] = []
    for idx in indices:
        for path in Path(".").glob(f"grid{idx:09d}_*.in"):
            text = path.read_text(encoding="utf-8")
            h_match = HDEN_RE.search(text)
            t_match = TEMP_RE.search(text)
            if h_match and t_match:
                hden_values.append(float(h_match.group(1)))
                temp_values.append(float(t_match.group(1)))
                break
    return hden_values, temp_values


def make_progress_bar(completed: int, total: int, width: int = 40) -> str:
    ratio = 0.0 if total == 0 else completed / total
    filled = min(width, int(round(ratio * width)))
    return f"[{'=' * filled}{'.' * (width - filled)}] {ratio:.2%}"


def format_eta(remaining: int, avg_duration: float) -> str:
    if remaining <= 0 or avg_duration <= 0:
        return "00:00:00"
    seconds = int(round(remaining * avg_duration))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def main() -> None:
    args = parse_arguments()
    specs = infer_grid_specs(args.input)
    if len(specs) < 2:
        print("Could not infer both density and temperature grid definitions.", file=sys.stderr)
        sys.exit(1)

    density_spec, temp_spec = specs[0], specs[1]
    density_count = count_points(*density_spec)
    temp_count = count_points(*temp_spec)
    total_points = density_count * temp_count

    completed_data = collect_completed_indices(args.prefix)
    completed_indices = [idx for idx, _ in completed_data]
    completed_points = len(completed_indices)

    print(f"Grid prefix: {args.prefix}")
    print(
        f"Density grid: {density_count} points (log {density_spec[0]} to {density_spec[1]} in steps of {density_spec[2]})"
    )
    print(
        f"Temperature grid: {temp_count} points (log {temp_spec[0]} to {temp_spec[1]} in steps of {temp_spec[2]})"
    )

    print(f"Completed models: {completed_points} / {total_points}")
    print(make_progress_bar(completed_points, total_points))
    remaining = total_points - completed_points
    print(f"Remaining points: {remaining}")

    if completed_points > 1:
        timestamps = [mtime for _, mtime in completed_data]
        timestamps.sort()
        durations = [t2 - t1 for t1, t2 in zip(timestamps[:-1], timestamps[1:]) if t2 > t1]
        avg_duration = sum(durations) / len(durations) if durations else 0.0
        eta_str = format_eta(remaining, avg_duration)
        print(f"Estimated time to completion: {eta_str}")
    else:
        print("Estimated time to completion: --")

    if args.detail and completed_points > 0:
        hden_vals, temp_vals = collect_completed_values(completed_indices)
        if hden_vals:
            unique_hden = sorted(set(round(v, 5) for v in hden_vals))
            print(f"log n_H completed ({len(unique_hden)} values): {unique_hden}")
        if temp_vals:
            unique_temp = sorted(set(round(v, 5) for v in temp_vals))
            print(f"log T completed   ({len(unique_temp)} values): {unique_temp}")


if __name__ == "__main__":
    main()
