#!/usr/bin/env python3
"""Extract the dominant cooling contributors per temperature from Cloudy grid outputs."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple


@dataclass
class CoolingComponent:
    label: str
    wavelength: float | None
    fraction: float


@dataclass
class CoolingRecord:
    path: Path
    temperature: float
    cooling_rate: float
    hydrogen_density: float
    components: List[CoolingComponent]
    line_total: float
    free_free_total: float
    compton_total: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize the strongest cooling components at each temperature from a "
            "Cloudy grid (save cooling ... last separate)."
        )
    )
    parser.add_argument(
        "grid_dir",
        type=Path,
        help="Directory containing grid?????????_prefix_cooling.txt files.",
    )
    parser.add_argument(
        "--pattern",
        default="grid*_cooling.txt",
        help="Glob pattern (relative to grid_dir) used to find cooling files.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Number of cooling components to retain for each temperature (default: 5).",
    )
    parser.add_argument(
        "--statistic",
        choices=("median", "mean", "max"),
        default="median",
        help="Statistic used to aggregate fractions across densities at the same temperature.",
    )
    parser.add_argument(
        "--round-digits",
        type=int,
        default=6,
        help="Decimal places used when grouping temperatures in Kelvin (default: 6).",
    )
    parser.add_argument(
        "--min-fraction",
        type=float,
        default=0.0,
        help="Discard components whose aggregated fraction is <= this threshold.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("top_cooling_components.csv"),
        help="Output CSV path (default: ./top_cooling_components.csv).",
    )
    parser.add_argument(
        "--totals-output",
        type=Path,
        default=Path("cooling_component_totals.csv"),
        help=(
            "Output CSV path containing per-model totals for line, free-free, "
            "and Compton cooling (default: ./cooling_component_totals.csv)."
        ),
    )
    return parser.parse_args()


def find_cooling_files(grid_dir: Path, pattern: str) -> List[Path]:
    files = sorted(p for p in grid_dir.glob(pattern) if p.is_file())
    if not files:
        raise FileNotFoundError(
            f"No files matching '{pattern}' found inside '{grid_dir}'."
        )
    return files


def read_data_line(path: Path) -> str:
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            return stripped
    raise ValueError(f"No data rows found in {path}")


def physical_conditions_path(cooling_path: Path) -> Path:
    candidate = cooling_path.with_name(
        cooling_path.name.replace("_cooling", "_physical_conditions")
    )
    if not candidate.exists():
        raise FileNotFoundError(
            f"Physical conditions file '{candidate.name}' missing for '{cooling_path.name}'"
        )
    return candidate


def read_hydrogen_density(path: Path) -> float:
    last_line: str | None = None
    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            last_line = stripped
    if last_line is None:
        raise ValueError(f"No data rows found in {path}")
    columns = last_line.split()
    if len(columns) < 3:
        raise ValueError(
            f"Unexpected layout in {path}: need >=3 columns, got '{last_line}'"
        )
    try:
        return float(columns[2])
    except ValueError as exc:
        raise ValueError(f"Could not parse hydrogen density in {path}") from exc


def classify_component(label: str) -> str:
    normalized = label.lower()
    if normalized.startswith("ff"):
        return "free_free"
    if normalized.startswith("comp"):
        return "compton"
    return "line"


def parse_cooling_components(path: Path) -> CoolingRecord:
    line = read_data_line(path)
    columns = line.split("\t")
    if len(columns) < 6:
        raise ValueError(f"Unexpected layout in {path}: fewer than 6 columns")

    try:
        temperature = float(columns[1])
        cooling_rate = float(columns[3])
    except ValueError as exc:
        raise ValueError(f"Could not parse primary columns in {path}") from exc

    hydrogen_density = read_hydrogen_density(physical_conditions_path(path))

    tail = columns[5:]
    if len(tail) % 2 == 1:
        tail = tail[:-1]

    components: List[CoolingComponent] = []
    line_total = 0.0
    free_free_total = 0.0
    compton_total = 0.0
    for raw_label, frac_str in zip(tail[0::2], tail[1::2]):
        label = raw_label.strip()
        fraction_text = frac_str.strip()
        if not label or not fraction_text:
            continue
        try:
            fraction = float(fraction_text)
        except ValueError:
            continue
        tokens = label.split()
        if not tokens:
            continue
        wavelength: float | None
        comp_label: str
        try:
            wavelength = float(tokens[-1])
            comp_label = " ".join(tokens[:-1]) if len(tokens) > 1 else tokens[-1]
        except ValueError:
            wavelength = None
            comp_label = label
        components.append(CoolingComponent(comp_label, wavelength, fraction))

        if fraction > 0.0:
            category = classify_component(comp_label)
            if category == "free_free":
                free_free_total += fraction
            elif category == "compton":
                compton_total += fraction
            else:
                line_total += fraction

    return CoolingRecord(
        path=path,
        temperature=temperature,
        cooling_rate=cooling_rate,
        hydrogen_density=hydrogen_density,
        components=components,
        line_total=line_total,
        free_free_total=free_free_total,
        compton_total=compton_total,
    )


def aggregate_components(
    records: Sequence[CoolingRecord],
    *,
    round_digits: int,
    top_k: int,
    statistic: Callable[[Iterable[float]], float],
    min_fraction: float,
) -> List[Dict[str, float | int | str | None]]:
    grouped: Dict[float, Dict[Tuple[str, float | None], List[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    temp_samples: Dict[float, List[float]] = defaultdict(list)

    for record in records:
        temp_key = round(record.temperature, round_digits)
        temp_samples[temp_key].append(record.temperature)
        bucket = grouped[temp_key]
        for component in record.components:
            if component.fraction <= 0.0:
                continue
            bucket[(component.label, component.wavelength)].append(component.fraction)

    rows: List[Dict[str, float | int | str | None]] = []
    for temp_key in sorted(grouped.keys()):
        samples = temp_samples[temp_key]
        mean_temp = statistics.fmean(samples)
        log_temp = math.log10(mean_temp)
        models_at_temp = len(samples)
        candidates: List[Tuple[float, str, float | None, int]] = []
        for (label, wavelength), fractions in grouped[temp_key].items():
            stat_value = statistic(fractions)
            if stat_value <= min_fraction:
                continue
            candidates.append((stat_value, label, wavelength, len(fractions)))
        candidates.sort(key=lambda item: item[0], reverse=True)
        for rank, (stat_value, label, wavelength, count) in enumerate(
            candidates[:top_k], start=1
        ):
            rows.append(
                {
                    "temperature_K": mean_temp,
                    "log10_temperature": log_temp,
                    "models_at_temperature": models_at_temp,
                    "component_rank": rank,
                    "component_label": label,
                    "wavelength": wavelength,
                    "fraction_statistic": stat_value,
                    "component_occurrences": count,
                }
            )
    return rows


def build_total_rows(records: Sequence[CoolingRecord]) -> List[Dict[str, float | str]]:
    rows: List[Dict[str, float | str]] = []
    for record in records:
        if record.temperature <= 0.0:
            raise ValueError(
                f"Non-positive temperature encountered in {record.path}: {record.temperature}"
            )
        line_rate = record.line_total * record.cooling_rate
        free_free_rate = record.free_free_total * record.cooling_rate
        compton_rate = record.compton_total * record.cooling_rate

        rows.append(
            {
                "temperature_K": record.temperature,
                "hydrogen_density": record.hydrogen_density,
                "total_line_cooling_fraction": record.line_total,
                "total_free_free_fraction": record.free_free_total,
                "total_compton_fraction": record.compton_total,
                "total_line_cooling_rate": line_rate,
                "total_free_free_cooling_rate": free_free_rate,
                "total_compton_cooling_rate": compton_rate,
            }
        )
    return rows


def write_csv(
    path: Path,
    rows: Sequence[Dict[str, float | int | str | None]],
    fieldnames: Sequence[str],
) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write to {path}")

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def statistic_from_name(name: str) -> Callable[[Iterable[float]], float]:
    if name == "median":
        return statistics.median
    if name == "mean":
        return statistics.fmean
    if name == "max":
        return max
    raise ValueError(f"Unsupported statistic '{name}'")


def main() -> None:
    args = parse_args()
    files = find_cooling_files(args.grid_dir, args.pattern)
    stat_fn = statistic_from_name(args.statistic)
    records = [parse_cooling_components(path) for path in files]
    rows = aggregate_components(
        records,
        round_digits=args.round_digits,
        top_k=args.top,
        statistic=stat_fn,
        min_fraction=args.min_fraction,
    )
    totals_rows = build_total_rows(records)

    top_components_fieldnames = [
        "temperature_K",
        "log10_temperature",
        "models_at_temperature",
        "component_rank",
        "component_label",
        "wavelength",
        "fraction_statistic",
        "component_occurrences",
    ]
    totals_fieldnames = [
        "temperature_K",
        "hydrogen_density",
        "total_line_cooling_fraction",
        "total_free_free_fraction",
        "total_compton_fraction",
        "total_line_cooling_rate",
        "total_free_free_cooling_rate",
        "total_compton_cooling_rate",
    ]

    write_csv(args.output, rows, top_components_fieldnames)
    write_csv(args.totals_output, totals_rows, totals_fieldnames)


if __name__ == "__main__":
    main()
