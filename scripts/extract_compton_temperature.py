#!/usr/bin/env python3
"""Extract Compton temperatures from a Cloudy grid .out file."""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List


T_COMP_PATTERN = re.compile(r"T\(Comp\):\s*([+-]?\d+(?:\.\d+)?[EeDd][+-]?\d+)")
ZONE_PATTERN = re.compile(
    r"^####\s+\d+\s+Te:([+-]?\d+(?:\.\d+)?[EeDd][+-]?\d+)\s+Hden:([+-]?\d+(?:\.\d+)?[EeDd][+-]?\d+)"
)


@dataclass
class ComptonRecord:
    index: int
    temperature: float
    hydrogen_density: float
    compton_temperature: float

    @property
    def log_temperature(self) -> float:
        return math.log10(self.temperature)

    @property
    def log_hydrogen_density(self) -> float:
        return math.log10(self.hydrogen_density)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse a Cloudy grid output (.out) file, extract the Compton temperature "
            "reported for each model, and save a CSV keyed by hydrogen density and gas "
            "temperature."
        )
    )
    parser.add_argument(
        "grid_output",
        type=Path,
        help="Path to the Cloudy grid .out file (e.g., hm_2012_cooling_grid_nograins.out).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("compton_temperature_vs_nh.csv"),
        help="Destination CSV file (default: ./compton_temperature_vs_nh.csv).",
    )
    return parser.parse_args()


def _to_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def extract_compton_records(path: Path) -> List[ComptonRecord]:
    records: List[ComptonRecord] = []
    pending_tcomp: float | None = None
    pending_needs_zone = False
    index = 0
    current_hden: float | None = None
    current_temp: float | None = None

    with path.open("r", encoding="ascii", errors="ignore") as handle:
        for line_number, line in enumerate(handle, start=1):
            upper = line.upper()
            if "* HDEN=" in upper:
                match = re.search(r"HDEN\s*=\s*([+-]?\d+(?:\.\d+)?)", upper)
                if match:
                    value = float(match.group(1))
                    if "LOG" in upper:
                        current_hden = 10.0 ** value
                    else:
                        current_hden = value
                continue

            if "* CONSTANT" in upper and "TEMP" in upper and (" LOG" in upper or " LINEAR" in upper):
                match = re.search(r"TEMP(?:ERATURE)?\s+([+-]?\d+(?:\.\d+)?)", upper)
                if match:
                    value = float(match.group(1))
                    if "LOG" in upper:
                        current_temp = 10.0 ** value
                    else:
                        current_temp = value
                continue

            match = T_COMP_PATTERN.search(line)
            if match:
                pending_tcomp = _to_float(match.group(1))
                if current_temp is not None and current_hden is not None:
                    records.append(
                        ComptonRecord(
                            index=index,
                            temperature=current_temp,
                            hydrogen_density=current_hden,
                            compton_temperature=pending_tcomp,
                        )
                    )
                    index += 1
                    pending_tcomp = None
                    pending_needs_zone = False
                else:
                    pending_needs_zone = True
                continue

            if pending_needs_zone and pending_tcomp is not None:
                zone_match = ZONE_PATTERN.search(line)
                if zone_match:
                    temperature = _to_float(zone_match.group(1))
                    hydrogen_density = _to_float(zone_match.group(2))
                    records.append(
                        ComptonRecord(
                            index=index,
                            temperature=temperature,
                            hydrogen_density=hydrogen_density,
                            compton_temperature=pending_tcomp,
                        )
                    )
                    index += 1
                    pending_tcomp = None
                    pending_needs_zone = False

    if pending_needs_zone and pending_tcomp is not None:
        raise RuntimeError(
            "Reached end of file while waiting for a zone line containing Te/Hden "
            "after encountering a T(Comp) entry."
        )
    if not records:
        raise RuntimeError(
            f"No Compton temperatures found in {path}. "
            "Ensure this is a Cloudy grid output with T(Comp) lines."
        )
    return records


def write_csv(path: Path, records: Iterable[ComptonRecord]) -> None:
    rows = list(records)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "model_index",
        "temperature_K",
        "log10_temperature",
        "hydrogen_density",
        "log10_hydrogen_density",
        "compton_temperature_K",
    ]
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in rows:
            writer.writerow(
                {
                    "model_index": record.index,
                    "temperature_K": f"{record.temperature:.6e}",
                    "log10_temperature": f"{record.log_temperature:.6f}",
                    "hydrogen_density": f"{record.hydrogen_density:.6e}",
                    "log10_hydrogen_density": f"{record.log_hydrogen_density:.6f}",
                    "compton_temperature_K": f"{record.compton_temperature:.6e}",
                }
            )


def main() -> None:
    args = parse_args()
    records = extract_compton_records(args.grid_output)
    write_csv(args.output, records)
    print(f"Wrote {len(records)} rows to {args.output}")


if __name__ == "__main__":
    main()
