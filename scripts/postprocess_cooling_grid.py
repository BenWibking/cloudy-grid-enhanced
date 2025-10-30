#!/usr/bin/env python3
"""Post-process Cloudy grid outputs into CIAOLoop-style cooling tables."""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


class GridPoint:
    def __init__(self, index: int, params: List[float], failed: bool, warned: bool):
        self.index = index
        self.params = params
        self.failed = failed
        self.warned = warned


class CoolingGridProcessor:
    def __init__(self, prefix: str, output_dir: Path):
        self.prefix = prefix
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.grid_path = Path(f"{prefix}_grid.grd")
        self.summary_path = Path(f"{prefix}_summary.txt")
        self.summary_pattern = f"grid?????????_{prefix}_summary.txt"
        self.output_pattern = f"grid?????????_{prefix}_*.out"

        if not self.grid_path.exists():
            raise FileNotFoundError(f"Grid file '{self.grid_path}' not found")

    def process(self) -> None:
        param_names, grid_points = self._parse_grid_file(self.grid_path)
        if not grid_points:
            raise RuntimeError("No grid entries found in save grid output")

        hden_idx, temp_idx = self._identify_parameter_columns(param_names)
        mmw_map = self._load_mean_molecular_weights()
        exec_times = self._collect_execution_times()
        log_warnings = self._collect_log_warnings()

        records: Dict[float, List[Tuple[float, float, float, float, float, float]]] = defaultdict(list)
        failed_details: List[Tuple[int, float, float, str]] = []
        warned_details: List[Tuple[int, float, float, str]] = []

        for point in grid_points:
            logn = point.params[hden_idx]
            logt = point.params[temp_idx]

            if point.failed:
                failed_details.append(
                    (point.index, logn, logt, self._log_file_for_index(point.index))
                )
                continue
            if point.warned:
                warned_details.append(
                    (point.index, logn, logt, self._log_file_for_index(point.index))
                )

            try:
                heating, cooling = self._read_heating_point(point.index)
            except FileNotFoundError:
                failed_details.append(
                    (point.index, logn, logt, self._log_file_for_index(point.index))
                )
                continue
            mmw = mmw_map.get(point.index)
            if mmw is None:
                raise RuntimeError(
                    f"Mean molecular weight not found for grid index {point.index}"
                )

            n_h = 10.0 ** logn
            scale_factor = n_h * n_h
            heating_scaled = heating / scale_factor
            cooling_scaled = cooling / scale_factor
            te_linear = 10.0 ** logt

            records[round(logn, 6)].append(
                (logn, logt, te_linear, heating_scaled, cooling_scaled, mmw)
            )

        self._write_output(records)
        warn_summary = self._aggregate_warnings(log_warnings, warned_details)

        self._write_statistics(exec_times, warn_summary, warned_details, failed_details)
        self._print_summary(warn_summary, warned_details, failed_details)

    def _parse_grid_file(self, path: Path) -> Tuple[List[str], List[GridPoint]]:
        param_names: List[str] = []
        points: List[GridPoint] = []

        with path.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith("#Index"):
                    parts = line.lstrip("#").split("\t")
                    param_names = [p.strip() for p in parts[6:-1]]
                    continue

                if line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 7:
                    continue

                idx = int(parts[0])
                failed = parts[1] == "T"
                warned = parts[2] == "T"
                param_values = [float(parts[6 + j]) for j in range(len(param_names))]
                points.append(GridPoint(idx, param_values, failed, warned))

        return param_names, points

    def _identify_parameter_columns(self, names: Iterable[str]) -> Tuple[int, int]:
        hden_idx = temp_idx = None
        for i, name in enumerate(names):
            lower = name.lower()
            if hden_idx is None and "hden" in lower:
                hden_idx = i
            if temp_idx is None and ("temper" in lower or "const" in lower):
                temp_idx = i
        if hden_idx is None or temp_idx is None:
            raise RuntimeError(
                f"Could not locate hden/temperature columns in save grid header ({names})"
            )
        return hden_idx, temp_idx

    def _read_heating_point(self, index: int) -> Tuple[float, float]:
        filename = Path(f"grid{index:09d}_{self.prefix}_heating.txt")
        if not filename.exists():
            filename = Path(f"grid{index:09d}_{self.prefix}_cooling.txt")
            if not filename.exists():
                raise FileNotFoundError(f"Heating/cooling file for index {index} missing")

        with filename.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                columns = line.split()
                if len(columns) < 4:
                    raise RuntimeError(
                        f"Unexpected format in '{filename}': '{line}'"
                    )
                heating = float(columns[2])
                cooling = float(columns[3])
                return heating, cooling

        raise RuntimeError(f"No data rows found in heating file '{filename}'")

    def _load_mean_molecular_weights(self) -> Dict[int, float]:
        mmw_map: Dict[int, float] = {}
        pattern = re.compile(r"MeanMolecularWeight\s+([0-9Ee+\-.]+)")

        for summary_file in Path(".").glob(self.summary_pattern):
            idx = int(summary_file.name[4:13])
            match = pattern.search(summary_file.read_text())
            if match:
                mmw_map[idx] = float(match.group(1))

        if self.summary_path.exists():
            text = self.summary_path.read_text()
            per_point = re.findall(
                r"grid(\d{9})_.*?MeanMolecularWeight\s+([0-9Ee+\-.]+)", text, re.S
            )
            if per_point:
                for idx_str, mmw_value in per_point:
                    mmw_map[int(idx_str)] = float(mmw_value)
            elif not mmw_map:
                matches = pattern.findall(text)
                for idx, mmw_value in enumerate(matches):
                    mmw_map[idx] = float(mmw_value)

        if not mmw_map:
            raise FileNotFoundError(
                "No mean molecular weight data found; expected save special output"
            )
        return mmw_map

    def _candidate_output_files(self) -> List[Path]:
        files = list(Path(".").glob(self.output_pattern))
        if not files:
            files = list(Path(".").glob(f"{self.prefix}*.out"))
        return files

    def _collect_execution_times(self) -> List[float]:
        pattern = re.compile(r"Cloudy ends:.*?ExecTime(?:\(s\))?\s*(?:=)?\s*([0-9.]+)")
        times: List[float] = []
        for out_file in self._candidate_output_files():
            match = pattern.findall(out_file.read_text())
            if match:
                times.extend(float(val) for val in match)
        return times

    def _collect_log_warnings(self) -> Dict[int, List[str]]:
        warnings: Dict[int, List[str]] = defaultdict(list)
        warn_pattern = re.compile(r"^ WARNING.*", re.MULTILINE)
        for out_file in self._candidate_output_files():
            text = out_file.read_text()
            matches = warn_pattern.findall(text)
            if not matches:
                continue
            idx_match = re.match(r"grid(\d{9})_", out_file.name)
            idx = int(idx_match.group(1)) if idx_match else -1
            warnings[idx].extend(line.strip() for line in matches)
        return warnings

    def _write_output(
        self,
        records: Dict[float, List[Tuple[float, float, float, float, float, float]]],
    ) -> None:
        sorted_densities = sorted(records.items(), key=lambda item: item[0])
        timestamp = datetime.now().strftime("%a %b %d %H:%M:%S %Y")

        for run_index, (logn_key, rows) in enumerate(sorted_densities, start=1):
            rows.sort(key=lambda row: row[1])
            file_path = self.output_dir / f"{self.prefix}_run{run_index}.dat"

            with file_path.open("w", encoding="utf-8") as handle:
                handle.write(f"# {timestamp}\n#\n")
                handle.write("# Cooling Map File\n#\n")
                handle.write("# Loop values:\n")
                handle.write(f"# hden {rows[0][0]:.6f}\n#\n")
                handle.write("# Data Columns:\n")
                handle.write("# Te [K]\n")
                handle.write("# Heating [erg s^-1 cm^3]\n")
                handle.write("# Cooling [erg s^-1 cm^3]\n")
                handle.write("# Mean Molecular Weight [amu]\n#\n")
                handle.write("#Te\t\tHeating\t\tCooling\t\tMMW\n")

                for _, _, te, heating, cooling, mmw in rows:
                    handle.write(
                        f"{te:.6e}\t{heating:.7e}\t{cooling:.7e}\t{mmw:.6f}\n"
                    )

    def _write_statistics(
        self,
        times: List[float],
        warn_summary: Dict[str, Dict[str, float]],
        warned: List[Tuple[int, float, float, str]],
        failed: List[Tuple[int, float, float, str]],
    ) -> None:
        stats_path = self.output_dir / f"{self.prefix}_stats.txt"
        with stats_path.open("w", encoding="utf-8") as handle:
            handle.write("# Execution time histogram\n")
            if times:
                edges = [10, 30, 60, 120, 300, 600]
                labels = [
                    "<10s",
                    "10-30s",
                    "30-60s",
                    "1-2min",
                    "2-5min",
                    "5-10min",
                    ">10min",
                ]
                counts = [0] * len(labels)
                for t in times:
                    if t < edges[0]:
                        counts[0] += 1
                    elif t < edges[1]:
                        counts[1] += 1
                    elif t < edges[2]:
                        counts[2] += 1
                    elif t < edges[3]:
                        counts[3] += 1
                    elif t < edges[4]:
                        counts[4] += 1
                    elif t < edges[5]:
                        counts[5] += 1
                    else:
                        counts[6] += 1
                for label, count in zip(labels, counts):
                    handle.write(f"{label}: {count}\n")
                avg = sum(times) / len(times)
                handle.write(f"Average ExecTime: {avg:.2f} s\n")
            else:
                handle.write("No ExecTime data found.\n")

            handle.write("\n# Warnings (from logs)\n")
            if warn_summary:
                for msg, stats in warn_summary.items():
                    logn_range = _format_range(stats["logn_min"], stats["logn_max"])
                    logt_range = _format_range(stats["logt_min"], stats["logt_max"])
                    handle.write(
                        f"{stats['count']} x {msg} (log n_H {logn_range}, log T {logt_range})\n"
                    )
            else:
                handle.write("None\n")

            handle.write("\n# Grid warnings\n")
            if warned:
                for idx, logn, logt, log_path in warned:
                    handle.write(
                        f"grid index {idx}: log n_H={logn:.3f}, log T={logt:.3f}, log={log_path}\n"
                    )
                handle.write(f"total warnings: {len(warned)}\n")
            else:
                handle.write("None\n")

            handle.write("\n# Failed runs\n")
            if failed:
                for idx, logn, logt, log_path in failed:
                    handle.write(
                        f"grid index {idx}: log n_H={logn:.3f}, log T={logt:.3f}, log={log_path}\n"
                    )
                handle.write(f"total failures: {len(failed)}\n")
            else:
                handle.write("None\n")

    def _print_summary(
        self,
        warn_summary: Dict[str, Dict[str, float]],
        warned: List[Tuple[int, float, float, str]],
        failed: List[Tuple[int, float, float, str]],
    ) -> None:
        if failed:
            print(f"FAILURES: {len(failed)} grid points did not complete (showing up to 10):")
            for idx, logn, logt, log_path in failed[:10]:
                print(
                    f" - grid index {idx}: log n_H={logn:.3f}, log T={logt:.3f}, log={log_path}"
                )
            if len(failed) > 10:
                print(f"   ... {len(failed)-10} more")
        else:
            print("No failed grid points.")

        if warned:
            print(
                f"WARNINGS: {len(warned)} grid points completed with warnings (showing up to 10):"
            )
            for idx, logn, logt, log_path in warned[:10]:
                print(
                    f" - grid index {idx}: log n_H={logn:.3f}, log T={logt:.3f}, log={log_path}"
                )
            if len(warned) > 10:
                print(f"   ... {len(warned)-10} more")

            print("\nWarning summary:")
            for msg, stats in warn_summary.items():
                logn_range = _format_range(stats["logn_min"], stats["logn_max"])
                logt_range = _format_range(stats["logt_min"], stats["logt_max"])
                print(
                    f" - {stats['count']} x {msg}"
                    f" (log n_H {logn_range}, log T {logt_range})"
                )
        else:
            print("No grid points reported warnings.")

    def _log_file_for_index(self, index: int) -> str:
        heating = Path(f"grid{index:09d}_{self.prefix}_heating.txt")
        if heating.exists():
            return heating.name
        summary = Path(f"grid{index:09d}_{self.prefix}_summary.txt")
        if summary.exists():
            return summary.name
        return "<log unavailable>"

    def _aggregate_warnings(
        self,
        log_warnings: Dict[int, List[str]],
        warned: List[Tuple[int, float, float, str]],
    ) -> Dict[str, Dict[str, float]]:
        summary: Dict[str, Dict[str, float]] = {}
        if not warned and not log_warnings:
            return summary

        warned_by_index = {idx: (logn, logt) for idx, logn, logt, _ in warned}

        if log_warnings:
            for idx, messages in log_warnings.items():
                for message in messages:
                    stats = summary.setdefault(
                        message,
                        {
                            "count": 0,
                            "logn_min": float("inf"),
                            "logn_max": float("-inf"),
                            "logt_min": float("inf"),
                            "logt_max": float("-inf"),
                        },
                    )
                    stats["count"] += 1
                    if idx in warned_by_index:
                        logn, logt = warned_by_index[idx]
                        stats["logn_min"] = min(stats["logn_min"], logn)
                        stats["logn_max"] = max(stats["logn_max"], logn)
                        stats["logt_min"] = min(stats["logt_min"], logt)
                        stats["logt_max"] = max(stats["logt_max"], logt)

        if not summary and warned_by_index:
            summary["<grid warning flag>"] = {
                "count": len(warned_by_index),
                "logn_min": min(logn for logn, _ in warned_by_index.values()),
                "logn_max": max(logn for logn, _ in warned_by_index.values()),
                "logt_min": min(logt for _, logt in warned_by_index.values()),
                "logt_max": max(logt for _, logt in warned_by_index.values()),
            }

        return summary


def _format_range(min_val: float, max_val: float) -> str:
    if math.isinf(min_val) or math.isinf(max_val):
        return "unknown"
    return f"{min_val:.3f}-{max_val:.3f}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Post-process Cloudy grid outputs into CIAOLoop-style cooling tables "
            "(assumes coolingScaleFactor = 1)."
        )
    )
    parser.add_argument(
        "--prefix",
        default="isrf_ism",
        help="Prefix used in save commands (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default="converted_cooling_tables",
        help="Directory for generated CIAOLoop-style tables",
    )
    args = parser.parse_args()

    processor = CoolingGridProcessor(args.prefix, Path(args.output_dir))
    processor.process()


if __name__ == "__main__":
    main()
def _format_range(min_val: float, max_val: float) -> str:
    if math.isinf(min_val) or math.isinf(max_val):
        return "unknown"
    return f"{min_val:.3f}-{max_val:.3f}"
