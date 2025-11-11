#!/usr/bin/env python3
"""Visualize the radiation fields referenced by a Cloudy grid input file.

The script understands a subset of continuum-defining commands that appear in
grid inputs (`table ism`, `table hm12`, and `cmb`). It reconstructs their
spectra directly from Cloudy's built-in datasets, combines them onto a common
energy grid, and produces a log-log plot of the individual components plus
their sum. No Cloudy run is required.

Limitations:
- Only the commands listed above are supported. The script raises an error if
  other continuum shapes are encountered.
- Extinction commands (e.g. `extinguish column = ...`) are detected and noted in
  the console output, but their effect is not applied to the plotted spectra.
"""

from __future__ import annotations

import argparse
import bisect
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# Physical constants (cgs)
PLANCK_ERG_S = 6.62607015e-27
BOLTZMANN = 1.380649e-16
LIGHT_SPEED = 2.99792458e10
EV_TO_ERG = 1.602176634e-12
RYD_TO_EV = 13.605693122994
RYD_TO_ERG = RYD_TO_EV * EV_TO_ERG
FR1RYD = RYD_TO_ERG / PLANCK_ERG_S
CMB_TEMP = 2.725  # K

# Black's ISM continuum data (see source/parse_table.cpp)
_TABLE_ISM_LOG_NU = [
    6.00,
    10.72,
    11.00,
    11.23,
    11.47,
    11.55,
    11.85,
    12.26,
    12.54,
    12.71,
    13.10,
    13.64,
    14.14,
    14.38,
    14.63,
    14.93,
    15.08,
    15.36,
    15.43,
    16.25,
    17.09,
    18.00,
    23.00,
]
_TABLE_ISM_LOG_NUFNU = [
    -16.708,
    -2.96,
    -2.47,
    -2.09,
    -2.11,
    -2.34,
    -3.66,
    -2.72,
    -2.45,
    -2.57,
    -3.85,
    -3.34,
    -2.30,
    -1.79,
    -1.79,
    -2.34,
    -2.72,
    -2.55,
    -2.62,
    -5.68,
    -6.45,
    -6.30,
    -11.3,
]
_NUM_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?")


def parse_numbers(text: str) -> List[float]:
    values: List[float] = []
    for match in _NUM_PATTERN.findall(text):
        token = match.replace("d", "e").replace("D", "E")
        values.append(float(token))
    return values


@dataclass
class ComponentSpec:
    kind: str
    label: str
    scale: float = 0.0
    redshift: float = 0.0


def logspace(min_energy: float, max_energy: float, n: int) -> List[float]:
    log_min = math.log10(min_energy)
    log_max = math.log10(max_energy)
    step = (log_max - log_min) / (n - 1)
    return [10 ** (log_min + i * step) for i in range(n)]


def parse_first_number(text: str) -> float | None:
    numbers = parse_numbers(text)
    return numbers[0] if numbers else None


def parse_table_scale(line: str) -> float:
    numbers = parse_numbers(line)
    number = numbers[0] if numbers else None
    if number is None:
        return 0.0
    lg_keyword = " log" in line.lower()
    if number > 0.0 and not lg_keyword:
        return math.log10(number)
    return number


def parse_cmb_redshift(line: str) -> float:
    value = parse_first_number(line)
    return max(value, 0.0) if value is not None else 0.0


def parse_hm12_params(line: str) -> Tuple[float, float]:
    numbers = parse_numbers(line)
    if not numbers:
        raise ValueError("table HM12 command is missing a redshift value")
    redshift = numbers[0]
    scale = 0.0
    if len(numbers) > 1:
        raw_scale = numbers[1]
        scale = math.log10(raw_scale) if raw_scale > 0.0 else raw_scale
    return redshift, scale


def parse_grid_file(path: Path) -> Tuple[List[ComponentSpec], List[str]]:
    components: List[ComponentSpec] = []
    extinction_cmds: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            lower = line.lower()
            tokens = lower.split()
            if not tokens:
                continue
            keyword = tokens[0]
            if keyword == "table" and len(tokens) >= 2:
                subtype = tokens[1]
                if subtype != "ism":
                    if subtype == "hm12":
                        if "quas" in lower:
                            raise NotImplementedError(
                                "The QUASAR option on TABLE HM12 is not supported."
                            )
                        redshift, scale = parse_hm12_params(line)
                        components.append(
                            ComponentSpec(
                                kind="hm12",
                                label=line,
                                scale=scale,
                                redshift=redshift,
                            )
                        )
                    else:
                        raise NotImplementedError(
                            f"Unsupported table command: '{line}'. Only TABLE ISM, TABLE HM12, and CMB are handled."
                        )
                else:
                    scale = parse_table_scale(line)
                    components.append(ComponentSpec(kind="table_ism", label=line, scale=scale))
            elif keyword == "cmb":
                redshift = parse_cmb_redshift(line)
                components.append(ComponentSpec(kind="cmb", label=line, redshift=redshift))
            elif keyword == "extinguish":
                extinction_cmds.append(line)
    if not components:
        raise ValueError("No supported radiation-field commands were found in the grid file.")
    return components, extinction_cmds


class RadiationSpectrum:
    label: str

    def plot_points(self) -> Tuple[List[float], List[float]]:
        raise NotImplementedError

    def energy_bounds(self) -> Tuple[float, float]:
        raise NotImplementedError

    def evaluate(self, energies: Sequence[float]) -> List[float]:
        raise NotImplementedError


class DiscreteSpectrum(RadiationSpectrum):
    def __init__(self, label: str, energies: List[float], fluxes: List[float]):
        self.label = label
        self._energies = energies
        self._fluxes = fluxes

    def plot_points(self) -> Tuple[List[float], List[float]]:
        return self._energies, self._fluxes

    def energy_bounds(self) -> Tuple[float, float]:
        return self._energies[0], self._energies[-1]

    def evaluate(self, energies: Sequence[float]) -> List[float]:
        results: List[float] = []
        logs = [math.log10(e) for e in self._energies]
        for energy in energies:
            if energy <= self._energies[0] or energy >= self._energies[-1]:
                results.append(0.0)
                continue
            idx = bisect.bisect_left(self._energies, energy)
            if self._energies[idx] == energy:
                results.append(self._fluxes[idx])
                continue
            e1 = self._energies[idx - 1]
            e2 = self._energies[idx]
            y1 = math.log10(self._fluxes[idx - 1])
            y2 = math.log10(self._fluxes[idx])
            fraction = (math.log10(energy) - logs[idx - 1]) / (logs[idx] - logs[idx - 1])
            results.append(10 ** (y1 + (y2 - y1) * fraction))
        return results


def _table_ism_base() -> Tuple[List[float], List[float]]:
    energies: List[float] = []
    fluxes: List[float] = []
    # custom low-frequency point
    low_energy_ryd = (10 ** _TABLE_ISM_LOG_NU[0]) / FR1RYD
    energies.append(low_energy_ryd * RYD_TO_EV)
    fluxes.append(10 ** (-21.21 - 6.0))
    for log_nu, log_nufnu in zip(_TABLE_ISM_LOG_NU[6:], _TABLE_ISM_LOG_NUFNU[6:]):
        nu_hz = 10 ** log_nu
        energy_ryd = nu_hz / FR1RYD
        energies.append(energy_ryd * RYD_TO_EV)
        fluxes.append(10 ** (log_nufnu - log_nu))
    return energies, fluxes


_TABLE_ISM_ENERGY, _TABLE_ISM_FLUX = _table_ism_base()


def _token_reader(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            for token in stripped.split():
                yield token


@dataclass
class HM12Dataset:
    redshifts: List[float]
    energies_ev: List[float]
    flux_grid: List[List[float]]

    @classmethod
    def from_ascii(cls, path: Path) -> "HM12Dataset":
        tokens = _token_reader(path)
        try:
            int(next(tokens))  # version, not used directly
            ndim = int(next(tokens))
            npar = int(next(tokens))
            for _ in range(npar):
                next(tokens)  # parameter labels (e.g., "z")
            nmods = int(next(tokens))
            ngrid = int(next(tokens))
            wavl_kind = next(tokens).lower()
            wavl_convert = float(next(tokens))
            flux_kind = next(tokens).lower()
            flux_convert = float(next(tokens))

            params: List[List[float]] = []
            for _ in range(nmods):
                row = [float(next(tokens)) for _ in range(npar)]
                params.append(row)
            raw_axis = [float(next(tokens)) * wavl_convert for _ in range(ngrid)]
        except StopIteration as exc:
            raise ValueError("Unexpected end of hm12_galaxy.ascii header") from exc

        redshifts = [row[0] for row in params]

        if ndim != 1 or npar < 1:
            raise ValueError("HM12 dataset format is not the expected 1-D grid")
        if wavl_kind != "lambda":
            raise NotImplementedError("Only wavelength-based HM12 tables are supported")
        if flux_kind != "f_nu":
            raise NotImplementedError("Only F_nu HM12 tables are supported")

        # Convert wavelengths (Angstrom) to frequency (Hz) and energy (eV)
        freqs_hz = [LIGHT_SPEED / (lam * 1e-8) for lam in raw_axis]
        reverse = freqs_hz[0] > freqs_hz[-1]
        if reverse:
            freqs_hz = list(reversed(freqs_hz))
        energies_ev = [PLANCK_ERG_S * nu / EV_TO_ERG for nu in freqs_hz]

        flux_grid: List[List[float]] = []
        for _ in range(nmods):
            try:
                row = [float(next(tokens)) * flux_convert for _ in range(ngrid)]
            except StopIteration as exc:
                raise ValueError("Unexpected end of HM12 flux table") from exc
            if reverse:
                row = list(reversed(row))
            flux_grid.append(row)

        # Ensure redshifts and fluxes share ascending order
        order = sorted(range(len(redshifts)), key=lambda idx: redshifts[idx])
        redshifts_sorted = [redshifts[idx] for idx in order]
        flux_sorted = [flux_grid[idx] for idx in order]

        return cls(redshifts=redshifts_sorted, energies_ev=energies_ev, flux_grid=flux_sorted)

    def interpolate(self, redshift: float) -> List[float]:
        if redshift <= self.redshifts[0]:
            base = self.flux_grid[0]
            return list(base)
        if redshift >= self.redshifts[-1]:
            base = self.flux_grid[-1]
            return list(base)

        idx = bisect.bisect_left(self.redshifts, redshift)
        z_low = self.redshifts[idx - 1]
        z_high = self.redshifts[idx]
        if z_high == z_low:
            return list(self.flux_grid[idx])
        weight = (redshift - z_low) / (z_high - z_low)
        low_flux = self.flux_grid[idx - 1]
        high_flux = self.flux_grid[idx]
        return [
            low_flux[i] + weight * (high_flux[i] - low_flux[i])
            for i in range(len(low_flux))
        ]


_HM12_DATA_CACHE: Optional[HM12Dataset] = None


def get_hm12_dataset(repo_root: Path) -> HM12Dataset:
    global _HM12_DATA_CACHE
    if _HM12_DATA_CACHE is None:
        data_path = repo_root / "data" / "hm12_galaxy.ascii"
        if not data_path.exists():
            raise FileNotFoundError(f"HM12 data file not found at {data_path}")
        _HM12_DATA_CACHE = HM12Dataset.from_ascii(data_path)
    return _HM12_DATA_CACHE


class TableISMSpectrum(DiscreteSpectrum):
    def __init__(self, label: str, scale: float):
        factor = 10 ** scale
        energies = list(_TABLE_ISM_ENERGY)
        fluxes = [f * factor for f in _TABLE_ISM_FLUX]
        super().__init__(label, energies, fluxes)


class HM12Spectrum(DiscreteSpectrum):
    def __init__(self, label: str, redshift: float, scale: float, dataset: HM12Dataset):
        factor = 10 ** scale
        flux = dataset.interpolate(redshift)
        scaled_flux = [f * factor for f in flux]
        super().__init__(label, dataset.energies_ev, scaled_flux)


def planck_flux(energy_ev: float, temperature: float) -> float:
    energy_erg = energy_ev * EV_TO_ERG
    if energy_erg <= 0.0:
        return 0.0
    freq = energy_erg / PLANCK_ERG_S
    exponent = energy_erg / (BOLTZMANN * temperature)
    if exponent > 700:
        return 0.0
    numerator = 2.0 * PLANCK_ERG_S * freq ** 3 / (LIGHT_SPEED ** 2)
    return numerator / math.expm1(exponent)


class CMBSpectrum(RadiationSpectrum):
    def __init__(self, label: str, redshift: float):
        self.label = label
        self._temperature = CMB_TEMP * (1.0 + redshift)
        self._plot_energy = logspace(1e-6, 1e2, 400)
        self._plot_flux = [planck_flux(e, self._temperature) for e in self._plot_energy]

    def plot_points(self) -> Tuple[List[float], List[float]]:
        return self._plot_energy, self._plot_flux

    def energy_bounds(self) -> Tuple[float, float]:
        return self._plot_energy[0], self._plot_energy[-1]

    def evaluate(self, energies: Sequence[float]) -> List[float]:
        return [planck_flux(e, self._temperature) for e in energies]


def build_spectra(specs: Sequence[ComponentSpec], repo_root: Path) -> List[RadiationSpectrum]:
    spectra: List[RadiationSpectrum] = []
    hm12_data: Optional[HM12Dataset] = None
    for spec in specs:
        if spec.kind == "table_ism":
            spectra.append(TableISMSpectrum(spec.label, spec.scale))
        elif spec.kind == "hm12":
            if hm12_data is None:
                hm12_data = get_hm12_dataset(repo_root)
            spectra.append(HM12Spectrum(spec.label, spec.redshift, spec.scale, hm12_data))
        elif spec.kind == "cmb":
            spectra.append(CMBSpectrum(spec.label, spec.redshift))
        else:
            raise NotImplementedError(f"Unsupported component kind: {spec.kind}")
    return spectra


def compute_combined_spectrum(
    spectra: Sequence[RadiationSpectrum], num_points: int = 500
) -> Tuple[List[float], List[float]]:
    min_energy = min(spec.energy_bounds()[0] for spec in spectra)
    max_energy = max(spec.energy_bounds()[1] for spec in spectra)
    grid = logspace(min_energy, max_energy, num_points)
    total = [0.0 for _ in grid]
    for spec in spectra:
        contribution = spec.evaluate(grid)
        total = [a + b for a, b in zip(total, contribution)]
    return grid, total


@dataclass
class SpectrumBundle:
    source: Path
    spectra: List[RadiationSpectrum]
    combined_energy: List[float]
    combined_flux: List[float]


def convert_energy_units(energies: Sequence[float], units: str) -> List[float]:
    if units == "eV":
        return list(energies)
    if units == "ryd":
        return [e / RYD_TO_EV for e in energies]
    if units == "keV":
        return [e / 1e3 for e in energies]
    raise ValueError(f"Unsupported energy units: {units}")


def plot_spectra(
    bundles: Sequence[SpectrumBundle],
    output_path: Path,
    energy_units: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    show_components = len(bundles) == 1
    if show_components:
        # Plot individual components for clarity when only one grid is supplied.
        for spec in bundles[0].spectra:
            energy, flux = spec.plot_points()
            plt.loglog(
                convert_energy_units(energy, energy_units),
                flux,
                label=spec.label,
                linewidth=1.0,
            )

    for bundle in bundles:
        plt.loglog(
            convert_energy_units(bundle.combined_energy, energy_units),
            bundle.combined_flux,
            linewidth=2.0,
            label=f"{bundle.source.stem} total",
        )
    unit_label = {
        "eV": "eV",
        "ryd": "Ryd",
        "keV": "keV",
    }[energy_units]
    plt.xlabel(f"Photon energy [{unit_label}]")
    plt.ylabel("Incident flux (arbitrary units)")
    if energy_units == "eV":
        xmin, xmax = 1e-3, 1e5
    elif energy_units == "ryd":
        xmin = 1e-3 / RYD_TO_EV
        xmax = 1e5 / RYD_TO_EV
    elif energy_units == "keV":
        xmin = 1e-6
        xmax = 1e2
    else:  # pragma: no cover (defensive, choices already enforced)
        xmin, xmax = None, None
    if xmin is not None and xmax is not None:
        plt.xlim(xmin, xmax)
    plt.ylim(1e-27, 1e-10)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "grid_files",
        type=Path,
        nargs="+",
        help="One or more Cloudy grid input files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path for the output figure (PNG)",
    )
    parser.add_argument(
        "--energy-units",
        default="eV",
        choices=("eV", "ryd", "keV"),
        help="Units to use on the energy axis",
    )
    parser.add_argument(
        "--skip-plot",
        action="store_true",
        help="Compute spectra but skip figure generation",
    )
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    grid_paths = [path.resolve() for path in args.grid_files]
    for grid_path in grid_paths:
        if not grid_path.exists():
            raise FileNotFoundError(f"Grid input file not found: {grid_path}")

    if args.output is None:
        if len(grid_paths) == 1:
            output_path = grid_paths[0].with_name(grid_paths[0].stem + "_radiation.png")
        else:
            output_path = Path.cwd() / "radiation_fields_combined.png"
    else:
        output_path = args.output.resolve()

    bundles: List[SpectrumBundle] = []
    extinction_notes: List[Tuple[Path, List[str]]] = []
    for grid_path in grid_paths:
        specs, extinction_cmds = parse_grid_file(grid_path)
        if extinction_cmds:
            extinction_notes.append((grid_path, extinction_cmds))
        spectra = build_spectra(specs, repo_root)
        combined_energy, combined_flux = compute_combined_spectrum(spectra)
        bundles.append(
            SpectrumBundle(
                source=grid_path,
                spectra=spectra,
                combined_energy=combined_energy,
                combined_flux=combined_flux,
            )
        )

    if args.skip_plot:
        print("Computed radiation fields; skipping plot as requested.")
    else:
        plot_spectra(
            bundles=bundles,
            output_path=output_path,
            energy_units=args.energy_units,
        )
        print(f"Saved radiation field plot to {output_path}")

    for path, commands in extinction_notes:
        print(f"Note: extinction commands in {path} are not applied in the plot:")
        for cmd in commands:
            print(f"  - {cmd}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
