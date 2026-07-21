from __future__ import annotations

import csv
import pickle
from pathlib import Path

import numpy as np

SOURCE_ROOT = Path(r"E:\data\seismic_amplitude_timeseries_out\2019")
OUTPUT_ROOT = Path(r"E:\data\seismic_amplitude_timeseries_out\2019.WIZ.RSAM")
STATION_NAME = "WIZ.NZ"
TARGET_SUFFIX = ".WIZ.RSAM.p"
START_DAY = 204


def load_pickle_array(path: Path) -> np.ndarray | None:
    try:
        with path.open("rb") as fh:
            obj = pickle.load(fh)
        arr = np.asarray(obj)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 2)
        return arr
    except Exception as exc:
        print(f"WARNING: could not read {path} -> {exc}")
        return None


def write_csv(path: Path, arr: np.ndarray) -> bool:
    if arr is None:
        return False

    if arr.size == 0:
        print(f"EMPTY SOURCE: {path}")
        return False

    if arr.ndim == 1:
        arr = arr.reshape(-1, 2)
    if arr.ndim != 2 or arr.shape[1] < 2:
        print(f"WARNING: unexpected shape in {path}: {arr.shape}")
        return False

    output_path = OUTPUT_ROOT / f"{path.stem}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["unix_timestamp", "rsam"])
        for row in arr:
            unix_ts = int(row[0])
            rsam_value = row[1]
            writer.writerow([unix_ts, rsam_value])

    print(f"CONVERTED: {path.name} -> {output_path.name}")
    return True


def main() -> None:
    if not SOURCE_ROOT.exists():
        print(f"ERROR: source folder not found: {SOURCE_ROOT}")
        return

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Source root: {SOURCE_ROOT}")
    print(f"Output root: {OUTPUT_ROOT}")
    print("=" * 70)

    source_files = []
    for day_dir in sorted(SOURCE_ROOT.iterdir()):
        if not day_dir.is_dir() or not day_dir.name.startswith("2019."):
            continue

        try:
            day_number = int(day_dir.name.split(".")[1])
        except (IndexError, ValueError):
            continue

        if day_number < START_DAY:
            continue

        station_dir = day_dir / STATION_NAME
        if not station_dir.exists() or not station_dir.is_dir():
            continue
        source_files.extend(sorted(station_dir.glob(f"*{TARGET_SUFFIX}")))

    if not source_files:
        print("No matching RSAM pickle files were found.")
        return

    converted = 0
    empty_or_failed = 0
    overwritten_existing = 0

    for source_file in source_files:
        output_path = OUTPUT_ROOT / f"{source_file.stem}.csv"
        if output_path.exists():
            output_path.unlink()
            overwritten_existing += 1
            print(f"REPROCESSING existing output: {output_path.name}")

        arr = load_pickle_array(source_file)
        if arr is None or arr.size == 0:
            empty_or_failed += 1
            continue
        if write_csv(source_file, arr):
            converted += 1

    print("=" * 70)
    print(f"Completed: {converted} file(s) converted, {empty_or_failed} empty/failed file(s) skipped, {overwritten_existing} existing file(s) rebuilt")
    print(f"Output folder: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
