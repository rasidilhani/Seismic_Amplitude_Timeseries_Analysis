r"""
Extracts seismic second-level pickle data for WIZ displacement and RSAM files.

Source files:
  *WIZ.displacement.p
  *WIZ.RSAM.p

Output:
  C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis\data\seconds_data

Years processed:
  2018 to 2022
  E:\Seismic_Amplitude_Timeseries_Analysis
"""

import csv
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

# ============================================================
#  SETTINGS
# ============================================================
DATA_ROOT = Path(r"E:\data\seismic_amplitude_timeseries_out")
OUTPUT_DIR = Path(r"E:\Seismic_Amplitude_Timeseries_Analysis")
STATION = "WIZ.NZ"
YEARS = list(range(2019, 2023))
TIMEZONE = "Pacific/Auckland"
CHUNK_SIZE = 250_000
# ============================================================


def load_p_file(path: Path):
    """Load a .p pickle file and return as a numpy array, or None on error."""
    try:
        with open(path, "rb") as f:
            arr = pickle.load(f)
            arr = np.asarray(arr)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 2)
            return arr
    except Exception as e:
        print(f"    WARNING: Could not read {path.name} — {e}")
        return None


def format_datetime_chunk(unix_chunk):
    """Convert a chunk of UNIX timestamps to NZ datetime strings without loading the whole array at once."""
    dt_index = pd.to_datetime(unix_chunk, unit="s", utc=True)
    dt_index = dt_index.tz_convert(TIMEZONE)
    return dt_index.strftime("%Y-%m-%d %H:%M:%S").astype(str)


def write_day_rows(output_handle, day_name: str, disp_arr, rsam_arr):
    """Write one day's displacement and RSAM rows to the CSV in manageable chunks."""
    if disp_arr is None and rsam_arr is None:
        return 0

    if disp_arr is not None and rsam_arr is not None:
        total_rows = min(len(disp_arr), len(rsam_arr))
    elif disp_arr is not None:
        total_rows = len(disp_arr)
    else:
        total_rows = len(rsam_arr)

    writer = csv.writer(output_handle)
    written = 0

    if disp_arr is not None and rsam_arr is not None:
        for start in range(0, total_rows, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, total_rows)
            disp_chunk = disp_arr[start:end]
            rsam_chunk = rsam_arr[start:end]

            unix_chunk = disp_chunk[:, 0].astype(np.int64)
            datetime_chunk = format_datetime_chunk(unix_chunk)
            disp_values = disp_chunk[:, 1]
            rsam_values = rsam_chunk[:, 1]

            for unix_ts, dt_str, disp_val, rsam_val in zip(
                unix_chunk.tolist(),
                datetime_chunk.tolist(),
                disp_values.tolist(),
                rsam_values.tolist(),
            ):
                writer.writerow([day_name, int(unix_ts), dt_str, disp_val, rsam_val])
                written += 1
    elif disp_arr is not None:
        for start in range(0, total_rows, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, total_rows)
            disp_chunk = disp_arr[start:end]
            unix_chunk = disp_chunk[:, 0].astype(np.int64)
            datetime_chunk = format_datetime_chunk(unix_chunk)
            disp_values = disp_chunk[:, 1]

            for unix_ts, dt_str, disp_val in zip(unix_chunk.tolist(), datetime_chunk.tolist(), disp_values.tolist()):
                writer.writerow([day_name, int(unix_ts), dt_str, disp_val, None])
                written += 1
    else:
        for start in range(0, total_rows, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, total_rows)
            rsam_chunk = rsam_arr[start:end]
            unix_chunk = rsam_chunk[:, 0].astype(np.int64)
            datetime_chunk = format_datetime_chunk(unix_chunk)
            rsam_values = rsam_chunk[:, 1]

            for unix_ts, dt_str, rsam_val in zip(unix_chunk.tolist(), datetime_chunk.tolist(), rsam_values.tolist()):
                writer.writerow([day_name, int(unix_ts), dt_str, None, rsam_val])
                written += 1

    return written


def find_target_files(station_path: Path):
    """Return displacement and RSAM files matching the requested naming pattern."""
    disp_files = [
        p for p in station_path.iterdir()
        if p.is_file() and p.name.endswith("WIZ.displacement.p")
    ]
    rsam_files = [
        p for p in station_path.iterdir()
        if p.is_file() and p.name.endswith("WIZ.RSAM.p")
    ]
    return sorted(disp_files), sorted(rsam_files)


def main():
    if not DATA_ROOT.exists():
        print(f"\nERROR: Data folder not found:\n  {DATA_ROOT}")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nSource data: {DATA_ROOT}")
    print(f"Output folder: {OUTPUT_DIR}\n")
    print("=" * 60)

    for year in YEARS:
        year_dir = DATA_ROOT / str(year)
        if not year_dir.exists():
            print(f"\n⊘  Year {year}: Folder not found at {year_dir}")
            continue

        day_folders = sorted(
            p for p in year_dir.iterdir()
            if p.is_dir() and p.name.startswith(f"{year}.")
        )

        if not day_folders:
            print(f"\n⊘  Year {year}: No day folders found")
            continue

        print(f"\n📅 Processing {year}...")
        print(f"   Found {len(day_folders)} day folder(s)\n")

        output_file = OUTPUT_DIR / f"WIZ_NZ_{year}.csv"
        ok_count = 0
        skip_count = 0
        total_rows = 0

        with output_file.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["day", "unix_timestamp", "datetime_nz", "displacement", "rsam"])

            for day_folder in day_folders:
                station_path = day_folder / STATION
                if not station_path.exists():
                    skip_count += 1
                    continue

                disp_files, rsam_files = find_target_files(station_path)
                if not disp_files and not rsam_files:
                    skip_count += 1
                    continue

                disp_arr = load_p_file(disp_files[0]) if disp_files else None
                rsam_arr = load_p_file(rsam_files[0]) if rsam_files else None
                written = write_day_rows(handle, day_folder.name, disp_arr, rsam_arr)

                if written > 0:
                    ok_count += 1
                    total_rows += written
                else:
                    skip_count += 1

        if ok_count == 0:
            print(f"   No data was processed for {year}.")
            continue

        print(f"   ✅ {year}: {ok_count} day(s) written, {skip_count} skipped")
        print(f"   📄 Saved to: {output_file.name}")
        print(f"      Total rows: {total_rows:,}")

    print("\n" + "=" * 60)
    print("✅ All years processed!")


if __name__ == "__main__":
    main()
