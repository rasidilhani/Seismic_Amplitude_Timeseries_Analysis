"""
Converts seismic .p files to CSV format for all years (2011-2022).

Folder structure:
  C:\\Users\\UserA1\\Documents\\GitHub\\Seismic_Amplitude_Timeseries_Analysis\data\2011\
      2007.001\WIZ.NZ\
          2007_001_WIZ_displacement_average.p
          2007_001_WIZ_RSAM_average.p
      2007.002\WIZ.NZ\
          ...

Output:
  C:\\Users\\UserA1\\Documents\\GitHub\\Seismic_Amplitude_Timeseries_Analysis\csv\
      WIZ_NZ_2007.csv
      WIZ_NZ_2008.csv
      ...

Requirements:
  pip install numpy pandas
"""

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import timezone

# ============================================================
#  SETTINGS — edit these if your paths differ
# ============================================================
DATA_ROOT   = Path(r"E:\data\seismic_amplitude_timeseries_out")
OUTPUT_DIR  = Path(r"C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis\csv")
STATION     = "WIZ.NZ"
YEARS       = list(range(2007, 2011))  # 2007 to 2010
TIMEZONE    = "Pacific/Auckland"   # change to "UTC" if you prefer UTC times
# ============================================================


def load_p_file(path: Path): 
    """Load a .p pickle file and return as numpy array, or None on error."""
    try:
        with open(path, "rb") as f:
            return np.array(pickle.load(f))
    except Exception as e:
        print(f"    WARNING: Could not read {path.name} — {e}")
        return None


def arr_to_df(arr, value_col: str) -> pd.DataFrame:
    """Convert [timestamp, value] array to a labelled DataFrame."""
    df = pd.DataFrame(arr, columns=["unix_timestamp", value_col])
    df["datetime_nz"] = (
        pd.to_datetime(df["unix_timestamp"], unit="s", utc=True)
        .dt.tz_convert(TIMEZONE)
        .dt.strftime("%Y-%m-%d %H:%M:%S")
    )
    return df[["unix_timestamp", "datetime_nz", value_col]]



def main():
    # ── Validate base directory ──────────────────────────────────────────────
    if not DATA_ROOT.exists():
        print(f"\nERROR: Data folder not found:\n  {DATA_ROOT}")
        print("Please check the DATA_ROOT path at the top of this script.")
        input("\nPress Enter to exit...")
        return

    # ── Create output directory ──────────────────────────────────────────────
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nSource data: {DATA_ROOT}")
    print(f"Output folder: {OUTPUT_DIR}\n")
    print("=" * 60)

    # ── Process each year ────────────────────────────────────────────────────
    for year in YEARS:
        year_dir = DATA_ROOT / str(year)

        # Check if year folder exists
        if not year_dir.exists():
            print(f"\n⊘  Year {year}: Folder not found at {year_dir}")
            continue

        # ── Find all YYYY.NNN day folders ────────────────────────────────────
        day_folders = sorted(
            p for p in year_dir.iterdir()
            if p.is_dir() and p.name.startswith(f"{year}.")
        )

        if not day_folders:
            print(f"\n⊘  Year {year}: No day folders found")
            continue

        print(f"\n📅 Processing {year}...")
        print(f"   Found {len(day_folders)} day folder(s)\n")

        all_data = []
        ok_count = 0
        skip_count = 0

        for day_folder in day_folders:
            station_path = day_folder / STATION

            # ── Check station folder exists ──────────────────────────────────
            if not station_path.exists():
                skip_count += 1
                continue

            # ── Find displacement and RSAM files ────────────────────────────
            disp_files = list(station_path.glob("*displacement_average.p"))
            rsam_files = list(station_path.glob("*RSAM_average.p"))

            if not disp_files and not rsam_files:
                skip_count += 1
                continue

            # ── Load arrays ──────────────────────────────────────────────────
            disp_arr = load_p_file(disp_files[0]) if disp_files else None
            rsam_arr = load_p_file(rsam_files[0]) if rsam_files else None

            # ── Convert to DataFrames ────────────────────────────────────────
            df_disp = arr_to_df(disp_arr, "displacement_avg_m") if disp_arr is not None else None
            df_rsam = arr_to_df(rsam_arr, "rsam_avg")           if rsam_arr is not None else None

            # ── Merge on timestamp ───────────────────────────────────────────
            if df_disp is not None and df_rsam is not None:
                df = pd.merge(
                    df_disp,
                    df_rsam[["unix_timestamp", "rsam_avg"]],
                    on="unix_timestamp", how="outer"
                ).sort_values("unix_timestamp").reset_index(drop=True)

                # Fill any missing datetime after outer join
                mask = df["datetime_nz"].isna()
                if mask.any():
                    df.loc[mask, "datetime_nz"] = (
                        pd.to_datetime(df.loc[mask, "unix_timestamp"], unit="s", utc=True)
                        .dt.tz_convert(TIMEZONE)
                        .dt.strftime("%Y-%m-%d %H:%M:%S")
                    )

            elif df_disp is not None:
                df = df_disp.copy()
                df["rsam_avg"] = None
            else:
                df = df_rsam.copy()
                df["displacement_avg_m"] = None

            # Reorder columns
            df = df[["unix_timestamp", "datetime_nz", "displacement_avg_m", "rsam_avg"]]
            df = df.where(pd.notna(df), other=None)
            df["day"] = day_folder.name  # Add day identifier

            # ── Accumulate data ──────────────────────────────────────────────
            all_data.append(df)
            ok_count += 1

        # ── Combine all data and save as CSV ─────────────────────────────────
        if ok_count == 0:
            print(f"   No data was processed for {year}.")
            continue

        combined_df = pd.concat(all_data, ignore_index=True)
        combined_df = combined_df[["day", "unix_timestamp", "datetime_nz", "displacement_avg_m", "rsam_avg"]]

        output_file = OUTPUT_DIR / f"WIZ_NZ_{year}.csv"
        combined_df.to_csv(output_file, index=False)

        print(f"   ✅ {year}: {ok_count} day(s) written, {skip_count} skipped")
        print(f"   📄 Saved to: {output_file.name}")
        print(f"      Total rows: {len(combined_df):,}")

    print("\n" + "=" * 60)
    print("✅ All years processed!")
    input("\nPress Enter to exit...")


if __name__ == "__main__":
    main()