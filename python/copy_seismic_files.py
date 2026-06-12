import os
import shutil
from pathlib import Path

# Define source and destination directories
SOURCE_DIR = r"E:\data\seismic_amplitude_timeseries_out"
DEST_DIR = r"C:\Users\UserA1\Documents\GitHub\Seismic_Amplitude_Timeseries_Analysis\data"

def copy_files_preserving_structure(source, destination):
    """
    Copy files from source to destination while preserving the folder structure.
    
    Args:
        source (str): Source directory path
        destination (str): Destination directory path
    """
    # Convert to Path objects for easier handling
    source_path = Path(source)
    dest_path = Path(destination)
    
    # Verify source exists
    if not source_path.exists():
        print(f"Error: Source directory does not exist: {source}")
        return
    
    # Create destination directory if it doesn't exist
    dest_path.mkdir(parents=True, exist_ok=True)
    
    # Walk through source directory
    file_count = 0
    for root, dirs, files in os.walk(source_path):
        # Calculate relative path from source
        rel_path = Path(root).relative_to(source_path)
        
        # Create corresponding directory in destination
        target_dir = dest_path / rel_path
        target_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy only files ending with _average.p
        for file in files:
            if file.endswith("_average.p"):
                source_file = Path(root) / file
                dest_file = target_dir / file
                
                try:
                    shutil.copy2(source_file, dest_file)
                    print(f"Copied: {source_file} -> {dest_file}")
                    file_count += 1
                except Exception as e:
                    print(f"Error copying {source_file}: {e}")
    
    print(f"\nTotal files copied: {file_count}")

if __name__ == "__main__":
    print(f"Starting file copy operation...")
    print(f"Source: {SOURCE_DIR}")
    print(f"Destination: {DEST_DIR}\n")
    
    copy_files_preserving_structure(SOURCE_DIR, DEST_DIR)
    
    print("\nFile copy operation completed!")
