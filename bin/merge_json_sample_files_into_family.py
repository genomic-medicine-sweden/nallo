
import argparse
import json
import os


def merge_json_files(json_files, output_file):
    """
    This script merges multiple JSON files into a single JSON file.
    Each input file's content is stored under a key derived from the first part of its filename (before the first dot).

    Usage:
        python merge_json.py --files_in file1.json file2.json file3.json -o output.json

    Arguments:
        --files_in: List of JSON files to merge.
        -o, --output: Name of the output JSON file (default: merged.json).

    Example:
        python merge_json.py --files_in data1.json data2.json -o merged.json
    """
    merged_data = {}

    for file in json_files:
        try:
            with open(file, 'r', encoding='utf-8') as f:
                data = json.load(f)
                file_key = os.path.basename(file).split('.')[0]  # Extract first part of filename
                merged_data[file_key] = data
        except Exception as e:
            print(f"Error reading {file}: {e}")

    try:
        with open(output_file, 'w', encoding='utf-8') as out_f:
            json.dump(merged_data, out_f, indent=4)
        print(f"Merged JSON saved to {output_file}")
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple JSON files into one JSON file with filenames as keys.")
    parser.add_argument("--files_in", nargs='+', required=True, help="List of JSON files to merge")
    parser.add_argument("-o", "--output", default="merged.json", help="Output JSON file name (default: merged.json)")

    args = parser.parse_args()
    merge_json_files(args.files_in, args.output)

