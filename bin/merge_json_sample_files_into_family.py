#!/usr/bin/env python3

# Released under the MIT license.

import argparse
import json
import sys
from pathlib import Path


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge JSON files per family",
        epilog="Example: python vcfparser.py --file_in vep.vcf --file_out vep.most_severe_csq.vcf --variant_csq variant_consequence.txt",
    )
    parser.add_argument(
        "--files_in",
        metavar="FILES_IN",
        type=list(Path),
        help="JSON files per sample.",
    )
    parser.add_argument(
        "--file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Merged JSON file.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.file_out, "w") as out_json:
        for file in args.files_in:
            if not file.is_file():
                print(f"The given input file {file} was not found!")
                sys.exit(2)
                with open(file, "rt") as in_json:
                read_json(in_vcf, out_vcf, var_csq)


if __name__ == "__main__":
    sys.exit(main())
