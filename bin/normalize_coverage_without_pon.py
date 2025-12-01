#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Tuple

VERSION = "0.1.0"

description = """
Normalize coverage to log2 ratios compared to its own median coverage.

Expected input format is a tab delimited bed file with four columns:

chr start end coverage
"""

COVERAGE_COL_INDEX = 3


def main() -> None:
    args = parse_args()

    headers, rows = read_coverage_table(args.input)
    cov_vals = [row[COVERAGE_COL_INDEX] for row in rows]
    print("Total number of bins:", len(cov_vals))
    median = calculate_median(cov_vals)
    print("Median coverage:", median)
    if median <= 0:
        raise ValueError(f"Mean coverage must be positive, got {median}")

    write_normalized_output(headers, rows, median, args.output)


def calculate_median(coverages: Iterable[float]) -> float:
    sorted_cov = sorted(coverages)
    n = len(sorted_cov)
    mid = n // 2
    if n % 2:
        return sorted_cov[mid]
    return (sorted_cov[mid - 1] + sorted_cov[mid]) / 2


def read_coverage_table(
    path: Path,
) -> Tuple[List[str], List[Tuple[str, int, int, float]]]:
    """Read a coverage table, returning headers and coverage rows."""
    headers: List[str] = []
    rows: List[Tuple[str, int, int, float]] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("@") or line.startswith("CONTIG"):
                headers.append(line)
                continue

            contig, start_str, end_str, coverage_str = line.split("\t")
            try:
                coverage = float(coverage_str)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid coverage value '{coverage_str}' in line: {line}"
                ) from exc

            rows.append((contig, int(start_str), int(end_str), coverage))

    if not rows:
        raise ValueError(f"No coverage rows found in {path}")

    return headers, rows


def write_normalized_output(
    headers: List[str],
    rows: List[Tuple[str, int, int, float]],
    median_cov: float,
    output: Path,
) -> None:
    with output.open("w", encoding="utf-8") as handle:
        for header_line in headers:
            handle.write(f"{header_line}\n")

        for contig, start, end, coverage in rows:
            # Skip those with 0 coverage
            if coverage == 0:
                continue
            else:
                log2_ratio = math.log2(coverage / median_cov)
            handle.write(f"{contig}\t{start}\t{end}\t{log2_ratio}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Four column coverage file (format: chr, start, end, value)",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to write normalized coverage TSV",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


if __name__ == "__main__":
    main()
