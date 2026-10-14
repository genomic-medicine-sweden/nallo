#!/usr/bin/env python3
"""
Collapse Mitorsaw's per-haplotype genotype (GT) field into a standard diploid genotype.

Mitorsaw reports one GT entry per detected mitochondrial haplotype (pseudo-ploidy), e.g.
`GT=1|1|1|1|1` (5 haplotypes, all carrying the ALT allele) or `GT=0|0|0|1|0` (1 of 5
haplotypes carries the ALT allele). This is not a standard diploid genotype and breaks
interoperability with the rest of the pipeline, which expects two-allele GTs (0/0, 0/1, 1/1).

This script rewrites GT based on the presence of the ALT allele across all reported
haplotypes:
  * all haplotypes REF (all 0)    -> 0/0
  * all haplotypes ALT (all 1)    -> 1/1
  * a mix of REF and ALT          -> 0/1 (heteroplasmic)
All other FORMAT fields (DP, AD, VAF, ...) are left untouched.
"""

import argparse
import gzip
import re
import sys

__version__ = "1.0.0"


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def collapse_gt(gt):
    alleles = re.split(r"[|/]", gt)
    known = [allele for allele in alleles if allele != "."]
    if not known:
        return "./."
    if all(allele == "0" for allele in known):
        return "0/0"
    if all(allele != "0" for allele in known):
        return "1/1"
    return "0/1"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", help="Input VCF file (plain or bgzipped)")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = parser.parse_args()

    with open_vcf(args.vcf) as vcf_in:
        for line in vcf_in:
            line = line.rstrip("\n")

            if line.startswith("#"):
                sys.stdout.write(line + "\n")
                continue

            fields = line.split("\t")
            if len(fields) < 10:
                sys.stdout.write(line + "\n")
                continue

            format_keys = fields[8].split(":")
            gt_index = format_keys.index("GT") if "GT" in format_keys else None

            if gt_index is None:
                sys.stdout.write(line + "\n")
                continue

            for sample_index in range(9, len(fields)):
                sample_values = fields[sample_index].split(":")
                sample_values[gt_index] = collapse_gt(sample_values[gt_index])
                fields[sample_index] = ":".join(sample_values)

            sys.stdout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
