#!/usr/bin/env python3
"""
Normalise an SV VCF before VEP annotation.

Remaps non-canonical SVTYPEs (TRA, DUP/INS, DEL/INV, INVDUP) to
VEP-recognised equivalents and records the original in ORIG_SVTYPE.
"""

import argparse
import gzip
import re

__version__ = "1.0.0"

SVTYPE_REMAP = {
    'TRA':        'BND',  # observed: dysgu, debreak
    'DUP/INS':    'DUP',
    'DEL/INV':    'DEL',
    'INVDUP':     'DUP',
    'INV/DEL':    'DEL',  # defensive
    'INV/INVDUP': 'DUP',  # defensive
}

_SVTYPE_RE = re.compile(r'SVTYPE=([^;]+)')


def remap_svtype(info):
    """Remap non-canonical SVTYPE in an INFO field string (tab-delimited VCF column 8)."""
    m = _SVTYPE_RE.search(info)
    if not m:
        return info
    svtype = m.group(1)
    canonical = SVTYPE_REMAP.get(svtype)
    if canonical is None:
        return info
    return info[: m.start()] + f'SVTYPE={canonical}' + info[m.end() :] + f';ORIG_SVTYPE={svtype}'


def main():
    parser = argparse.ArgumentParser(
        description='Normalise SV VCF for VEP: remap non-canonical SVTYPEs'
    )
    parser.add_argument('vcf', help='Input VCF file')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    orig_svtype_header_added = False

    opener = gzip.open(args.vcf, 'rt') if args.vcf.endswith('.gz') else open(args.vcf)
    for line in opener:
        line = line.rstrip('\n')

        if line.startswith('#'):
            if not orig_svtype_header_added and line.startswith('#CHROM'):
                print('##INFO=<ID=ORIG_SVTYPE,Number=1,Type=String,Description="Original SVTYPE before VEP normalisation">')
                orig_svtype_header_added = True
            print(line)
            continue

        fields = line.split('\t')
        if len(fields) >= 8:
            fields[7] = remap_svtype(fields[7])
            print('\t'.join(fields))
        else:
            print(line)


if __name__ == '__main__':
    main()
