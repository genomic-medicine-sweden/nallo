#!/usr/bin/env python3

import argparse
import re

SVTYPE_REMAP = {
    'TRA':        'BND',  # observed: sniffles1, dysgu, debreak
    'DUP/INS':    'DUP',  # observed: sniffles1
    'DEL/INV':    'DEL',  # observed: sniffles1 (SVLEN negative → net deletion)
    'INVDUP':     'DUP',  # observed: sniffles1 (original copy intact, second copy added)
    'INV/DEL':    'DEL',  # defensive
    'INV/INVDUP': 'DUP',  # defensive
}


def remap_svtype(info):
    m = re.search(r'SVTYPE=([^;]+)', info)
    if not m:
        return info
    svtype = m.group(1)
    canonical = SVTYPE_REMAP.get(svtype)
    if canonical is None:
        return info
    info = re.sub(r'SVTYPE=[^;]+', f'SVTYPE={canonical}', info)
    info = info + f';ORIG_SVTYPE={svtype}'
    return info


def main():
    parser = argparse.ArgumentParser(
        description='Normalise SV VCF for VEP: remap non-canonical SVTYPEs and apply caller-specific fixes'
    )
    parser.add_argument('vcf', help='Input VCF file')
    parser.add_argument('--caller', default='', help='Caller identifier (e.g. sniffles1, sniffles2, severus)')
    args = parser.parse_args()

    strandbias_added = False
    orig_svtype_header_added = False

    for line in open(args.vcf):
        line = line.rstrip('\n')

        if line.startswith('#'):
            if not orig_svtype_header_added and line.startswith('#CHROM'):
                print('##INFO=<ID=ORIG_SVTYPE,Number=1,Type=String,Description="Original SVTYPE before VEP normalisation">')
                orig_svtype_header_added = True
            # Sniffles v1 uses STRANDBIAS in FILTER column but omits it from the header
            if args.caller == 'sniffles1' and not strandbias_added and line.startswith('##FILTER'):
                print('##FILTER=<ID=STRANDBIAS,Description="Variant supported by reads from a single strand only">')
                strandbias_added = True
            print(line)
            continue

        # Sniffles v1: strip synthetic END from INS/BND/DUP+INS records
        # (END on TRA is the CHR2 coordinate and is meaningful — do not strip)
        if args.caller == 'sniffles1':
            if 'SVTYPE=INS' in line or 'SVTYPE=BND' in line or 'SVTYPE=DUP/INS' in line:
                parts = line.split(';END=')
                line = parts[0] + parts[-1].lstrip('0123456789')
            # Strip synthetic SVLEN=1 from BND and TRA (interchromosomal events)
            if 'SVTYPE=BND' in line or 'SVTYPE=TRA' in line:
                parts = line.split(';SVLEN=1')
                line = parts[0] + parts[-1]

        # Generic: remap compound SVTYPEs to VEP-canonical types, store original in ORIG_SVTYPE
        fields = line.split('\t')
        if len(fields) >= 8:
            fields[7] = remap_svtype(fields[7])
            print('\t'.join(fields))
        else:
            print(line)


if __name__ == '__main__':
    main()
