import argparse
import csv
import datetime
from multiprocessing import Pool
from pathlib import Path
import re
import struct
from zoneinfo import ZoneInfo

import pysam
import vcfpy


def bc_flip(b: str) -> str:
    return b.translate(bc_flip.TRANS)


bc_flip.TRANS = str.maketrans('ATGC', 'TACG')


def read_manifest(*, manifest_file, refseq_file):
    refseq = pysam.FastaFile(str(manifest_file.parent / refseq_file) if not Path(refseq_file).is_absolute() else str(refseq_file))

    def is_snv(row):
        snp = row.get('SNP')
        if not isinstance(snp, str):
            return False
        return re.fullmatch(r'\[[ATGC]/[ATGC]\]', snp) is not None

    def get_snp(row):
        seq = row.get('TopGenomicSeq')
        if not isinstance(seq, str):
            return None
        rem = re.fullmatch(r'.*\[([ATGC])/([ATGC])\].*', seq)
        if rem is None:
            return None
        return rem.group(1, 2)

    def get_loc(row):
        chrom_val = row.get('Chr', '')
        pos_val = row.get('MapInfo', '')
        if not (isinstance(chrom_val, str) and chrom_val.isdigit() and isinstance(pos_val, str) and pos_val.isdigit()):
            return None
        chrom, pos = int(chrom_val), int(pos_val)
        if chrom < 1 or chrom > 22 or pos <= 0:
            return None
        return (chrom, pos)

    def tobe_flip(s1, s2):
        return {
            ('TOP', '+'): False,
            ('BOT', '-'): False,
            ('TOP', '-'): True,
            ('BOT', '+'): True,
        }[(s1, s2)]

    rs_tbl = {}
    loc_tbl = {}
    with open(manifest_file, 'rt') as ifd:
        reader = csv.reader(ifd, delimiter=',')
        for row in reader:
            if row and row[0] == '[Assay]':
                break
        reader = csv.DictReader(ifd, delimiter=',')
        for row in reader:
            if not is_snv(row):
                continue
            snp = get_snp(row)
            if snp is None:
                continue
            loc = get_loc(row)
            if loc is None:
                continue
            flip = tobe_flip(row['IlmnStrand'], row['RefStrand'])
            if flip:
                snp = tuple(map(bc_flip, snp))

            try:
                ref = refseq.fetch(str(loc[0]), loc[1] - 1, loc[1]).upper()
            except KeyError:
                ref = refseq.fetch(f'chr{loc[0]}', loc[1] - 1, loc[1]).upper()

            if ref not in snp:
                continue

            snp = tuple(sorted(snp, key=lambda b: 0 if b == ref else 1))
            if loc not in loc_tbl:
                loc_tbl[loc] = {'rs': []}
            rs = row['Name']
            loc_tbl[loc]['rs'].append(rs)
            rs_tbl[rs] = {'loc': loc, 'flip': flip, 'snp': snp}

    for loc_info in loc_tbl.values():
        alts = []
        for rs in loc_info['rs']:
            rs_info = rs_tbl[rs]
            ref, alt = rs_info['snp']
            if alt not in alts:
                alts.append(alt)
        loc_info['snp'] = tuple([ref] + alts)

    chrom_tbl = {}
    for chrom in range(1, 23):
        try:
            length = refseq.get_reference_length(str(chrom))
        except KeyError:
            length = refseq.get_reference_length(f'chr{chrom}')
        chrom_tbl[chrom] = {'length': length}

    return (loc_tbl, rs_tbl, chrom_tbl)


def read_finalreport(*, final_report_file, loc_tbl, rs_tbl):
    with open(final_report_file, 'rt') as ifd:
        reader = csv.reader(ifd, delimiter=',')
        for row in reader:
            if row and row[0] == '[Data]':
                break
        sids = tuple(next(reader)[1:])
        num_sids = len(sids)
        num_locs = len(loc_tbl)

        # GT 用: 2 bytes / call (ASCII '..', '00', '01', etc.)
        results_gt = bytearray(b'..' * num_locs * num_sids)
        # GCSCORE 用: 2 bytes unsigned short (H) / call (0〜10000 にスケーリング、欠損は 0xFFFF)
        results_gc = bytearray(b'\xff\xff' * num_locs * num_sids)

        for row in reader:
            rs = row[0]
            if rs not in rs_tbl:
                continue
            rs_info = rs_tbl[rs]
            loc = rs_info['loc']
            if loc not in loc_tbl:
                continue

            flip = rs_info['flip']
            loc_info = loc_tbl[loc]
            loc_idx = loc_info['idx']
            loc_offset_gt = loc_idx * num_sids * 2
            loc_offset_gc = loc_idx * num_sids * 2
            snp = loc_info['snp']

            for sid_idx, call_cell in enumerate(row[1:]):
                if '|' in call_cell:
                    call, gcscore_str = call_cell.split('|', 1)
                else:
                    call, gcscore_str = call_cell, ''

                if call == '--' or len(call) < 2:
                    continue

                call_bases = (call[0], call[1])
                if flip:
                    call_bases = (bc_flip(call_bases[0]), bc_flip(call_bases[1]))

                try:
                    gt1 = ''.join(map(str, sorted(map(lambda b: snp.index(b), call_bases))))
                except ValueError:
                    continue

                # GenCall Score を 10000 倍した整数 (0〜10000) として保持
                gc_val = 0
                if gcscore_str:
                    try:
                        gc_val = int(round(float(gcscore_str) * 10000))
                    except ValueError:
                        gc_val = 0

                offset_gt = loc_offset_gt + sid_idx * 2
                offset_gc = loc_offset_gc + sid_idx * 2

                gt0 = str(results_gt[offset_gt : offset_gt + 2], encoding='ASCII')

                # 実データ検証に基づく調停ロジック
                if gt0 == '..':
                    results_gt[offset_gt : offset_gt + 2] = bytearray(gt1, encoding='ASCII')
                    struct.pack_into('<H', results_gc, offset_gc, gc_val)
                    continue
                if gt0 == gt1:
                    continue
                if gt0 != '00' and gt1 == '00':
                    continue

                results_gt[offset_gt : offset_gt + 2] = bytearray(gt1, encoding='ASCII')
                struct.pack_into('<H', results_gc, offset_gc, gc_val)

    return (sids, results_gt, results_gc)


def save_vcf(*, chrom, chrom_length, loc_tbl, sids, results_gt, results_gc, source_name, vcf_filename):
    header_lines = [
        vcfpy.HeaderLine('fileformat', 'VCFv4.2'),
        vcfpy.HeaderLine('fileDate', datetime.datetime.now(tz=ZoneInfo('Asia/Tokyo')).isoformat()),
        vcfpy.HeaderLine('source', source_name),
        vcfpy.InfoHeaderLine.from_mapping({'ID': 'AC', 'Number': 'A', 'Type': 'Integer', 'Description': 'Allele count in genotypes, for each ALT allele, in the same order as listed'}),
        vcfpy.InfoHeaderLine.from_mapping({'ID': 'AN', 'Number': '1', 'Type': 'Integer', 'Description': 'Total number of alleles in called genotypes'}),
        vcfpy.InfoHeaderLine.from_mapping({'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele frequency in genotypes, for each ALT allele, in the same order as listed'}),
        vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'}),
        vcfpy.FormatHeaderLine.from_mapping({'ID': 'GCSCORE', 'Number': 1, 'Type': 'Float', 'Description': 'Illumina GenCall Score'}),
        vcfpy.ContigHeaderLine.from_mapping({'ID': str(chrom), 'length': chrom_length}),
    ]
    header = vcfpy.Header(lines=header_lines, samples=vcfpy.SamplesInfos(sids))

    num_sids = len(sids)
    with open(vcf_filename, 'wb') as ofd:
        writer = vcfpy.Writer.from_stream(ofd, header, use_bgzf=True)

        for loc, loc_info in sorted(loc_tbl.items(), key=lambda kv: kv[1]['idx']):
            chrom_str = str(loc[0])
            pos = loc[1]
            snp_id = ';'.join(loc_info['rs']) if loc_info['rs'] else None
            ref = loc_info['snp'][0]
            alt = [vcfpy.Substitution(type_='SNV', value=a) for a in loc_info['snp'][1:]]

            loc_offset = loc_info['idx'] * num_sids * 2

            calls = []
            num_alts = len(alt)
            ac_counts = [0] * num_alts
            total_an = 0

            for idx, sid in enumerate(sids):
                off = loc_offset + idx * 2
                g1_char = chr(results_gt[off])
                g2_char = chr(results_gt[off + 1])

                # GCSCORE の取り出し
                raw_gc = struct.unpack_from('<H', results_gc, off)[0]
                gc_val = round(raw_gc / 10000.0, 4) if raw_gc != 0xFFFF else None

                call_data = {'GT': f'{g1_char}/{g2_char}'}
                if gc_val is not None:
                    call_data['GCSCORE'] = gc_val

                calls.append(vcfpy.Call(sample=sid, data=call_data))

                for g_char in (g1_char, g2_char):
                    if g_char != '.':
                        total_an += 1
                        allele_idx = int(g_char)
                        if allele_idx > 0:
                            ac_counts[allele_idx - 1] += 1

            info_dict = {'AN': total_an, 'AC': ac_counts}
            if total_an > 0:
                info_dict['AF'] = [round(ac / total_an, 6) for ac in ac_counts]
            else:
                info_dict['AF'] = [0.0] * num_alts

            record = vcfpy.Record(
                CHROM=chrom_str,
                POS=pos,
                ID=snp_id,
                REF=ref,
                ALT=alt,
                QUAL=None,
                FILTER=['PASS'],
                INFO=info_dict,
                FORMAT=['GT', 'GCSCORE'],
                calls=calls,
            )
            writer.write_record(record)

        writer.close()

    # 出力された .vcf.gz に対する tabix インデックス作成
    pysam.tabix_index(str(vcf_filename), force=True, preset='vcf')


def proc(args):
    loc_tbl_chrom, rs_tbl, chrom, chrom_length, final_report_file, source_name, vcf_filename = args
    sids, results_gt, results_gc = read_finalreport(
        final_report_file=final_report_file, loc_tbl=loc_tbl_chrom, rs_tbl=rs_tbl
    )
    save_vcf(
        chrom=chrom,
        chrom_length=chrom_length,
        loc_tbl=loc_tbl_chrom,
        sids=sids,
        results_gt=results_gt,
        results_gc=results_gc,
        source_name=source_name,
        vcf_filename=vcf_filename,
    )

def main():
    parser = argparse.ArgumentParser(description='Convert Illumina FinalReport to chromosome-separated VCF.gz files.')
    parser.add_argument('-m', '--manifest', required=True, type=Path, help='Path to Illumina manifest CSV file (e.g. ASA-24v1-0_A1.csv)')
    parser.add_argument('-r', '--refseq', required=True, type=Path, help='Path to reference FASTA file')
    parser.add_argument('-f', '--final-report', required=True, type=Path, help='Path to GenomeStudio FinalReport CSV file')
    parser.add_argument('-o', '--out-prefix', default='output', help='Prefix for output VCF files (default: output)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of parallel processes (default: 4)')
    args = parser.parse_args()

    loc_tbl, rs_tbl, chrom_tbl = read_manifest(manifest_file=args.manifest, refseq_file=args.refseq)

    # 染色体ごとに loc_tbl をあらかじめ分割してインデックスを 0 から振り直す
    tasks = []
    source_name = args.final_report.name

    for chrom in range(1, 23):
        chrom_loc_items = sorted(
            filter(lambda kv: kv[0][0] == chrom, loc_tbl.items()),
            key=lambda kv: kv[0],
        )
        if not chrom_loc_items:
            continue

        loc_tbl_chrom = {loc: loc_info | {'idx': idx} for idx, (loc, loc_info) in enumerate(chrom_loc_items)}
        vcf_filename = f'{args.out_prefix}-chr{chrom}.vcf.gz'

        tasks.append((
            loc_tbl_chrom,
            rs_tbl,
            chrom,
            chrom_tbl[chrom]['length'],
            args.final_report,
            source_name,
            vcf_filename,
        ))

    with Pool(processes=args.threads) as p:
        p.map(proc, tasks)


if __name__ == '__main__':
    main()
