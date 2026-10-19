#!/usr/bin/env python3
"""Value-preserving submission exports. This does not rerun any analysis."""
import argparse
import csv
import hashlib
import math
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_table(path, delimiter=','):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.reader(stream, delimiter=delimiter))
    header, rows = rows[0], rows[1:]
    assert all(len(row) == len(header) for row in rows), path
    assert len(set(header)) == len(header), path
    return header, rows


def write_table(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


class Exporter:
    def __init__(self, root, assay, overwrite):
        self.root = root
        self.out = root / assay / 'processed'
        self.audit = root / '00_inventory' / (assay + '_export_audit')
        self.overwrite = overwrite
        self.records = []
        self.columns = []
        with (root / '00_inventory/raw_data_manifest.tsv').open() as stream:
            self.sample_ids = {r['sample_id'] for r in csv.DictReader(stream, delimiter='\t')}

    def export(self, source, name, annotations, samples=(), delimiter=',',
               headerless=None, exclude_unmapped=False, note=''):
        if headerless:
            with source.open() as stream:
                rows = list(csv.reader(stream, delimiter=delimiter))
            header = headerless
            assert all(len(r) == len(header) for r in rows), source
        else:
            header, rows = read_table(source, delimiter)
        pairs = list(annotations) + list(samples)
        indices = [header.index(old) for old, new in pairs]
        labels = [new for old, new in pairs]
        assert len(set(labels)) == len(labels)
        assert set(new for old, new in samples) <= self.sample_ids
        selected, seen, skipped = [], set(), 0
        for row in rows:
            values = [row[i] for i in indices]
            if exclude_unmapped and values[0] == '*':
                skipped += 1
                continue
            assert values[0] and values[0] not in seen, (source, values[0])
            seen.add(values[0])
            for value in values[len(annotations):]:
                assert math.isfinite(float(value)) and float(value) >= 0, (source, value)
            selected.append(values)
        dest = self.out / name
        if dest.exists() and not self.overwrite:
            raise FileExistsError(dest)
        write_table(dest, labels, selected)
        got_header, got_rows = read_table(dest, '\t')
        assert got_header == labels and got_rows == selected, dest
        self.records.append([name, str(source), digest(source), digest(dest),
                             len(selected), len(samples), skipped, note])
        self.columns.extend([name, old, new] for old, new in pairs)
        print(name, len(selected), 'rows', len(samples), 'samples', flush=True)

    def finish(self):
        write_table(self.audit / 'files.tsv',
                    ['file', 'source_file', 'source_sha256', 'output_sha256',
                     'rows', 'samples', 'unmapped_rows_excluded', 'notes'], self.records)
        write_table(self.audit / 'columns.tsv',
                    ['file', 'source_column', 'submission_column'], self.columns)


def export_small(root, overwrite):
    e = Exporter(root, '02_GEO_smallRNAseq', overwrite)
    base = Path('/wanglab1/gangx/data/LYW')
    times = (0, 4, 8, 12, 24, 72)
    configs = [('DROSHA', 'R', '240704', '240704', '240705'),
               ('DGCR8', 'G', '240722', '240722', '240722'),
               ('DICER1', '', '250522', '250523', '250523')]
    for target, prefix, date, count_date, srna_date in configs:
        raw_map, norm_map = [], []
        for t in times:
            for rep in (1, 2):
                sid = f'small_{target}_{t}h_rep{rep}'
                raw = f'{prefix}{t}-{rep}' if prefix else f'{t}h-{rep}'
                norm = f'{prefix}{t}.{rep}' if prefix else f'X{t}h.{rep}'
                raw_map.append((raw, sid))
                norm_map.append((norm, sid))
        for category in ('mature', 'hairpin', 'small_ncRNA'):
            source = base / 'microRNA-seq' / (
                f'{srna_date}_miRNA_snoRNA.csv' if category == 'small_ncRNA'
                else f'{count_date}_{category}-miRNA.csv')
            header, _ = read_table(source)
            ann = [('', 'feature_id')]
            ann += [(x, x) for x in ('chromosome', 'geneID', 'geneType',
                                    'transcriptType', 'geneSymbol', 'Length') if x in header]
            e.export(source, f'{target}_{category}_raw_counts.tsv', ann, raw_map,
                     note='Original weighted counts; decimal values retained. Selected manuscript samples only.')
            result = base / date / 'result'
            norm = ('sRNA_normalized_counts.csv' if category == 'small_ncRNA'
                    else ('mature_normalized_counts-2.csv'
                          if target == 'DICER1' and category == 'mature'
                          else f'{category}_normalized_counts.csv'))
            if category != 'small_ncRNA' or target == 'DICER1':
                e.export(result / norm, f'{target}_{category}_normalized_counts.tsv',
                         [('', 'feature_id')], norm_map,
                         note='Existing normalized values, not recomputed; original expression filtering retained.')
            for t in times[1:]:
                if target == 'DICER1':
                    suffix = {'mature': 'mature-2', 'hairpin': 'hairpin',
                              'small_ncRNA': 'sRNA'}[category]
                    filename = f'{t}h_vs_0h_{suffix}_DESeq2.csv'
                else:
                    suffix = '' if category == 'small_ncRNA' else category + '_'
                    filename = f'{prefix}{t}_vs_{prefix}0_{suffix}DESeq2.csv'
                source = result / filename
                header, _ = read_table(source)
                id_col = 'Row.names' if 'Row.names' in header else ''
                ann = [(id_col, 'feature_id')]
                ann += [(x, x) for x in header if x and x != id_col]
                e.export(source, f'{target}_{t}h_vs_0h_{category}_DESeq2.tsv', ann,
                         note=('Existing statistics and NA values retained. DICER1 original model included 1/2 h; '
                               'these time points are outside the selected submission scope.'
                               if target == 'DICER1' else 'Existing statistics and NA values retained.'))
    basal_names = ['WT_WT-1', 'WT_WT-2', 'G_G13-1', 'G_G13-2',
                   'I_I13-1', 'I_I13-2', 'R_R1-1', 'R_R1-2']
    basal_ids = [f'small_basal_{target}_rep{rep}'
                 for target in ('WT', 'DGCR8', 'DICER1', 'DROSHA') for rep in (1, 2)]
    e.export(Path('/home/xinhuil/lyw-GROseq/260227/merged_miRNA_mature_2.txt'),
             'basal_mature_raw_counts.tsv', [('feature_id', 'feature_id'), ('Length', 'Length')],
             list(zip(basal_names, basal_ids)), delimiter='\t',
             headerless=['feature_id', 'Length'] + basal_names, exclude_unmapped=True,
             note='Original samtools idxstats mapped-alignment counts; not weighted timecourse counts. Terminal * row excluded.')
    e.finish()


def gro_mapping(target):
    names = (['G0_G0-1', 'G0_G0-2', 'G4_G4-1', 'G4_G4-2', 'G8-1', 'G8-2', 'G12-1', 'G12-2']
             if target == 'DGCR8' else
             ['R0-1', 'R0-4-2_R0-2', 'R4-1', 'R0-4-2_R4-2',
              'R8_R8-1', 'R8_R8-2', 'R12_R12-1', 'R12_R12-2'])
    ids = [f'GRO_{target}_{t}h_rep{rep}' for t in (0, 4, 8, 12) for rep in (1, 2)]
    return list(zip(names, ids))


def export_gro(root, overwrite):
    e = Exporter(root, '03_GEO_GROseq', overwrite)
    for target in ('DGCR8', 'DROSHA'):
        gene_tables = []
        for file, label in [('countTable.txt', 'raw_counts'),
                            ('homer_genes_norm.txt', 'normalized_counts')]:
            source = root / '00_inventory/source_tables' / f'GRO_{target}' / file
            header, rows = read_table(source, '\t')
            by_gene = {r[7].split('|')[0]: r for r in rows}
            assert len(by_gene) == len(rows), 'Ambiguous condensed gene symbols'
            gene_tables.append(by_gene)
            names = {Path(col.split(' ')[0]).name: col for col in header[8:]}
            assert len(names) == 8
            samples = [(names[name], sid) for name, sid in gro_mapping(target)]
            e.export(source, f'{target}_GRO_gene_{label}.tsv',
                     [(header[0], 'transcript_id')] + [(h, h) for h in header[1:8]],
                     samples, delimiter='\t',
                     note='HOMER mm10 condensed gene counts. Normalized table is per 1e6 mapped reads, not TPM. Columns matched by name.')
        raw, norm = gene_tables
        assert raw.keys() == norm.keys(), target
        write_table(e.audit / f'{target}_raw_normalized_gene_crosswalk.tsv',
                    ['gene_symbol', 'raw_representative_transcript_id',
                     'normalized_representative_transcript_id', 'same_coordinate_annotation'],
                    [[gene, raw[gene][0], norm[gene][0],
                      str(raw[gene][1:8] == norm[gene][1:8])]
                     for gene in sorted(raw)])
    e.finish()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--assay', choices=('small', 'gro', 'all'), default='all')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if args.assay in ('small', 'all'):
        export_small(args.root, args.overwrite)
    if args.assay in ('gro', 'all'):
        export_gro(args.root, args.overwrite)
