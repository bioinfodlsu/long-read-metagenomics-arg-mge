#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd


def read_fasta(fasta_file):
    name, chunks = None, []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                header = line[1:].split(None, 1)
                name, chunks = (header[0] if header else ''), []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(chunks)


def extract_regions(arg_coords, query_length):
    regions = []
    pos = 0
    for start, end in sorted(arg_coords, key=lambda x: x[0]):
        start, end = int(start), int(end)
        if start > pos:
            regions.append((pos, start - 1))
        pos = max(pos, end + 1)

    if pos <= query_length - 1:
        regions.append((pos, query_length - 1))

    return regions


def format_fragments(fasta_id, seq, regions):
    out = []
    for idx, (start, end) in enumerate(regions, 1):
        frag = seq[start:end + 1]
        if len(frag) >= 50:
            out.append(f">{fasta_id}_part{idx} {start}:{end}\n")
            for i in range(0, len(frag), 60):
                out.append(f"{frag[i:i + 60]}\n")
    return ''.join(out)


def extract_sequences(df_summary, fasta_file, output_file):
    todo = {}
    for name, coords, length in zip(df_summary["query_name"], df_summary["arg_coords"],
                                    df_summary["query_length"]):
        regions = extract_regions(coords, int(length))
        if regions:
            todo[name] = regions

    if not todo:
        return

    written = set()
    with open(output_file, 'w') as out:
        for name, seq in read_fasta(fasta_file):
            regions = todo.get(name)
            if regions is None or name in written:
                continue
            written.add(name)
            out.write(format_fragments(name, seq, regions))


def parse_lastal(input_file):
    cols = {k: [] for k in ('target_name', 'query_name', 'target_from', 'target_al_length',
                            'target_length', 'query_from', 'query_al_length', 'query_length',
                            'target_strand', 'query_strand', 'score', 'E-value', 'EG2_value')}

    with open(input_file, 'r') as f:
        for line in f:
            if not line.startswith('a '):
                continue

            score = evalue = eg2value = None
            for part in line.split():
                if part.startswith('score='):
                    score = float(part.split('=')[1])
                elif part.startswith('E='):
                    evalue = float(part.split('=')[1])
                elif part.startswith('EG2='):
                    eg2value = float(part.split('=')[1])

            s1 = next(f, '').strip()
            s2 = next(f, '').strip()
            if not (s1.startswith('s ') and s2.startswith('s ')):
                continue
            s1_parts, s2_parts = s1.split(), s2.split()

            cols['target_name'].append(s1_parts[1])
            cols['target_from'].append(int(s1_parts[2]))
            cols['target_al_length'].append(int(s1_parts[3]))
            cols['target_strand'].append(s1_parts[4])
            cols['target_length'].append(int(s1_parts[5]))
            cols['query_name'].append(s2_parts[1])
            cols['query_from'].append(int(s2_parts[2]))
            cols['query_al_length'].append(int(s2_parts[3]))
            cols['query_strand'].append(s2_parts[4])
            cols['query_length'].append(int(s2_parts[5]))
            cols['score'].append(score)
            cols['E-value'].append(evalue)
            cols['EG2_value'].append(eg2value)

    df_arg = pd.DataFrame(cols)
    df_arg['E-value'] = df_arg['E-value'].astype(float)
    df_arg['query_to'] = df_arg['query_from'] + df_arg['query_al_length'] - 1
    df_arg['target_to'] = df_arg['target_from'] + df_arg['target_al_length'] - 1
    df_arg = correct_negative_strand(df_arg)
    split = df_arg["target_name"].str.split("|")
    df_arg["ARO accession"] = split.str[2]
    df_arg["AMR gene"] = split.str[-1]
    df_arg["Model Name"] = split.str[1]
    df_arg["reference_coverage"] = df_arg["target_al_length"] / df_arg["target_length"]

    return df_arg


def correct_negative_strand(df_arg):
    df_arg = df_arg.assign(
        start=np.where(
            df_arg["query_strand"] == "-",
            df_arg["query_length"] - df_arg["query_to"] - 1,
            df_arg["query_from"]
        ),
        end=np.where(
            df_arg["query_strand"] == "-",
            df_arg["query_length"] - df_arg["query_from"] - 1,
            df_arg["query_to"]
        ),
    )

    return df_arg


def keep_best_by_blocks(df, touching_is_overlap=True):
    codes, _ = pd.factorize(df["query_name"], sort=True)
    by_read = np.argsort(codes, kind="stable")
    bounds = np.flatnonzero(np.diff(codes[by_read])) + 1
    firsts = np.concatenate(([0], bounds)).tolist()
    lasts = np.concatenate((bounds, [len(by_read)])).tolist()

    starts = df["start"].to_numpy(dtype=np.int64)
    ends = df["end"].to_numpy()
    evalues = df["E-value"].to_numpy(dtype=float)

    keep = []
    for a, b in zip(firsts, lasts):
        rows = by_read[a:b]
        if b - a == 1:
            keep.append(rows[0])
            continue
        rows = rows[np.argsort(starts[rows], kind="quicksort")]

        best = rows[0]
        cur_end, best_ev = ends[best], evalues[best]
        for r in rows[1:]:
            overlaps = starts[r] <= cur_end if touching_is_overlap else starts[r] < cur_end
            if overlaps:
                cur_end = max(cur_end, ends[r])
                if evalues[r] < best_ev:
                    best, best_ev = r, evalues[r]
            else:
                keep.append(best)
                best, cur_end, best_ev = r, ends[r], evalues[r]
        keep.append(best)

    return df.iloc[keep]


def add_metadata(df_arg, metadata_file):
    df_metadata = pd.read_csv(metadata_file, sep="\t")
    df_metadata["Drug Class"] = df_metadata["Drug Class"].str.replace("antibiotic", "", regex=False)
    df_arg_with_metadata = df_arg.merge(
            df_metadata[["ARO Accession", "AMR Gene Family", "Drug Class", "Resistance Mechanism", "CARD Short Name"]],
            left_on="ARO accession",
            right_on="ARO Accession",
            how="left"
        ).drop(columns=["ARO Accession"])
    return df_arg_with_metadata


def summarize_args(df):
    df = df.sort_values(["query_name", "start"], kind="mergesort").reset_index(drop=True)
    names = df["query_name"].to_numpy()
    if len(df) == 0:
        return pd.DataFrame(columns=["query_name", "query_length", "arg_coords", "arg_accessions",
                                     "arg_evalues", "arg_genes", "arg_class", "arg_mechanism",
                                     "arg_strand", "arg_reference_coverage", "num_args"])

    bounds = np.flatnonzero(names[1:] != names[:-1]) + 1
    first = np.concatenate(([0], bounds))
    last = np.concatenate((bounds, [len(df)]))

    def per_read(col):
        values = df[col].tolist()
        return [values[a:b] for a, b in zip(first, last)]

    starts = per_read("start")
    ends = per_read("end")
    return pd.DataFrame({
        "query_name": names[first],
        "query_length": df["query_length"].to_numpy()[first],
        "arg_coords": [list(zip(map(int, st), map(int, en))) for st, en in zip(starts, ends)],
        "arg_accessions": per_read("ARO accession"),
        "arg_evalues": per_read("E-value"),
        "arg_genes": per_read("AMR gene"),
        "arg_class": per_read("Drug Class"),
        "arg_mechanism": per_read("Resistance Mechanism"),
        "arg_strand": per_read("query_strand"),
        "arg_reference_coverage": per_read("reference_coverage"),
        "num_args": (last - first),
    })


def write_summary_output(df, output_file):
    """Write ARG summary to a CSV file."""
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Extract non-ARG regions from sequences based on ARG table data"
    )
    parser.add_argument("--input", type=str, help="Path to the input FASTA file")
    parser.add_argument("--maf", type=str, default="sample_ARG.maf",
                       help="Path to the ARG maf file")
    parser.add_argument("--output", type=str, default="output_sequences.fasta",
                       help="Path to the output FASTA file")
    parser.add_argument("--metadata", type=str, default="ReferenceGeneCatalog.txt",
                       help="Path to the metadata TSV file")
    parser.add_argument("--summary", type=str, default="ARG_summary.csv",
                       help="Path to the summary output file")
    parser.add_argument("--max-concurrent", type=int, default=50,
                       help="Maximum concurrent async operations (default: 50)")

    args = parser.parse_args()

    df_arg = parse_lastal(args.maf)
    df_arg_best = keep_best_by_blocks(df_arg).reset_index(drop=True)
    del df_arg
    df_arg_best = add_metadata(df_arg_best, args.metadata)
    df_arg_best.drop_duplicates(subset=["query_name", "target_name", "start", "end"], inplace=True)

    df_arg_summary = summarize_args(df_arg_best)
    del df_arg_best

    write_summary_output(df_arg_summary, args.summary)

    extract_sequences(df_arg_summary, args.input, args.output)


if __name__ == "__main__":
    main()
