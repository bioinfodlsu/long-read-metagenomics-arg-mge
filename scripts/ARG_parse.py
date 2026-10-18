#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import ast
import asyncio
import aiofiles
import io


def load_fasta_into_memory(fasta_file):
    
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def has_non_arg_regions(row):
    return len(extract_regions(row["arg_coords"], int(row["query_length"]))) > 0


def prefilter_arg_summary(df):
        
    mask = df.apply(has_non_arg_regions, axis=1)
    df_filtered = df[mask].reset_index(drop=True)
    
    skipped = len(df) - len(df_filtered)
   
    return df_filtered


def extract_regions(arg_coords, query_length):
    # Non-ARG regions of a read as (start, end), 0-based and inclusive like the ARG coordinates,
    # so a region never shares a base with an ARG:
    #   before the first ARG [0, start - 1], between ARGs [end + 1, next start - 1],
    #   after the last ARG [end + 1, query_length - 1]
    regions = []
    pos = 0  # first base not covered by any ARG so far
    for start, end in sorted(arg_coords, key=lambda x: x[0]):
        start, end = int(start), int(end)
        if start > pos:
            regions.append((pos, start - 1))
        pos = max(pos, end + 1)  # max() keeps a nested ARG from reopening a covered stretch

    if pos <= query_length - 1:
        regions.append((pos, query_length - 1))

    return regions


def process_row(row, sequences_dict):
    fasta_id = row["query_name"]
    arg_coords = row["arg_coords"]
    query_length = int(row["query_length"])
    
    if fasta_id not in sequences_dict:
        return ""
    
    original_seq = sequences_dict[fasta_id]
    
    regions = extract_regions(arg_coords, query_length)
    
    if not regions:
        return ""
    
    output = io.StringIO()
    count = 0
    
    for idx, (start, end) in enumerate(regions, 1):
        seq = original_seq[start:end+1]
        if len(seq) >= 50:
            output.write(f">{fasta_id}_part{idx} {start}:{end}\n")
            for i in range(0, len(seq), 60):
                output.write(f"{seq[i:i+60]}\n")
            count += 1
    
    return output.getvalue()


async def process_and_write_async(row, sequences_dict, output_file, semaphore):
    result = process_row(row, sequences_dict)
    
    if result:
        async with semaphore:  
            async with aiofiles.open(output_file, 'a') as f:
                await f.write(result)


async def extract_sequences_async(df, sequences_dict, output_file, max_concurrent=50):
    open(output_file, 'w').close()
    semaphore = asyncio.Semaphore(max_concurrent)
    
    tasks = [
        process_and_write_async(row, sequences_dict, output_file, semaphore)
        for _, row in df.iterrows()
    ]
    
    await asyncio.gather(*tasks)


def keep_best_by_blocks(df, touching_is_overlap=True):
    df = df.sort_values("start")
    keep = []

    cur_start = None
    cur_end = None
    cur_best_idx = None
    cur_best_eval = np.inf

    def flush():
        if cur_best_idx is not None:
            keep.append(cur_best_idx)

    for idx, row in df.iterrows():
        s, e, ev = row["start"], row["end"], row["E-value"]

        if cur_start is None:
            cur_start, cur_end = s, e
            cur_best_idx, cur_best_eval = idx, ev
            continue

        overlaps = s <= cur_end if touching_is_overlap else s < cur_end

        if overlaps:
            cur_end = max(cur_end, e)
            if ev < cur_best_eval:
                cur_best_idx, cur_best_eval = idx, ev
        else:
            flush()
            cur_start, cur_end = s, e
            cur_best_idx, cur_best_eval = idx, ev

    flush()
    return df.loc[keep]


def correct_negative_strand(df_arg):
    df_arg = df_arg.assign(
        start=np.where(
            df_arg["query_strand"] == "-",
            df_arg["query_length"] - df_arg["query_to"] -1,  # revese complement strand
            df_arg["query_from"]
        ),
        end=np.where(
            df_arg["query_strand"] == "-",
            df_arg["query_length"] - df_arg["query_from"] -1,  # revese complement strand
            df_arg["query_to"]
        ),
    )
    
    return df_arg
def parse_lastal(input_file):

    alignments = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('#'):
            i += 1
            continue
        
        if line.startswith('a '):
            a_parts = line.split()
            score = None
            evalue = None
            eg2value = None
            
            for j, part in enumerate(a_parts):
                if part.startswith('score='):
                    score = float(part.split('=')[1])
                elif part.startswith('E='):
                    evalue = float(part.split('=')[1])
                elif part.startswith('EG2='):
                    eg2value = float(part.split('=')[1])
            
            i += 1
            if i >= len(lines):
                break
            
            s_line1 = lines[i].strip()
            i += 1
            if i >= len(lines):
                break
            
            s_line2 = lines[i].strip()
            
            if s_line1.startswith('s '):
                s1_parts = s_line1.split()
                target_name = s1_parts[1]
                target_from = int(s1_parts[2])
                target_al_length = int(s1_parts[3])
                target_strand = s1_parts[4]
                target_length = int(s1_parts[5])

            
            if s_line2.startswith('s '):
                s2_parts = s_line2.split()
                query_name = s2_parts[1]
                query_from = int(s2_parts[2])
                query_al_length = int(s2_parts[3])
                query_strand = s2_parts[4]
                query_length = int(s2_parts[5])
                
            
            query_to = query_from + query_al_length -1
            target_to = target_from + target_al_length -1
            
            alignment = {
                'target_name': target_name,
                'query_name': query_name,
                'target_from': target_from,
                'target_to': target_to,
                'target_al_length': target_al_length,
                'target_length': target_length,
                'query_from': query_from,
                'query_to': query_to,
                'query_al_length': query_al_length,
                'query_length': query_length,
                'target_strand': target_strand,
                'query_strand': query_strand,
                'score': score,
                'E-value': evalue,
                'EG2_value': eg2value
            }
            
            alignments.append(alignment)
        
        i += 1
    
    df_arg = pd.DataFrame(alignments)
    df_arg = correct_negative_strand(df_arg)
    df_arg["ARO accession"] = df_arg["target_name"].str.split("|").str[2]
    df_arg["AMR gene"] = df_arg["target_name"].str.split("|").str[-1]
    df_arg["Model Name"] = df_arg["target_name"].str.split("|").str[1]
    df_arg["reference_coverage"] = df_arg["target_al_length"] / df_arg["target_length"]
    df_arg.sort_values(by="query_name").reset_index(drop=True)
    
    return df_arg

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


def compute_arg_info(df, query_name=None):
    df = df.sort_values("start").reset_index(drop=True)
    
    arg_coords = list(zip(df["start"].astype(int), df["end"].astype(int)))
    
    result = {
        "query_name": query_name,  # passed in explicitly
        "query_length": df.iloc[0]["query_length"],
        "arg_coords": arg_coords,
        "arg_accessions": df["ARO accession"].tolist(),
        "arg_evalues": df["E-value"].tolist(),
        "arg_genes": df["AMR gene"].tolist(),
        "arg_class": df["Drug Class"].tolist(),
        "arg_mechanism": df["Resistance Mechanism"].tolist(),
        "arg_strand": df["query_strand"].tolist(),
        "arg_reference_coverage": df["reference_coverage"].tolist(),
        "num_args": len(df)
    }
    return pd.DataFrame([result])

def summarize_args(df):
    df_summary = (df
                  .groupby("query_name", group_keys=False)
                  .apply(lambda g: compute_arg_info(g, query_name=g.name), include_groups=False)
                  .reset_index(drop=True))
    return df_summary

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
    # parser.add_argument("--argsum", type=str, default=None,
    #                    help="Path to a precomputed ARG summary file")
    parser.add_argument("--output", type=str, default="output_sequences.fasta",
                       help="Path to the output FASTA file")
    parser.add_argument("--metadata", type=str, default="ReferenceGeneCatalog.txt",
                       help="Path to the metadata TSV file")
    parser.add_argument("--summary", type=str, default="ARG_summary.csv",
                       help="Path to the summary output file")
    parser.add_argument("--max-concurrent", type=int, default=50,
                       help="Maximum concurrent async operations (default: 50)")
    
    args = parser.parse_args()
   
    sequences_dict = load_fasta_into_memory(args.input)
    
    
    df_arg = parse_lastal(args.maf)

    kept_indices = (df_arg
        .groupby("query_name", group_keys=False)
        .apply(keep_best_by_blocks, include_groups=False)
        .index)
    df_arg_best = df_arg.loc[kept_indices].reset_index(drop=True)
    df_arg_best = add_metadata(df_arg_best, args.metadata)
    df_arg_best.drop_duplicates(subset=["query_name", "target_name", "start", "end"], inplace=True)


    df_arg_summary = summarize_args(df_arg_best)
    
    write_summary_output(df_arg_summary, args.summary)
    
    # Pre-filter sequences with no non-ARG regions
    df_filtered = prefilter_arg_summary(df_arg_summary)
    if len(df_filtered) == 0:
        return
    
    # Extract sequences using AsyncIO
    asyncio.run(extract_sequences_async(df_filtered, sequences_dict, 
                                        args.output, args.max_concurrent))



if __name__ == "__main__":
    main()