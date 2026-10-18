#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import ast
import asyncio
import aiofiles
import io
import re
import os




def correct_negative_strand(df):
    df = df.assign(
        start=np.where(
            df["query_strand"] == "-",
            df["query_length"] - df["query_to"]-1,  # revese complement strand
            df["query_from"]
        ),
        end=np.where(
            df["query_strand"] == "-",
            df["query_length"] - df["query_from"]-1,  # revese complement strand
            df["query_to"]
        ),
    )
    
    return df

def read_fragment_offsets(header_txt_path):
    # {fragment name: start on the read} from the "{read}_partN start:end" headers ARG_parse.py writes.
    # Lines may come with or without the leading '>' (the Snakefile strips it).
    offsets = {}
    with open(header_txt_path, 'r') as header_file:
        for line in header_file:
            header_parts = line.strip().lstrip('>').split()
            if len(header_parts) < 2:
                continue
            coords = header_parts[1].split(":")
            if len(coords) != 2:
                continue
            try:
                offsets[header_parts[0]] = int(coords[0])
            except ValueError:
                continue
    return offsets


def correct_for_length(df, header_txt_path, start_col="start", end_col="end"):
    # Fragment coordinates -> read coordinates: a hit at fragment position k sits at read base start + k.
    # Every fragment must have an offset; a missing one would leave its hits in fragment coordinates.
    offsets = read_fragment_offsets(header_txt_path)
    missing = sorted(set(df["query_name"]) - set(offsets))
    if missing:
        raise ValueError(
            f"{len(missing)} of {df['query_name'].nunique()} aligned fragments have no 'start:end' offset in "
            f"{header_txt_path} (e.g. {missing[:3]}); check that it lists the headers of the fragment FASTA "
            f"written by ARG_parse.py for this sample"
        )
    offset = df["query_name"].map(offsets)
    df[start_col] = df["start"] + offset
    df[end_col] = df["end"] + offset
    return df

def extract_info(target_name, type, fasta_file = None):
    if pd.isna(target_name):
        return "unknown"

    if type =="tn3":
        name = str(target_name).lower()

        if "transposase" in name or "tnp" in name:
            return "transposase"
        elif "resolvase" in name or "tnr" in name or "tnpR" in name or "res" in name:
            return "resolvase"
        else:
            return "unknown"
    elif type == "tn7":
        name = str(target_name).lower()
        if "tnsa" in name:
            return "tnsA"
        elif "tnsb" in name or "tnib" in name:
            return "tnsB"
        elif "tnsc" in name:
            return "tnsC"
        elif "tnsd" in name or "tniq" in name:
            return "tnsD"
        elif "tnse" in name:
            return "tnsE"
        else:
            return "unknown"
    elif type == "tn554":
        name = str(target_name).lower()
        if "tnpa" in name:
            return "tnpA"
        elif "tnpb" in name or "tnpe" in name:
            return "tnpB"
        elif "tnpc" in name or "tnpf" in name or "trac" in name:
            return "tnpC"
        else:
            return "unknown"

    elif type == "integron":
        target_name = str(target_name)

        with open(fasta_file) as f:
            for line in f:
                if target_name in line:
                    match = re.search(r'gene=([^|\s]+)', line)
                    if match:
                        return match.group(1)
        return "unknown"
    else:
        return "unknown"

def classify_integrase(integrase_name):
    if pd.isna(integrase_name):
        return "unknown"

    name = str(integrase_name).lower()

    if "inti1" in name:
        return "intI1"
    elif "inti2" in name:
        return "intI2"
    elif "inti3" in name:
        return "intI3"
    elif name == "int":
        return "intI"
    else:
        return integrase_name
    
def parse_lastal(input_file, fasta_path, type, integrase_headers=None):
    if not input_file or not os.path.exists(input_file):
        return pd.DataFrame()

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
                elif part.startswith(('EG=', 'EG2=')):
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
            
            # Store alignment info
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
                'EG2-value': eg2value
            }
            
            alignments.append(alignment)
        
        i += 1
    
    if not alignments:
        return pd.DataFrame()
    
    df = pd.DataFrame(alignments)
    df = correct_negative_strand(df)
    df.sort_values(by="query_name").reset_index(drop=True)
    df["read"] = df["query_name"].str.split('_').str[0]
    df = correct_for_length(df, fasta_path, start_col=f"{type}_start", end_col=f"{type}_end")
    df["function"] = df["target_name"].apply(extract_info, args=(type, integrase_headers))
    if type == "integron":
        df["integrase"] = df["function"].apply(classify_integrase)
        df= df[(df['integrase'] != 'intron') & (df['integrase'] != 'fused')]
    df["reference_coverage"] = df["target_al_length"] / df["target_length"]
    return df

def keep_best_by_blocks(df, start_col="start", end_col="end"):
    if df.empty:
        return df
    
    df = df.sort_values(start_col)
    keep = []

    cur_start = None
    cur_end = None
    cur_best_idx = None
    cur_best_eval = np.inf

    def flush():
        if cur_best_idx is not None:
            keep.append(cur_best_idx)

    for idx, row in df.iterrows():
        s, e, ev = row[start_col], row[end_col], row["EG2-value"]

        if cur_start is None:
            cur_start, cur_end = s, e
            cur_best_idx, cur_best_eval = idx, ev
            continue

        overlaps = s <= cur_end

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

def compute_info(g, type, start_col="start", end_col="end", read=None):
    g = g.sort_values(start_col).reset_index(drop=True)
    
    coords = []
    genes = []
    names = []
    evalues = []
    functions = []
    ref_coverages = []

    for _, row in g.iterrows():
        coords.append((row[start_col], row[end_col]))
        evalues.append(row["EG2-value"])
        genes.append(row["target_name"])
        ref_coverages.append(row["reference_coverage"])
        if type != "integron":
            names.append(extract_info(row["target_name"], type))
            functions.append(row["function"])
        else:
            names.append(row["function"])
            functions.append(row["integrase"])
    
    # Create single row with aggregated info
    data = {
        "read": read if read is not None else g.iloc[0]["read"],
        f"{type}_coords": coords,
        f"{type}_genes": genes,
        f"{type}_names": names,
        f"{type}_evalues": evalues,
        f"{type}_reference_coverages": ref_coverages,
        f"{type}_functions": functions
    }

    # if type != "integron":
    #     data[f"{type}_functions"] = functions

    # Return as a Series (not DataFrame) so groupby.apply returns a DataFrame with one row per group
    return pd.Series(data)


def parse_column(col_str):
    try:
        if pd.isna(col_str):
            return []
        parsed = ast.literal_eval(col_str)
        return parsed if isinstance(parsed, list) else [parsed]
    except:
        return []
    
def read_arg_is_csv(arg_csv_path):
    df_arg_is = pd.read_csv(arg_csv_path)

    df_arg_is['arg_coords'] = df_arg_is['arg_coords'].apply(parse_column)
    
    for col in ['arg_accessions', 'arg_evalues', 'arg_genes', 'arg_class', 'arg_mechanism', 'arg_reference_coverage', 'IS_coords', 'IS_accessions', 'IS_evalues', 'IS_family', 'IS_names', 'IS_reference_coverages',
                'Before_ARG', 'After_ARG', 'Between_ARG', 'before_IS_fam', 'after_IS_fam', 'between_IS_fam', 'before_IS_name', 'after_IS_name', 'between_IS_name']:
        if col in df_arg_is.columns:
            df_arg_is[col] = df_arg_is[col].apply(parse_column)
    
    return df_arg_is


def merge_tn_summaries(df, summaries):
    df_merged = df.copy()

    for df_tn, prefix, cols in summaries:
        if df_tn.empty:
            continue
        merge_cols = ['read'] + cols
        df_merged = pd.merge(
            df_merged,
            df_tn[merge_cols],
            left_on='query_name',
            right_on='read',
            how='left'
        )
        if 'read' in df_merged.columns:
            df_merged = df_merged.drop(columns=['read'])

    return df_merged

def main():
    parser = argparse.ArgumentParser(description="Parse LAST alignment output for transposon/integron detection.")
    parser.add_argument("--input_tn3", help="Path to LAST alignment output file for tn3.", default=None)
    parser.add_argument("--input_tn7", help="Path to LAST alignment output file for tn7.", default=None)
    parser.add_argument("--input_tn554", help="Path to LAST alignment output file for tn554.", default=None)
    parser.add_argument("--input_integron", help="Path to LAST alignment output file for integrons.")
    parser.add_argument("--integrase_headers", help="Path to integrase headers file.", default=None)
    parser.add_argument('--fasta_headers', type=str, required=True, help='Path to fasta headers text file')
    parser.add_argument('--output', type=str, required=True, help='Path to output CSV file for merged results')
    parser.add_argument('--arg_is', type=str, required=False, help='Path to ARG and IS alignment output file')

    args = parser.parse_args()

    df_tn3 = parse_lastal(args.input_tn3, args.fasta_headers, "tn3", None)
    df_tn7 = parse_lastal(args.input_tn7, args.fasta_headers, "tn7", None)
    df_tn554 = parse_lastal(args.input_tn554, args.fasta_headers, "tn554", None)
    df_integron = parse_lastal(args.input_integron, args.fasta_headers, "integron", args.integrase_headers)

    df_tn3 = (
        df_tn3.groupby("query_name", group_keys=False)
        .apply(lambda g: keep_best_by_blocks(g, start_col="tn3_start", end_col="tn3_end"))
        .reset_index(drop=True)
    ) if not df_tn3.empty else df_tn3

    df_tn7 = (
        df_tn7.groupby("query_name", group_keys=False)
        .apply(lambda g: keep_best_by_blocks(g, start_col="tn7_start", end_col="tn7_end"))
        .reset_index(drop=True)
    ) if not df_tn7.empty else df_tn7

    df_tn554 = (
        df_tn554.groupby("query_name", group_keys=False)
        .apply(lambda g: keep_best_by_blocks(g, start_col="tn554_start", end_col="tn554_end"))
        .reset_index(drop=True)
    ) if not df_tn554.empty else df_tn554

    df_integron = (
        df_integron.groupby("query_name", group_keys=False)
        .apply(lambda g: keep_best_by_blocks(g, start_col="integron_start", end_col="integron_end"))
        .reset_index(drop=True)
    ) if not df_integron.empty else df_integron
    # # Select best function per read (skip if empty)
    # df_tn3 = (
    #     df_tn3
    #     .sort_values("EG2-value", ascending=False)
    #     .groupby(["query_name", "function"], as_index=False)
    #     .head(1)
    # ) if not df_tn3.empty else pd.DataFrame()
    
    # df_tn7 = (
    #     df_tn7
    #     .sort_values("EG2-value", ascending=False)
    #     .groupby(["query_name", "function"], as_index=False)
    #     .head(1)
    # ) if not df_tn7.empty else pd.DataFrame()
    
    # df_tn554 = (
    #     df_tn554
    #     .sort_values("EG2-value", ascending=False)
    #     .groupby(["query_name", "function"], as_index=False)
    #     .head(1)
    # ) if not df_tn554.empty else pd.DataFrame()
    
    # df_integron_best_function = (
    #     df_integron
    #     .sort_values("EG2-value", ascending=False)
    #     .groupby(["query_name", "integrase"], as_index=False)
    #     .head(1)
    # ) if not df_integron.empty else pd.DataFrame()
    
    # print(df_tn3.columns)

    df_tn3_summary = (
        df_tn3.groupby("read", group_keys=False)
        .apply(lambda g: compute_info(g, type="tn3", start_col="tn3_start", end_col="tn3_end", read = g.name), include_groups=False)
        .reset_index(drop=True)
    ) if not df_tn3.empty else pd.DataFrame()
    
    df_tn7_summary = (
        df_tn7.groupby("read", group_keys=False)
        .apply(lambda g: compute_info(g, type="tn7", start_col="tn7_start", end_col="tn7_end", read = g.name), include_groups=False)
        .reset_index(drop=True) 
    ) if not df_tn7.empty else pd.DataFrame()
    
    df_tn554_summary = (
        df_tn554.groupby("read", group_keys=False)
        .apply(lambda g: compute_info(g, type="tn554", start_col="tn554_start", end_col="tn554_end", read = g.name), include_groups=False)
        .reset_index(drop=True)
    ) if not df_tn554.empty else pd.DataFrame()
    
    df_integron_summary = (
        df_integron.groupby("read", group_keys=False)
        .apply(lambda g: compute_info(g, type="integron", start_col="integron_start", end_col="integron_end", read = g.name), include_groups=False)
        .reset_index(drop=True)
    ) if not df_integron.empty else pd.DataFrame()
    

    df_arg_is = read_arg_is_csv(args.arg_is) if args.arg_is else pd.DataFrame()
    # print (df_tn3_summary.columns)
    # print (df_integron_summary.columns)

    summaries = [
        (df_tn3_summary, 'tn3', ['tn3_coords','tn3_functions', 'tn3_genes','tn3_evalues', 'tn3_reference_coverages']),
        (df_tn7_summary, 'tn7', ['tn7_coords','tn7_functions', 'tn7_genes', 'tn7_evalues', 'tn7_reference_coverages']),
        (df_tn554_summary, 'tn554', ['tn554_coords','tn554_functions', 'tn554_genes', 'tn554_evalues', 'tn554_reference_coverages']),
        (df_integron_summary, 'integrase', ['integron_coords', 'integron_functions', 'integron_genes','integron_evalues', 'integron_reference_coverages'])
    ]

    df_merged = merge_tn_summaries(df_arg_is, summaries)
    df_merged.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()