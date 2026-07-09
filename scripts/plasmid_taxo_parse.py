#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import sys
import os
import gc
from collections import defaultdict

def process_plasmid_streaming(plasmid_file, typing_file, finder_file):
    
    df_typing = pd.read_csv(
        typing_file,
        usecols=["NUCCORE_ACC", "relaxase_type(s)", "observed_host_range_ncbi_name", "observed_host_range_ncbi_rank"]
    )
    df_typing = df_typing.drop_duplicates(subset=["NUCCORE_ACC"], keep='first')

    typing_dict = df_typing.set_index("NUCCORE_ACC").to_dict('index')
    del df_typing
    gc.collect()
    
    df_finder = pd.read_csv(finder_file, usecols=["NUCCORE_ACC", "typing"])
    df_finder = df_finder.drop_duplicates(subset=["NUCCORE_ACC"], keep='first')
    finder_dict = df_finder.set_index("NUCCORE_ACC").to_dict('index')
    del df_finder
    gc.collect()
    plasmid_summary = {}

    def process_query(query_name, alignments):
        if not alignments:
            return
        alignments.sort(key=lambda x: x['start'])
        keep_alignments = []
        cur_end = alignments[0]['start']
        cur_best = alignments[0]
        cur_best_eval = alignments[0]['eval'] if alignments[0]['eval'] is not None else float('inf')

        for aln in alignments[1:]:
            aln_eval = aln['eval'] if aln['eval'] is not None else float('inf')
            if aln['start'] <= cur_end:
                cur_end = max(cur_end, aln['end'])
                if aln_eval < cur_best_eval:
                    cur_best = aln
                    cur_best_eval = aln_eval
            else:
                keep_alignments.append(cur_best)
                cur_end = aln['end']
                cur_best = aln
                cur_best_eval = aln_eval
        keep_alignments.append(cur_best)

        plasmid_summary[query_name] = {
            'plasmid_accessions':      [a['target'] for a in keep_alignments],
            'plasmid_query_lengths':   [a['query_length'] for a in keep_alignments],
            'plasmid_query_al_lengths':[a['query_al_length'] for a in keep_alignments],
            'plasmid_coords':          [(a['start'], a['end']) for a in keep_alignments],
            'plasmid_reference_coverages':   [a['target_al_length'] / a['target_length'] for a in keep_alignments],
            'plasmid_evalues':         [a['eval'] for a in keep_alignments],
            'relaxase_types':          [typing_dict.get(a['target'], {}).get('relaxase_type(s)') for a in keep_alignments],
            'replicon_types':          [finder_dict.get(a['target'], {}).get('typing') for a in keep_alignments],
            'host_ranges':             [typing_dict.get(a['target'], {}).get('observed_host_range_ncbi_name') for a in keep_alignments],
            'host_range_ranks':        [typing_dict.get(a['target'], {}).get('observed_host_range_ncbi_rank') for a in keep_alignments],
        }

    current_query = None
    current_alignments = []
    processed = 0

    with open(plasmid_file, 'r', buffering=1024*1024) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('a '):
                eg2value = None
                for part in line.split():
                    if part.startswith(('EG=', 'EG2=')):
                        eg2value = float(part.split('=', 1)[1])
                        break

                s_line1 = f.readline().strip()
                s_line2 = f.readline().strip()
                if not s_line1 or not s_line2:
                    break

                s1_parts = s_line1.split()
                target_name    = s1_parts[1]
                target_al_length = int(s1_parts[3])
                target_length  = int(s1_parts[5])

                s2_parts = s_line2.split()
                query_name     = s2_parts[1]
                query_from     = int(s2_parts[2])
                query_al_length = int(s2_parts[3])
                query_strand   = s2_parts[4]
                query_length   = int(s2_parts[5])

                if query_strand == "-":
                    start = query_length - (query_from + query_al_length - 1) - 1
                    end   = query_length - query_from - 1
                else:
                    start = query_from
                    end   = query_from + query_al_length - 1

                if query_name != current_query:
                    if current_query is not None:
                        process_query(current_query, current_alignments)
                        processed += 1
                    current_query = query_name
                    current_alignments = []

                current_alignments.append({
                    'target': target_name,
                    'start': start, 'end': end,
                    'eval': eg2value,
                    'target_al_length': target_al_length,
                    'target_length': target_length,
                    'query_al_length': query_al_length,
                    'query_length': query_length,
                })

    # flush last query
    if current_query is not None:
        process_query(current_query, current_alignments)

    return plasmid_summary
    

def process_taxonomy_streaming(reads_file, report_file):
    
    ranks_order = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
    rank_names = {
        'D': 'Domain', 'P': 'Phylum', 'C': 'Class', 
        'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'
    }
    
    lineage_map = {}
    current_lineage = {k: None for k in rank_names.values()}
    
    with open(report_file, 'r', buffering=1024*1024) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            rank_code = parts[3]
            tax_id = int(parts[4])
            name = parts[5].strip()
            
            if rank_code in ranks_order:
                rank_name = rank_names[rank_code]
                current_lineage[rank_name] = name
                current_idx = ranks_order.index(rank_code)
                for lower_rank_code in ranks_order[current_idx + 1:]:
                    current_lineage[rank_names[lower_rank_code]] = None
            
            lineage_map[tax_id] = current_lineage.copy()
    
    taxonomy_data = {}
    
    with open(reads_file, 'r', buffering=1024*1024) as f:
        count = 0
        for line in f:
            count += 1
            
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            read_name = parts[1]
            raw_taxonomy = parts[2]
            
            if '(taxid ' in raw_taxonomy:
                try:
                    start = raw_taxonomy.index('(taxid ') + 7
                    end = raw_taxonomy.index(')', start)
                    tax_id = int(raw_taxonomy[start:end])
                    
                    lineage = lineage_map.get(tax_id, {})
                    taxonomy_data[read_name] = {
                        'taxonomy_id': tax_id,
                        'Domain': lineage.get('Domain'),
                        'Phylum': lineage.get('Phylum'),
                        'Class': lineage.get('Class'),
                        'Order': lineage.get('Order'),
                        'Family': lineage.get('Family'),
                        'Genus': lineage.get('Genus'),
                        'Species': lineage.get('Species')
                    }
                except (ValueError, IndexError):
                    continue
    
    
    del lineage_map
    gc.collect()
    
    return taxonomy_data

def merge_all_data_streaming(summary_file, plasmid_dict, taxonomy_dict, output_file):
    
    plasmid_cols = ['plasmid_accessions', 'plasmid_coords', 'plasmid_evalues', 'plasmid_reference_coverages',
                    'relaxase_types', 'replicon_types', 'host_ranges', 'host_range_ranks']
    taxonomy_cols = ['taxonomy_id', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    
    chunk_size = 50000
    first_chunk = True
    
    for chunk in pd.read_csv(summary_file, chunksize=chunk_size):        
        for col in plasmid_cols:
            chunk[col] = chunk['query_name'].map(lambda x: plasmid_dict.get(x, {}).get(col))
        
        for col in taxonomy_cols:
            chunk[col] = chunk['query_name'].map(lambda x: taxonomy_dict.get(x, {}).get(col))
        
        mode = 'w' if first_chunk else 'a'
        header = first_chunk
        chunk.to_csv(output_file, mode=mode, index=False, header=header)
        
        first_chunk = False
        
        del chunk
        gc.collect()
    

def main():
    parser = argparse.ArgumentParser(description="Process LASTAL alignments and merge metadata.")
    
    parser.add_argument("--plasmid", required=True, help="Path to LASTAL alignment file (.maf)")
    parser.add_argument("--typing", required=True, help="Path to plasmid typing metadata CSV")
    parser.add_argument("--finder", required=True, help="Path to plasmid finder metadata CSV")
    parser.add_argument("--tax-reads", required=True, help="Path to per-read taxonomy file (.tab)")
    parser.add_argument("--tax-report", required=True, help="Path to taxonomy report file (.txt)")
    parser.add_argument("--summary", required=True, help="Path to the existing summary CSV file")
    parser.add_argument("--output", "-o", required=True, help="Path to save the final merged CSV")

    args = parser.parse_args()

    
    plasmid_dict = process_plasmid_streaming(args.plasmid, args.typing, args.finder)
    
    taxonomy_dict = process_taxonomy_streaming(args.tax_reads, args.tax_report)
    
    merge_all_data_streaming(args.summary, plasmid_dict, taxonomy_dict, args.output)
    
    del plasmid_dict
    del taxonomy_dict
    gc.collect()
    

if __name__ == "__main__":
    main()