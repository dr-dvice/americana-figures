# Version 2.0 - dual whole/ntop alignment score columns
# Output columns: cluster, new_name, target_types, whole_alignment_score,
#                 ntop_alignment_score, source_dataset, num_alignments, status
# NOTE: interactive resume requires a v2 (whole/ntop) CSV; convert older outputs
#       first with migrate_output_format.py.
import pandas as pd
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple


# Map a filename token to the source paper abbreviation shown in source_dataset.
# Edit/extend this for other target datasets.
PAPER_ABBREV = {
    'fly2018': 'davie2018',
    'fly2020': 'kurman2020',
    'fly2021': 'ozel2021',
}


def paper_from_source(filename: str) -> str:
    """Convert a source filename to its paper abbreviation (unchanged if no token matches)."""
    for token, abbrev in PAPER_ABBREV.items():
        if token in filename:
            return abbrev
    return filename


def find_common_prefix(strings: List[str]) -> str:
    """Find the longest common prefix among a list of strings"""
    if not strings:
        return ""

    # Start with the first string as the prefix
    prefix = strings[0]

    for string in strings[1:]:
        # Reduce the prefix until it matches the start of the current string
        while not string.startswith(prefix):
            prefix = prefix[:-1]
            if not prefix:
                return ""

    return prefix


def colorize_score(score: float, text: str) -> str:
    """Apply color gradient to text based on score (0.0-1.0)"""
    # ANSI color codes
    RED = '\033[91m'
    YELLOW = '\033[93m'
    GREEN = '\033[92m'
    RESET = '\033[0m'

    if score >= 0.7:
        return f"{GREEN}{text}{RESET}"
    elif score >= 0.4:
        return f"{YELLOW}{text}{RESET}"
    else:
        return f"{RED}{text}{RESET}"


def load_alignment_data(file_paths: List[str], args) -> Dict[str, pd.DataFrame]:
    """Load and validate multiple alignment files"""
    dataframes = {}

    for file_path in file_paths:
        path = Path(file_path)
        if not path.exists():
            print(f"Error: Input file '{file_path}' not found", file=sys.stderr)
            sys.exit(1)

        try:
            df = pd.read_csv(file_path, index_col=0)
            dataframes[path.stem] = df
            if not args.quiet:
                print(f"Loaded {path.name}: {df.shape}")
        except Exception as e:
            print(f"Error reading {file_path}: {e}", file=sys.stderr)
            sys.exit(1)

    return dataframes

def get_cluster_candidates(dataframes: Dict[str, pd.DataFrame], cluster_id: str, args) -> Tuple[List[Dict], List[Dict]]:
    """Get top candidates for a cluster from all dataframes, separated by file type"""
    normal_candidates = []
    ntop_candidates = []

    for file_name, df in dataframes.items():
        source_rows = [idx for idx in df.index if str(idx).startswith(args.source_prefix)]
        target_columns = [col for col in df.columns if col.startswith(args.target_prefix)]

        if cluster_id not in source_rows:
            continue

        candidates = []
        for target_col in target_columns:
            score = df.loc[cluster_id, target_col]
            if pd.notna(score) and score > 0:
                target_type = target_col.replace(args.target_prefix, '')
                candidates.append({
                    'target_type': target_type,
                    'score': float(score),
                    'file': file_name
                })

        # Sort by score and take top candidates
        candidates.sort(key=lambda x: x['score'], reverse=True)
        top_candidates = candidates[:5]  # Top 5 from each file

        # Categorize by file type (accept 'top', 'ntop', etc.)
        if 'top' in file_name.lower():
            ntop_candidates.extend(top_candidates)
        else:
            normal_candidates.extend(top_candidates)

    # Sort both lists by score
    normal_candidates.sort(key=lambda x: x['score'], reverse=True)
    ntop_candidates.sort(key=lambda x: x['score'], reverse=True)

    return normal_candidates, ntop_candidates

def get_label_scores(cluster_label: str, target_type: str,
                     normal_candidates: List[Dict], ntop_candidates: List[Dict],
                     quiet: bool = False) -> Tuple:
    """Look up both the whole and ntop alignment scores for a target type.

    Returns (whole_score, ntop_score); either may be None if the target type is
    absent from that file type. A valid label should have both scores, so when
    exactly one side is missing a warning is printed (unless quiet).
    """
    whole = next((c['score'] for c in normal_candidates if c['target_type'] == target_type), None)
    ntop = next((c['score'] for c in ntop_candidates if c['target_type'] == target_type), None)

    if not quiet and (whole is None) != (ntop is None):
        missing = 'ntop' if ntop is None else 'whole'
        present = 'whole' if missing == 'ntop' else 'ntop'
        print(f"WARNING: {cluster_label} label '{target_type}' has no {missing} score "
              f"(only {present} present); writing NA for {missing}.", file=sys.stderr)

    return whole, ntop


def _scores_for_types(cluster_label: str, target_types: List[str],
                      normal_candidates: List[Dict], ntop_candidates: List[Dict],
                      quiet: bool = False) -> Tuple[List, List]:
    """Build parallel (whole_scores, ntop_scores) lists for a list of target types."""
    whole_scores, ntop_scores = [], []
    for t in target_types:
        w, n = get_label_scores(cluster_label, t, normal_candidates, ntop_candidates, quiet)
        whole_scores.append(w)
        ntop_scores.append(n)
    return whole_scores, ntop_scores


def _fmt_scores(scores: List) -> str:
    """Format a list of scores (float or None) as a comma-separated string, None -> NA."""
    return ', '.join(f"{s:.3f}" if s is not None else "NA" for s in scores)


def _parse_scores(cell) -> List:
    """Parse a comma-separated score cell back into a list of floats/None ('NA' -> None)."""
    if pd.isna(cell) or str(cell).strip() == '':
        return []
    return [None if s.strip() == 'NA' else float(s.strip()) for s in str(cell).split(',')]

def interactive_selection(cluster_num: str, normal_candidates: List[Dict], ntop_candidates: List[Dict], args) -> Tuple[str, List[str], List[str], List, List]:
    """Interactive selection of cluster label with separate normal (whole) and ntop columns.

    The whole/ntop distinction of which option is clicked no longer affects the
    output: for each selected label, both the whole and ntop scores are looked up.

    Returns:
        Tuple of (selected_name, selected_types, source_files, whole_scores,
        ntop_scores). selected_types is the literal list of chosen target-type
        labels (a label may itself contain '/', e.g. a mixed 'Dm8/Dm11'), and
        whole_scores/ntop_scores are positionally parallel to it. selected_name is
        the '/'-joined display name. For "keep original" / no candidates, returns
        (original_name, [], [], [], []).
    """
    cluster_label = f"{args.source_prefix}{cluster_num}"

    def scores_for(types):
        return _scores_for_types(cluster_label, types, normal_candidates, ntop_candidates, args.quiet)

    print(f"\n--- Cluster {cluster_label} ---")

    # Prepare display lists (top 10 from each type)
    normal_display = normal_candidates[:10]
    ntop_display = ntop_candidates[:10]

    if not normal_display and not ntop_display:
        print("No candidates found.")
        return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], [], [], []

    # Find common prefix in filenames to shorten display
    all_files = [c['file'] for c in normal_display] + [c['file'] for c in ntop_display]
    common_prefix = find_common_prefix(all_files)

    # Display candidates in two columns
    header_left = "NORMAL FILES"
    header_right = "NTOP FILES"
    content_width = 35  # Width for "XX. CellType (score: X.XXX)"
    print(f"{header_left}\t\t\t\t\t{header_right}")
    print("-" * 50 + "\t" + "-" * 50)

    max_rows = max(len(normal_display), len(ntop_display))
    option_num = 1
    option_map = {}

    for i in range(max_rows):
        left_text = ""
        right_text = ""

        # Normal file entry
        if i < len(normal_display):
            candidate = normal_display[i]
            short_file = candidate['file'].replace(common_prefix, '')
            score_text = f"{candidate['score']:.3f}"
            colored_score = colorize_score(candidate['score'], score_text)
            left_content = f"{option_num:2d}. {candidate['target_type']} (score: {colored_score})"
            left_text = f"{left_content:<{content_width}}\t[{short_file}]"
            option_map[option_num] = ('normal', candidate)
            option_num += 1

        # Ntop file entry
        if i < len(ntop_display):
            candidate = ntop_display[i]
            short_file = candidate['file'].replace(common_prefix, '')
            score_text = f"{candidate['score']:.3f}"
            colored_score = colorize_score(candidate['score'], score_text)
            right_content = f"{option_num:2d}. {candidate['target_type']} (score: {colored_score})"
            right_text = f"{right_content:<{content_width}}\t[{short_file}]"
            option_map[option_num] = ('ntop', candidate)
            option_num += 1

        print(f"{left_text}\t{right_text}")

    # Add option to keep original name
    keep_original_num = option_num
    print(f"\n{keep_original_num:2d}. Keep original name: {args.source_prefix.rstrip('_')}_{cluster_num}")

    while True:
        try:
            choice = input(f"Select option(s) (1-{keep_original_num}), comma-separated for multiple, or 'auto' for automatic selection: ").strip().lower()

            if choice == 'auto':
                # Use automatic selection logic combining both lists
                all_candidates = normal_candidates + ntop_candidates
                best_candidates = [c for c in all_candidates if c['score'] >= args.min_threshold]
                if not best_candidates:
                    print(f"AUTO: No candidates above min threshold ({args.min_threshold}), keeping original name")
                    return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], [], [], []

                # Apply threshold logic
                good_candidates = [c for c in best_candidates if c['score'] >= args.threshold]
                if good_candidates:
                    top_score = good_candidates[0]['score']
                    final_candidates = [
                        c for c in good_candidates
                        if c['score'] >= max(args.threshold, top_score - args.top_threshold)
                    ]
                    # Remove duplicates while preserving order
                    seen_types = set()
                    unique_candidates = []
                    for c in final_candidates:
                        if c['target_type'] not in seen_types:
                            unique_candidates.append(c)
                            seen_types.add(c['target_type'])

                    selected_names = [c['target_type'] for c in unique_candidates]
                    selected_files = [c['file'] for c in unique_candidates]
                    whole_scores, ntop_scores = scores_for(selected_names)
                    result = '/'.join(selected_names)
                    print(f"AUTO: Selected {result} (whole: {[round(s, 3) if s is not None else None for s in whole_scores]}, "
                          f"ntop: {[round(s, 3) if s is not None else None for s in ntop_scores]})")
                    return result, selected_names, selected_files, whole_scores, ntop_scores
                else:
                    best_candidate = best_candidates[0]
                    whole_scores, ntop_scores = scores_for([best_candidate['target_type']])
                    print(f"AUTO: Selected {best_candidate['target_type']} (score: {best_candidate['score']:.3f}, above min but below main threshold)")
                    return best_candidate['target_type'], [best_candidate['target_type']], [best_candidate['file']], whole_scores, ntop_scores

            # Handle multiple selections (comma-separated)
            if ',' in choice:
                choice_nums = [int(x.strip()) for x in choice.split(',')]
                selected_candidates = []
                selected_files = []

                for choice_num in choice_nums:
                    if choice_num in option_map:
                        file_type, candidate = option_map[choice_num]
                        selected_candidates.append(candidate['target_type'])
                        selected_files.append(candidate['file'])
                    elif choice_num == keep_original_num:
                        return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], [], [], []
                    else:
                        print(f"Invalid selection: {choice_num}. Please enter numbers between 1 and {keep_original_num}")
                        break
                else:
                    # All selections were valid
                    if selected_candidates:
                        whole_scores, ntop_scores = scores_for(selected_candidates)
                        return '/'.join(selected_candidates), selected_candidates, selected_files, whole_scores, ntop_scores

            else:
                # Single selection
                choice_num = int(choice)
                if choice_num in option_map:
                    file_type, candidate = option_map[choice_num]
                    whole_scores, ntop_scores = scores_for([candidate['target_type']])
                    return candidate['target_type'], [candidate['target_type']], [candidate['file']], whole_scores, ntop_scores
                elif choice_num == keep_original_num:
                    return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], [], [], []
                else:
                    print(f"Please enter a number between 1 and {keep_original_num}")

        except EOFError:
            print("\nReached end of input. Keeping original name.")
            return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], [], [], []
        except ValueError:
            if not choice:
                continue
            print("Please enter valid number(s), comma-separated for multiple selections, or 'auto'")

def extract_species_from_prefix(prefix: str) -> str:
    """Extract species abbreviation from prefix (e.g., 'sa_' -> 'sa')"""
    return prefix.rstrip('_')

def build_all_results(cluster_mappings: Dict, original_name_clusters: Dict, args) -> List[Dict]:
    """Build the list of output-row dicts (new whole/ntop column layout).

    Column order: cluster, new_name, target_types, whole_alignment_score,
    ntop_alignment_score, source_dataset, num_alignments, status.
    """
    species_abbrev = extract_species_from_prefix(args.source_prefix)
    all_results = []

    all_cluster_nums = set(cluster_mappings.keys()) | set(original_name_clusters.keys())

    for cluster_num in all_cluster_nums:
        cluster_id = f"{args.source_prefix}{cluster_num}"
        if cluster_num in cluster_mappings:
            mapping = cluster_mappings[cluster_num]
            all_results.append({
                'cluster': cluster_id,
                'new_name': mapping['name'],
                'target_types': ', '.join(mapping['types']),
                'whole_alignment_score': _fmt_scores(mapping['whole_scores']),
                'ntop_alignment_score': _fmt_scores(mapping['ntop_scores']),
                'source_dataset': ', '.join(paper_from_source(f) for f in mapping.get('source_files', [])),
                'num_alignments': len(mapping['types']),
                'status': 'assigned',
            })
        elif cluster_num in original_name_clusters:
            cluster_info = original_name_clusters[cluster_num]
            reason = cluster_info['reason']
            original_name = f"{species_abbrev}_{cluster_num}"

            if reason == 'below_min_threshold':
                all_results.append({
                    'cluster': cluster_id,
                    'new_name': original_name,
                    'target_types': cluster_info['best_target'],
                    'whole_alignment_score': _fmt_scores([cluster_info.get('whole')]),
                    'ntop_alignment_score': _fmt_scores([cluster_info.get('ntop')]),
                    'source_dataset': '',
                    'num_alignments': 1,
                    'status': 'below_min_threshold',
                })
            elif reason == 'no_candidates':
                all_results.append({
                    'cluster': cluster_id,
                    'new_name': original_name,
                    'target_types': 'No alignments',
                    'whole_alignment_score': 'NA',
                    'ntop_alignment_score': 'NA',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'no_candidates',
                })
            elif reason == 'user_selected':
                all_results.append({
                    'cluster': cluster_id,
                    'new_name': original_name,
                    'target_types': 'User selected original',
                    'whole_alignment_score': 'NA',
                    'ntop_alignment_score': 'NA',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'user_selected_original',
                })

    return all_results


def save_results_csv(all_results: List[Dict], csv_output_path: Path, args):
    """Sort rows by numeric cluster number and write to CSV."""
    if not all_results:
        return
    results_df = pd.DataFrame(all_results)

    def extract_cluster_num(cluster_name):
        num_str = cluster_name.replace(args.source_prefix, '')
        try:
            return int(num_str)
        except ValueError:
            return float('inf')

    results_df['sort_key'] = results_df['cluster'].apply(extract_cluster_num)
    results_df = results_df.sort_values('sort_key').drop('sort_key', axis=1)
    results_df.to_csv(csv_output_path, index=False)


def write_progress_to_csv(cluster_mappings: Dict, original_name_clusters: Dict, csv_output_path: Path, args):
    """Write current progress to CSV file (used for incremental saves in interactive mode)."""
    save_results_csv(build_all_results(cluster_mappings, original_name_clusters, args), csv_output_path, args)

def find_unmatched_targets(dataframes: Dict[str, pd.DataFrame], args, threshold: float = 0.3) -> List[Dict]:
    """
    Find target cell types that don't have alignment scores > threshold across all files.

    Args:
        dataframes: Dictionary of dataframes with alignment data
        args: Command line arguments containing source/target prefixes
        threshold: Minimum alignment score threshold (default: 0.3)

    Returns:
        List of dictionaries containing unmatched target info
    """
    # Track maximum alignment score for each target cell type
    target_max_scores = {}

    for file_name, df in dataframes.items():
        # Get target columns (fly cell types)
        target_columns = [col for col in df.columns if col.startswith(args.target_prefix)]

        for target_col in target_columns:
            # Extract the cell type name without prefix
            target_type = target_col.replace(args.target_prefix, '')

            # Get all alignment scores for this target across all source clusters
            scores = df[target_col].dropna()

            if len(scores) > 0:
                max_score_in_file = scores.max()

                # Update global maximum for this target type
                if target_type not in target_max_scores:
                    target_max_scores[target_type] = {
                        'max_score': max_score_in_file,
                        'file': file_name,
                        'source_cluster': df[target_col].idxmax()
                    }
                elif max_score_in_file > target_max_scores[target_type]['max_score']:
                    target_max_scores[target_type] = {
                        'max_score': max_score_in_file,
                        'file': file_name,
                        'source_cluster': df[target_col].idxmax()
                    }

    # Find targets with max score below threshold
    unmatched_targets = []
    for target_type, info in target_max_scores.items():
        if info['max_score'] < threshold:
            unmatched_targets.append({
                'target_type': target_type,
                'max_score': info['max_score'],
                'best_file': info['file'],
                'best_source_cluster': info['source_cluster']
            })

    # Sort by max score (ascending) to show worst matches first
    unmatched_targets.sort(key=lambda x: x['max_score'])

    return unmatched_targets

def main():
    parser = argparse.ArgumentParser(
        description='Process alignment table to map clusters to cell types',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input_table.csv -o output_dir
  %(prog)s file1.csv file2.csv --source-prefix sa_ --target-prefix fl_ --threshold 0.7
  %(prog)s alignment1.csv alignment2.csv --interactive --min-threshold 0.4
  %(prog)s *.csv --output-dir results --csv-output mappings.csv
        """
    )
    
    # Required arguments
    parser.add_argument('input_files', nargs='+', help='Input CSV files with alignment data (one or more files)')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', default='.', 
                       help='Output directory (default: current directory)')
    parser.add_argument('--csv-output', default='cluster_mappings.csv',
                       help='Output CSV filename (default: cluster_mappings.csv)')
    parser.add_argument('--source-prefix', default='sa_',
                       help='Prefix for source clusters (default: sa_)')
    parser.add_argument('--target-prefix', default='fl_',
                       help='Prefix for target cell types (default: fl_)')
    parser.add_argument('--threshold', type=float, default=0.6,
                       help='Minimum alignment score threshold (default: 0.6)')
    parser.add_argument('--min-threshold', type=float, default=0.3,
                       help='Absolute minimum threshold - clusters below this keep original names (default: 0.3)')
    parser.add_argument('--top-threshold', type=float, default=0.1,
                       help='Score difference from top score to include (default: 0.1)')
    parser.add_argument('--interactive', '-i', action='store_true',
                       help='Interactive mode: manually select labels for each cluster')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress detailed output')
    parser.add_argument('--find-unmatched-targets', action='store_true',
                       help='Identify target cell types with no alignment score > 0.3 across all files')

    args = parser.parse_args()

    # Validate thresholds
    if args.min_threshold > args.threshold:
        print(f"Error: min-threshold ({args.min_threshold}) cannot be higher than threshold ({args.threshold})", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load alignment data from all files
    dataframes = load_alignment_data(args.input_files, args)

    if not args.quiet:
        print(f"\nLoaded {len(dataframes)} alignment files")
        for name, df in dataframes.items():
            source_clusters = len([idx for idx in df.index if str(idx).startswith(args.source_prefix)])
            target_types = len([col for col in df.columns if col.startswith(args.target_prefix)])
            print(f"  {name}: {source_clusters} {args.source_prefix.rstrip('_')} clusters, {target_types} {args.target_prefix.rstrip('_')} cell types")

    # Extract mappings for each source cluster
    cluster_mappings = {}
    unassigned_clusters = {}
    original_name_clusters = {}

    # Load existing progress in interactive mode
    csv_output_path = output_dir / args.csv_output
    already_processed = set()

    if args.interactive and csv_output_path.exists():
        try:
            existing_df = pd.read_csv(csv_output_path)
            already_processed = set([row['cluster'] for _, row in existing_df.iterrows()])

            # CRITICAL: Load the actual data from existing rows into our dictionaries
            species_abbrev = extract_species_from_prefix(args.source_prefix)
            for _, row in existing_df.iterrows():
                cluster_num = row['cluster'].replace(args.source_prefix, '')

                if row['status'] == 'assigned':
                    # Parse the data back into cluster_mappings (new whole/ntop layout)
                    target_types = [t.strip() for t in str(row['target_types']).split(',')]
                    whole_scores = _parse_scores(row['whole_alignment_score'])
                    ntop_scores = _parse_scores(row['ntop_alignment_score'])
                    source_files = [f.strip() for f in str(row['source_dataset']).split(',')] if pd.notna(row['source_dataset']) and row['source_dataset'] else []

                    cluster_mappings[cluster_num] = {
                        'name': row['new_name'],
                        'types': target_types,
                        'whole_scores': whole_scores,
                        'ntop_scores': ntop_scores,
                        'source_files': source_files,
                    }
                elif row['status'] in ['below_min_threshold', 'no_candidates', 'user_selected_original']:
                    # Reconstruct original_name_clusters
                    reason_map = {
                        'below_min_threshold': 'below_min_threshold',
                        'no_candidates': 'no_candidates',
                        'user_selected_original': 'user_selected'
                    }
                    cluster_info = {'reason': reason_map[row['status']]}

                    if row['status'] == 'below_min_threshold':
                        cluster_info['best_target'] = row['target_types']
                        whole = _parse_scores(row['whole_alignment_score'])
                        ntop = _parse_scores(row['ntop_alignment_score'])
                        cluster_info['whole'] = whole[0] if whole else None
                        cluster_info['ntop'] = ntop[0] if ntop else None

                    original_name_clusters[cluster_num] = cluster_info

            if not args.quiet:
                print(f"\nFound existing progress file with {len(already_processed)} clusters already processed")
                print("Will resume from where you left off...\n")
        except Exception as e:
            if not args.quiet:
                print(f"Could not load existing progress: {e}")
                print("Starting fresh...\n")
            # Reset to empty if loading failed
            cluster_mappings = {}
            original_name_clusters = {}
            already_processed = set()

    # Get all unique source clusters from all files
    all_source_clusters = set()
    for df in dataframes.values():
        source_rows = [idx for idx in df.index if str(idx).startswith(args.source_prefix)]
        all_source_clusters.update(source_rows)

    if not all_source_clusters:
        print(f"Error: No rows found with prefix '{args.source_prefix}' in any file", file=sys.stderr)
        sys.exit(1)

    # Validate that at least one file has target columns
    has_targets = False
    for df in dataframes.values():
        target_columns = [col for col in df.columns if col.startswith(args.target_prefix)]
        if target_columns:
            has_targets = True
            break

    if not has_targets:
        print(f"Error: No columns found with prefix '{args.target_prefix}' in any file", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"\nProcessing {len(all_source_clusters)} {args.source_prefix.rstrip('_')} clusters...")

    species_abbrev = extract_species_from_prefix(args.source_prefix)

    for source_cluster in sorted(all_source_clusters, key=lambda x: int(str(x).replace(args.source_prefix, '')) if str(x).replace(args.source_prefix, '').isdigit() else float('inf')):
        cluster_num = str(source_cluster).replace(args.source_prefix, '')

        # Skip if already processed (in interactive mode with existing CSV)
        if source_cluster in already_processed:
            if not args.quiet:
                print(f"Skipping {source_cluster} (already processed)")
            continue

        # Get candidates from all files
        normal_candidates, ntop_candidates = get_cluster_candidates(dataframes, source_cluster, args)

        if args.interactive:
            # Interactive mode
            selected_name, selected_types, selected_files, whole_scores, ntop_scores = interactive_selection(cluster_num, normal_candidates, ntop_candidates, args)
            if selected_name == f"{species_abbrev}_{cluster_num}":
                original_name_clusters[cluster_num] = {'reason': 'user_selected'}
            else:
                cluster_mappings[cluster_num] = {
                    'name': selected_name,
                    'types': selected_types,
                    'whole_scores': whole_scores,
                    'ntop_scores': ntop_scores,
                    'source_files': selected_files if selected_files else [],
                }

            # Write progress to CSV after each selection
            write_progress_to_csv(cluster_mappings, original_name_clusters, csv_output_path, args)
            if not args.quiet:
                print(f"Progress saved to {csv_output_path}")
                print(f"Completed: {len(cluster_mappings) + len(original_name_clusters)} / {len(all_source_clusters)} clusters\n")
        else:
            # Automatic mode with new threshold logic
            all_candidates = normal_candidates + ntop_candidates
            all_candidates.sort(key=lambda x: x['score'], reverse=True)

            if not all_candidates:
                original_name_clusters[cluster_num] = {'reason': 'no_candidates'}
                if not args.quiet:
                    print(f"{args.source_prefix}{cluster_num} -> {species_abbrev}_{cluster_num} (no candidates)")
                continue

            # Filter by minimum threshold first
            above_min_threshold = [c for c in all_candidates if c['score'] >= args.min_threshold]

            if not above_min_threshold:
                # Below minimum threshold - keep original name
                best_target = all_candidates[0]['target_type']
                best_whole, best_ntop = get_label_scores(
                    f"{args.source_prefix}{cluster_num}", best_target,
                    normal_candidates, ntop_candidates, args.quiet)
                original_name_clusters[cluster_num] = {
                    'reason': 'below_min_threshold',
                    'best_score': all_candidates[0]['score'],
                    'best_target': best_target,
                    'whole': best_whole,
                    'ntop': best_ntop,
                }
                if not args.quiet:
                    print(f"{args.source_prefix}{cluster_num} -> {species_abbrev}_{cluster_num} (best: {all_candidates[0]['target_type']} score: {all_candidates[0]['score']:.3f}, below min threshold {args.min_threshold})")
                continue

            # Apply threshold logic to above-minimum candidates
            good_candidates = [c for c in above_min_threshold if c['score'] >= args.threshold]

            if good_candidates:
                # Sort by score and apply top-threshold logic
                good_candidates.sort(key=lambda x: x['score'], reverse=True)
                top_score = good_candidates[0]['score']
                final_candidates = [
                    c for c in good_candidates
                    if c['score'] >= max(args.threshold, top_score - args.top_threshold)
                ]

                # Remove duplicates while preserving order
                seen_types = set()
                unique_candidates = []
                for c in final_candidates:
                    if c['target_type'] not in seen_types:
                        unique_candidates.append(c)
                        seen_types.add(c['target_type'])

                target_types = [c['target_type'] for c in unique_candidates]
                source_files = [c['file'] for c in unique_candidates]
                whole_scores, ntop_scores = _scores_for_types(
                    f"{args.source_prefix}{cluster_num}", target_types,
                    normal_candidates, ntop_candidates, args.quiet)

                cluster_name = '/'.join(target_types)
                cluster_mappings[cluster_num] = {
                    'name': cluster_name,
                    'types': target_types,
                    'whole_scores': whole_scores,
                    'ntop_scores': ntop_scores,
                    'source_files': source_files
                }

                if not args.quiet:
                    print(f"{args.source_prefix}{cluster_num} -> {cluster_name} "
                          f"(whole: {[round(s, 3) if s is not None else None for s in whole_scores]}, "
                          f"ntop: {[round(s, 3) if s is not None else None for s in ntop_scores]})")
            else:
                # Above min threshold but below main threshold - assign best candidate
                best_candidate = above_min_threshold[0]
                whole_scores, ntop_scores = _scores_for_types(
                    f"{args.source_prefix}{cluster_num}", [best_candidate['target_type']],
                    normal_candidates, ntop_candidates, args.quiet)
                cluster_mappings[cluster_num] = {
                    'name': best_candidate['target_type'],
                    'types': [best_candidate['target_type']],
                    'whole_scores': whole_scores,
                    'ntop_scores': ntop_scores,
                    'source_files': [best_candidate['file']]
                }

                if not args.quiet:
                    print(f"{args.source_prefix}{cluster_num} -> {best_candidate['target_type']} (score: {best_candidate['score']:.3f}, below threshold {args.threshold})")

    # Print final results
    if not args.quiet:
        print("\n" + "="*60)
        print("ASSIGNED CLUSTER MAPPINGS")
        print("="*60)

        for cluster_num in sorted(cluster_mappings.keys(), key=lambda x: (0, int(x)) if x.isdigit() else (1, x)):
            mapping = cluster_mappings[cluster_num]
            print(f"{args.source_prefix}{cluster_num} -> {mapping['name']} "
                  f"(whole: [{_fmt_scores(mapping['whole_scores'])}], ntop: [{_fmt_scores(mapping['ntop_scores'])}])")

        print("\n" + "="*60)
        print("CLUSTERS KEEPING ORIGINAL NAMES")
        print("="*60)

        for cluster_num in sorted(original_name_clusters.keys(), key=lambda x: (0, int(x)) if x.isdigit() else (1, x)):
            cluster_info = original_name_clusters[cluster_num]
            reason = cluster_info['reason']
            if reason == 'below_min_threshold':
                print(f"{args.source_prefix}{cluster_num} -> {species_abbrev}_{cluster_num} (best: {cluster_info['best_target']} score: {cluster_info['best_score']:.3f}, below min threshold)")
            elif reason == 'no_candidates':
                print(f"{args.source_prefix}{cluster_num} -> {species_abbrev}_{cluster_num} (no candidates)")
            elif reason == 'user_selected':
                print(f"{args.source_prefix}{cluster_num} -> {species_abbrev}_{cluster_num} (user selected)")

        # Generate R code
        print("\n" + "="*60)
        print("R CODE FOR CLUSTER RENAMING")
        print("="*60)

        if cluster_mappings or original_name_clusters:
            print("# Add this code to your analysis script after clustering")
            print("# Replace the numbered cluster identities with cell type names")
            print()
            print("cluster_names <- c(")

            # Combine all cluster numbers and sort them numerically
            all_cluster_nums = set()
            all_cluster_nums.update(cluster_mappings.keys())
            all_cluster_nums.update(original_name_clusters.keys())

            # Sort numerically
            sorted_clusters = sorted(all_cluster_nums, key=lambda x: int(x) if x.isdigit() else float('inf'))

            # Build list of cluster assignments
            cluster_assignments = []
            for cluster_num in sorted_clusters:
                if cluster_num in cluster_mappings:
                    mapping = cluster_mappings[cluster_num]
                    if mapping["name"]:  # Only include non-empty names
                        cluster_assignments.append(f'  "{cluster_num}" = "{mapping["name"]}"')
                elif cluster_num in original_name_clusters:
                    cluster_assignments.append(f'  "{cluster_num}" = "{species_abbrev}_{cluster_num}"')

            # Print assignments with commas except for the last one
            for i, assignment in enumerate(cluster_assignments):
                if i < len(cluster_assignments) - 1:
                    print(assignment + ",")
                else:
                    print(assignment)

            print(")")
            print()
            print("# Apply the renaming")
            print("so2 <- RenameIdents(so2, cluster_names)")
            print()
            print("# Update plots with new names")
            print("DimPlot(so2, reduction = \"tsne\", label=T, label.size = 2, label.box = T)")
            print("DimPlot(so2, reduction = \"umap\", label=T, label.size = 2)")

    # Save results - include all clusters
    csv_output_path = output_dir / args.csv_output
    all_results = build_all_results(cluster_mappings, original_name_clusters, args)

    if all_results:
        save_results_csv(all_results, csv_output_path, args)
        print(f"\nResults saved to {csv_output_path}")

        assigned_count = len([r for r in all_results if r['status'] == 'assigned'])
        original_count = len([r for r in all_results if r['status'] in ['below_min_threshold', 'no_candidates', 'user_selected_original']])
        print(f"Found {assigned_count} clusters with alignments above threshold")
        print(f"Found {original_count} clusters keeping original names ({species_abbrev}_N format)")

        if args.interactive:
            print("Interactive mode was used for manual selection")
    else:
        print(f"\nNo clusters found to process")

    # Find unmatched target cell types if requested
    if args.find_unmatched_targets:
        if not args.quiet:
            print("\n" + "="*60)
            print("ANALYZING UNMATCHED TARGET CELL TYPES")
            print("="*60)
            print(f"Finding {args.target_prefix.rstrip('_')} cell types with max alignment score < 0.3...")

        unmatched_targets = find_unmatched_targets(dataframes, args, threshold=0.3)

        if unmatched_targets:
            # Save to file
            unmatched_output_path = output_dir / 'unmatched_targets.csv'
            unmatched_df = pd.DataFrame(unmatched_targets)
            unmatched_df.to_csv(unmatched_output_path, index=False)

            if not args.quiet:
                print(f"\nFound {len(unmatched_targets)} {args.target_prefix.rstrip('_')} cell types with no strong alignment (< 0.3):")
                print()
                for target_info in unmatched_targets:
                    print(f"  {args.target_prefix}{target_info['target_type']}")
                    print(f"    Best score: {target_info['max_score']:.3f}")
                    print(f"    Best match: {target_info['best_source_cluster']} (in {target_info['best_file']})")
                    print()

                print(f"\nUnmatched targets saved to {unmatched_output_path}")
        else:
            if not args.quiet:
                print(f"\nAll {args.target_prefix.rstrip('_')} cell types have at least one alignment score >= 0.3")

if __name__ == "__main__":
    main()
