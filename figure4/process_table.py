# Version 1.0, 99% CERTIFIED CLAUDE AI GENERATED
# A helper script that displays the best scores across multiple SAMap alignment score tables
# and creates an interactive CLI that allows user to select chosen mappings per cluster
# also outputs chosen results and other info to a spreadsheet
import pandas as pd
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple


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
        top_candidates = candidates[:10]

        # Categorize by file type (accept 'top', 'ntop', etc.)
        if 'top' in file_name.lower():
            ntop_candidates.extend(top_candidates)
        else:
            normal_candidates.extend(top_candidates)

    # Sort both lists by score
    normal_candidates.sort(key=lambda x: x['score'], reverse=True)
    ntop_candidates.sort(key=lambda x: x['score'], reverse=True)

    return normal_candidates, ntop_candidates

def _build_ntop_info(selected_types_info: List[Tuple[str, Dict]], normal_candidates: List[Dict], ntop_candidates: List[Dict]) -> Dict:
    """Build NTOP tracking information for selected candidates

    Args:
        selected_types_info: List of tuples (file_type, candidate_dict) for each selected type
        normal_candidates: List of candidates from normal files
        ntop_candidates: List of candidates from NTOP files

    Returns:
        Dictionary with lists of NTOP tracking info for each selected type:
        {
            'ntop_chosen_list': ['yes', 'no', 'yes'],
            'alternative_file_type_list': ['WHOLE', 'NTOP', 'NA'],
            'alternative_score_list': ['0.878', '0.789', 'NA']
        }
    """
    if not selected_types_info:
        return {}

    ntop_chosen_list = []
    alternative_file_type_list = []
    alternative_score_list = []

    # Process each selected type individually
    for file_type, selected_candidate in selected_types_info:
        selected_target_type = selected_candidate['target_type']
        ntop_chosen = (file_type == 'ntop')

        # Find alternative score from the opposite file type
        alternative_score = None
        if ntop_chosen:
            # Selected from NTOP, look for same target_type in normal files
            for cand in normal_candidates:
                if cand['target_type'] == selected_target_type:
                    alternative_score = cand['score']
                    break
        else:
            # Selected from normal, look for same target_type in NTOP files
            for cand in ntop_candidates:
                if cand['target_type'] == selected_target_type:
                    alternative_score = cand['score']
                    break

        # Add to lists
        ntop_chosen_list.append('yes' if ntop_chosen else 'no')

        if alternative_score is not None:
            alternative_file_type_list.append('WHOLE' if ntop_chosen else 'NTOP')
            alternative_score_list.append(f"{alternative_score:.3f}")
        else:
            alternative_file_type_list.append('NA')
            alternative_score_list.append('NA')

    return {
        'ntop_chosen_list': ntop_chosen_list,
        'alternative_file_type_list': alternative_file_type_list,
        'alternative_score_list': alternative_score_list
    }

def interactive_selection(cluster_num: str, normal_candidates: List[Dict], ntop_candidates: List[Dict], args) -> Tuple[str, List[str], Dict]:
    """Interactive selection of cluster label with separate normal and ntop columns

    Returns:
        Tuple of (selected_name, list_of_source_files, ntop_info_dict)
        ntop_info_dict contains: {
            'ntop_chosen': bool,
            'selected_type': str,
            'selected_score': float,
            'alternative_score': float or None,
            'alternative_type': str ('normal' or 'ntop')
        }
    """
    print(f"\n--- Cluster {args.source_prefix}{cluster_num} ---")

    # Prepare display lists (top 10 from each type)
    normal_display = normal_candidates[:10]
    ntop_display = ntop_candidates[:10]

    if not normal_display and not ntop_display:
        print("No candidates found.")
        return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], {}

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
                    return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], {}

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
                    selected_scores = [c['score'] for c in unique_candidates]
                    selected_files = [c['file'] for c in unique_candidates]
                    result = '/'.join(selected_names)
                    print(f"AUTO: Selected {result} (scores: {[round(s, 3) for s in selected_scores]})")
                    return result, selected_files, {}
                else:
                    best_candidate = best_candidates[0]
                    print(f"AUTO: Selected {best_candidate['target_type']} (score: {best_candidate['score']:.3f}, above min but below main threshold)")
                    return best_candidate['target_type'], [best_candidate['file']], {}

            # Handle multiple selections (comma-separated)
            if ',' in choice:
                choice_nums = [int(x.strip()) for x in choice.split(',')]
                selected_candidates = []
                selected_files = []
                selected_types_info = []

                for choice_num in choice_nums:
                    if choice_num in option_map:
                        file_type, candidate = option_map[choice_num]
                        selected_candidates.append(candidate['target_type'])
                        selected_files.append(candidate['file'])
                        selected_types_info.append((file_type, candidate))
                    elif choice_num == keep_original_num:
                        return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], {}
                    else:
                        print(f"Invalid selection: {choice_num}. Please enter numbers between 1 and {keep_original_num}")
                        break
                else:
                    # All selections were valid
                    if selected_candidates:
                        # For multiple selections, use the first selection to determine ntop_info
                        ntop_info = _build_ntop_info(selected_types_info, normal_candidates, ntop_candidates)
                        # Also extract scores from selected_types_info
                        selected_scores = [cand['score'] for _, cand in selected_types_info]
                        ntop_info['selected_scores'] = selected_scores
                        return '/'.join(selected_candidates), selected_files, ntop_info

            else:
                # Single selection
                choice_num = int(choice)
                if choice_num in option_map:
                    file_type, candidate = option_map[choice_num]
                    ntop_info = _build_ntop_info([(file_type, candidate)], normal_candidates, ntop_candidates)
                    # Include the selected score
                    ntop_info['selected_scores'] = [candidate['score']]
                    return candidate['target_type'], [candidate['file']], ntop_info
                elif choice_num == keep_original_num:
                    return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], {}
                else:
                    print(f"Please enter a number between 1 and {keep_original_num}")

        except EOFError:
            print("\nReached end of input. Keeping original name.")
            return f"{args.source_prefix.rstrip('_')}_{cluster_num}", [], {}
        except ValueError:
            if not choice:
                continue
            print("Please enter valid number(s), comma-separated for multiple selections, or 'auto'")

def extract_species_from_prefix(prefix: str) -> str:
    """Extract species abbreviation from prefix (e.g., 'sa_' -> 'sa')"""
    return prefix.rstrip('_')

def write_progress_to_csv(cluster_mappings: Dict, original_name_clusters: Dict, csv_output_path: Path, args):
    """Write current progress to CSV file

    Args:
        cluster_mappings: Dictionary of assigned cluster mappings
        original_name_clusters: Dictionary of clusters keeping original names
        csv_output_path: Path to output CSV file
        args: Command line arguments
    """
    species_abbrev = extract_species_from_prefix(args.source_prefix)
    all_results = []

    # Combine all cluster numbers
    all_cluster_nums = set()
    all_cluster_nums.update(cluster_mappings.keys())
    all_cluster_nums.update(original_name_clusters.keys())

    for cluster_num in all_cluster_nums:
        if cluster_num in cluster_mappings:
            mapping = cluster_mappings[cluster_num]
            ntop_info = mapping.get('ntop_info', {})

            # Format NTOP tracking columns
            if ntop_info:
                ntop_chosen = ', '.join(ntop_info.get('ntop_chosen_list', []))
                alternative_file_type = ', '.join(ntop_info.get('alternative_file_type_list', []))
                alternative_score = ', '.join(ntop_info.get('alternative_score_list', []))
            else:
                ntop_chosen = 'N/A'
                alternative_file_type = 'N/A'
                alternative_score = 'N/A'

            all_results.append({
                'cluster': f"{args.source_prefix}{cluster_num}",
                'new_name': mapping['name'],
                'target_types': ', '.join(mapping['types']),
                'alignment_scores': ', '.join([f"{s:.3f}" for s in mapping['scores']]),
                'source_dataset': ', '.join(mapping.get('source_files', [])),
                'num_alignments': len(mapping['types']),
                'status': 'assigned',
                'ntop_chosen': ntop_chosen,
                'alternative_file_type': alternative_file_type,
                'alternative_score': alternative_score
            })
        elif cluster_num in original_name_clusters:
            cluster_info = original_name_clusters[cluster_num]
            reason = cluster_info['reason']
            original_name = f"{species_abbrev}_{cluster_num}"

            if reason == 'below_min_threshold':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': cluster_info['best_target'],
                    'alignment_scores': f"{cluster_info['best_score']:.3f}",
                    'source_dataset': '',
                    'num_alignments': 1,
                    'status': 'below_min_threshold',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })
            elif reason == 'no_candidates':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': 'No alignments',
                    'alignment_scores': '0.000',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'no_candidates',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })
            elif reason == 'user_selected':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': 'User selected original',
                    'alignment_scores': '',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'user_selected_original',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })

    if all_results:
        results_df = pd.DataFrame(all_results)
        # Sort by numeric cluster number for proper ordering
        def extract_cluster_num(cluster_name):
            num_str = cluster_name.replace(args.source_prefix, '')
            try:
                return int(num_str)
            except ValueError:
                return float('inf')

        results_df['sort_key'] = results_df['cluster'].apply(extract_cluster_num)
        results_df = results_df.sort_values('sort_key').drop('sort_key', axis=1)
        results_df.to_csv(csv_output_path, index=False)

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
                    # Parse the data back into cluster_mappings
                    target_types = [t.strip() for t in str(row['target_types']).split(',')]
                    scores_str = [s.strip() for s in str(row['alignment_scores']).split(',')]
                    scores = [float(s) for s in scores_str]
                    source_files = [f.strip() for f in str(row['source_dataset']).split(',')] if pd.notna(row['source_dataset']) and row['source_dataset'] else []

                    # Reconstruct ntop_info
                    ntop_info = {}
                    if pd.notna(row['ntop_chosen']) and row['ntop_chosen'] != 'N/A':
                        ntop_info['ntop_chosen_list'] = [v.strip() for v in str(row['ntop_chosen']).split(',')]
                        ntop_info['alternative_file_type_list'] = [v.strip() for v in str(row['alternative_file_type']).split(',')]
                        ntop_info['alternative_score_list'] = [v.strip() for v in str(row['alternative_score']).split(',')]
                        ntop_info['selected_scores'] = scores

                    cluster_mappings[cluster_num] = {
                        'name': row['new_name'],
                        'types': target_types,
                        'scores': scores,
                        'source_files': source_files,
                        'ntop_info': ntop_info
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
                        cluster_info['best_score'] = float(row['alignment_scores'])

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
            selected_name, selected_files, ntop_info = interactive_selection(cluster_num, normal_candidates, ntop_candidates, args)
            if selected_name == f"{species_abbrev}_{cluster_num}":
                original_name_clusters[cluster_num] = {'reason': 'user_selected'}
            else:
                # Extract types and scores
                selected_types = selected_name.split('/')

                # Use scores from ntop_info if available (from interactive selection)
                if ntop_info and 'selected_scores' in ntop_info:
                    selected_scores = ntop_info['selected_scores']
                else:
                    # Fallback: try to match scores (for backwards compatibility)
                    all_candidates = normal_candidates + ntop_candidates
                    selected_scores = []
                    for stype in selected_types:
                        for cand in all_candidates:
                            if cand['target_type'] == stype:
                                selected_scores.append(cand['score'])
                                break

                # Use selected_files if provided, otherwise empty
                selected_files_final = selected_files if selected_files else []

                cluster_mappings[cluster_num] = {
                    'name': selected_name,
                    'types': selected_types,
                    'scores': selected_scores if selected_scores else [0.0] * len(selected_types),
                    'source_files': selected_files_final,
                    'ntop_info': ntop_info
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
                original_name_clusters[cluster_num] = {
                    'reason': 'below_min_threshold',
                    'best_score': all_candidates[0]['score'],
                    'best_target': all_candidates[0]['target_type']
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
                scores = [c['score'] for c in unique_candidates]
                source_files = [c['file'] for c in unique_candidates]

                cluster_name = '/'.join(target_types)
                cluster_mappings[cluster_num] = {
                    'name': cluster_name,
                    'types': target_types,
                    'scores': scores,
                    'source_files': source_files
                }

                if not args.quiet:
                    print(f"{args.source_prefix}{cluster_num} -> {cluster_name} (scores: {[round(s, 3) for s in scores]})")
            else:
                # Above min threshold but below main threshold - assign best candidate
                best_candidate = above_min_threshold[0]
                cluster_mappings[cluster_num] = {
                    'name': best_candidate['target_type'],
                    'types': [best_candidate['target_type']],
                    'scores': [best_candidate['score']],
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
            score_str = ', '.join([f"{s:.3f}" for s in mapping['scores']])
            print(f"{args.source_prefix}{cluster_num} -> {mapping['name']} (scores: {score_str})")

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
    all_results = []

    # Create a combined list of all cluster numbers for proper sorting
    all_cluster_nums = set()
    all_cluster_nums.update(cluster_mappings.keys())
    all_cluster_nums.update(original_name_clusters.keys())

    for cluster_num in all_cluster_nums:
        if cluster_num in cluster_mappings:
            mapping = cluster_mappings[cluster_num]
            ntop_info = mapping.get('ntop_info', {})

            # Format NTOP tracking columns
            if ntop_info:
                ntop_chosen = ', '.join(ntop_info.get('ntop_chosen_list', []))
                alternative_file_type = ', '.join(ntop_info.get('alternative_file_type_list', []))
                alternative_score = ', '.join(ntop_info.get('alternative_score_list', []))
            else:
                # No NTOP info (e.g., auto mode)
                ntop_chosen = 'N/A'
                alternative_file_type = 'N/A'
                alternative_score = 'N/A'

            all_results.append({
                'cluster': f"{args.source_prefix}{cluster_num}",
                'new_name': mapping['name'],
                'target_types': ', '.join(mapping['types']),
                'alignment_scores': ', '.join([f"{s:.3f}" for s in mapping['scores']]),
                'source_dataset': ', '.join(mapping.get('source_files', [])),
                'num_alignments': len(mapping['types']),
                'status': 'assigned',
                'ntop_chosen': ntop_chosen,
                'alternative_file_type': alternative_file_type,
                'alternative_score': alternative_score
            })
        elif cluster_num in original_name_clusters:
            cluster_info = original_name_clusters[cluster_num]
            reason = cluster_info['reason']
            original_name = f"{species_abbrev}_{cluster_num}"

            if reason == 'below_min_threshold':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': cluster_info['best_target'],
                    'alignment_scores': f"{cluster_info['best_score']:.3f}",
                    'source_dataset': '',
                    'num_alignments': 1,
                    'status': 'below_min_threshold',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })
            elif reason == 'no_candidates':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': 'No alignments',
                    'alignment_scores': '0.000',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'no_candidates',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })
            elif reason == 'user_selected':
                all_results.append({
                    'cluster': f"{args.source_prefix}{cluster_num}",
                    'new_name': original_name,
                    'target_types': 'User selected original',
                    'alignment_scores': '',
                    'source_dataset': '',
                    'num_alignments': 0,
                    'status': 'user_selected_original',
                    'ntop_chosen': 'N/A',
                    'alternative_file_type': 'N/A',
                    'alternative_score': 'N/A'
                })

    if all_results:
        results_df = pd.DataFrame(all_results)
        # Sort by numeric cluster number for proper ordering
        def extract_cluster_num(cluster_name):
            num_str = cluster_name.replace(args.source_prefix, '')
            try:
                return int(num_str)
            except ValueError:
                return float('inf')  # Put non-numeric cluster names at the end

        results_df['sort_key'] = results_df['cluster'].apply(extract_cluster_num)
        results_df = results_df.sort_values('sort_key').drop('sort_key', axis=1)
        results_df.to_csv(csv_output_path, index=False)
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
