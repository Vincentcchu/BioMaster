#!/usr/bin/env python3
"""
Automate running BioMaster classifier on all h5ad files across tissue types.

This script:
1. Finds all h5ad files in data/*/h5ad_unlabelled/
2. For each file, creates a temporary config
3. Runs the classifier
4. Tracks progress and creates summaries
"""

import os
import sys
import yaml
import subprocess
from pathlib import Path
import json
from datetime import datetime

# Configuration
DATA_BASE_DIR = "/cs/student/projects2/aisd/2024/shekchu/projects/data"
CONFIG_TEMPLATE = "config.yaml"
RUN_SCRIPT = "run_tracked2.py"
GOAL = 'Annotate cell types for each cell, and store the predicted label in adata.obs["cell_type"]. No placeholders—supply fully validated, complete details for markers for Lymphoid, Myeloid, Endothelial, Fibroblasts, Epithelial_non_malignant, Malignant, Benign, Other_non_neoplastic. Save the updated AnnData object to a .h5ad file.'

def load_config_template(config_path):
    """Load the template config file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def find_all_h5ad_files(data_base_dir, h5ad_subdir="h5ad_unlabelled_clustered"):
    """
    Find all h5ad files organized by tissue type.
    
    Returns: list of tuples (tissue_name, h5ad_file_path, dataset_name)
    """
    all_files = []
    data_path = Path(data_base_dir)
    
    if not data_path.exists():
        print(f"Error: Data directory not found: {data_base_dir}")
        return []
    
    # Iterate through tissue folders
    for tissue_dir in sorted(data_path.iterdir()):
        if not tissue_dir.is_dir():
            continue
        
        # Skip non-tissue folders
        if tissue_dir.name in ['metadata_processing', '__pycache__', '.git', '.ipynb_checkpoints']:
            continue
        
        # Look for h5ad_unlabelled subdirectory
        h5ad_dir = tissue_dir / h5ad_subdir
        if not h5ad_dir.exists():
            continue
        
        # Find all h5ad files
        for h5ad_file in sorted(h5ad_dir.glob("*.h5ad")):
            dataset_name = h5ad_file.stem  # filename without extension
            all_files.append((tissue_dir.name, str(h5ad_file), dataset_name))
    
    return all_files

def create_run_id(tissue_name, dataset_name, task="cell_type_annotation"):
    """Generate a unique run ID"""
    # Clean up dataset name (remove "Data_" prefix if exists)
    clean_name = dataset_name.replace("Data_", "")
    
    # Remove tissue name suffix if present (case-insensitive)
    # e.g., "Choudhury2022_Brain" -> "Choudhury2022"
    tissue_suffix = f"_{tissue_name.capitalize()}"
    if clean_name.endswith(tissue_suffix):
        clean_name = clean_name[:-len(tissue_suffix)]
    
    # Remove underscores and convert to lowercase
    clean_name = clean_name.replace("_", "").lower()
    
    return f"{tissue_name}_{clean_name}_{task}_clustered"

def create_temp_config(template_config, tissue_name, h5ad_path, dataset_name, output_path):
    """Create a temporary config file for this specific run"""
    config = template_config.copy()
    
    # Update run ID
    run_id = create_run_id(tissue_name, dataset_name)
    config['biomaster']['id'] = run_id
    
    # Update file path with description
    file_description = f"{h5ad_path}: A unannotated single-cell RNA-seq {tissue_name} cancer dataset containing gene expression profiles with cell barcodes and gene names. The AnnData object includes precomputed clustering information stored in the `obs['cluster']` column."
    config['data']['files'] = [file_description]
    
    # Keep the goal from template or use default
    if 'goal' not in config['data']:
        config['data']['goal'] = GOAL
    
    # Write temporary config
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    return run_id

def run_classifier(config_path, run_id, dry_run=False):
    """Run the classifier script"""
    if dry_run:
        print(f"  [DRY-RUN] Would run: python {RUN_SCRIPT} {config_path}")
        return True, 0, "dry-run"
    
    try:
        print(f"  [Running] python {RUN_SCRIPT} {config_path}")
        result = subprocess.run(
            ['python', RUN_SCRIPT, config_path],
            capture_output=False,  # Show output in real-time
            text=True
        )
        return result.returncode == 0, result.returncode, "completed"
    except Exception as e:
        print(f"  [ERROR] Failed to run: {e}")
        return False, -1, str(e)

def save_batch_summary(results, output_file="./output/BATCH_SUMMARY.json"):
    """Save summary of all runs"""
    summary = {
        "batch_start": results[0]['start_time'] if results else None,
        "batch_end": datetime.now().isoformat(),
        "total_runs": len(results),
        "successful": sum(1 for r in results if r['success']),
        "failed": sum(1 for r in results if not r['success']),
        "runs": results
    }
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=4)
    
    return output_file

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Run classifier on all h5ad files across tissues")
    parser.add_argument('--data-dir', type=str, default=DATA_BASE_DIR,
                        help='Base data directory containing tissue folders')
    parser.add_argument('--config-template', type=str, default=CONFIG_TEMPLATE,
                        help='Template config file to use')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without actually running')
    parser.add_argument('--tissues', nargs='+', default=None,
                        help='Specific tissues to process (default: all)')
    parser.add_argument('--skip-existing', action='store_true',
                        help='Skip runs if output summary already exists')
    
    args = parser.parse_args()
    
    # Load template config
    if not os.path.exists(args.config_template):
        print(f"Error: Template config not found: {args.config_template}")
        sys.exit(1)
    
    template_config = load_config_template(args.config_template)
    
    # Find all h5ad files
    print(f"[INFO] Scanning for h5ad files in: {args.data_dir}")
    all_files = find_all_h5ad_files(args.data_dir)
    
    if not all_files:
        print("[WARN] No h5ad files found!")
        sys.exit(1)
    
    # Filter by tissue if specified
    if args.tissues:
        all_files = [(t, f, d) for t, f, d in all_files if t in args.tissues]
        print(f"[INFO] Filtering to tissues: {args.tissues}")
    
    print(f"[INFO] Found {len(all_files)} file(s) across {len(set(t for t, _, _ in all_files))} tissue(s)")
    print()
    
    # Group by tissue for display
    by_tissue = {}
    for tissue, filepath, dataset in all_files:
        by_tissue.setdefault(tissue, []).append(dataset)
    
    for tissue, datasets in sorted(by_tissue.items()):
        print(f"  {tissue}: {len(datasets)} file(s)")
    print()
    
    if args.dry_run:
        print("[INFO] DRY RUN MODE - No actual processing will occur\n")
    
    # Process each file
    results = []
    temp_config_path = "temp_config.yaml"
    
    for idx, (tissue, h5ad_path, dataset_name) in enumerate(all_files, 1):
        print(f"\n{'='*70}")
        print(f"Processing {idx}/{len(all_files)}: {tissue}/{dataset_name}")
        print(f"{'='*70}")
        
        run_id = create_run_id(tissue, dataset_name)
        
        # Check if already processed
        summary_file = f"./output/{run_id}_RUN_SUMMARY.json"
        if args.skip_existing and os.path.exists(summary_file):
            print(f"  [SKIP] Output already exists: {summary_file}")
            results.append({
                'tissue': tissue,
                'dataset': dataset_name,
                'run_id': run_id,
                'success': True,
                'status': 'skipped',
                'start_time': None
            })
            continue
        
        start_time = datetime.now()
        
        # Create temporary config
        create_temp_config(template_config, tissue, h5ad_path, dataset_name, temp_config_path)
        print(f"  Run ID: {run_id}")
        print(f"  File: {h5ad_path}")
        
        # Run classifier
        success, return_code, status = run_classifier(temp_config_path, run_id, args.dry_run)
        
        # Record result
        results.append({
            'tissue': tissue,
            'dataset': dataset_name,
            'run_id': run_id,
            'file_path': h5ad_path,
            'success': success,
            'return_code': return_code,
            'status': status,
            'start_time': start_time.isoformat()
        })
        
        if success:
            print(f"  ✓ Completed successfully")
        else:
            print(f"  ✗ Failed with code {return_code}")
    
    # Clean up temp config
    if os.path.exists(temp_config_path):
        os.remove(temp_config_path)
    
    # Save batch summary
    if not args.dry_run:
        summary_file = save_batch_summary(results)
        print(f"\n{'='*70}")
        print(f"BATCH SUMMARY")
        print(f"{'='*70}")
        print(f"Total runs:    {len(results)}")
        print(f"Successful:    {sum(1 for r in results if r['success'])}")
        print(f"Failed:        {sum(1 for r in results if not r['success'])}")
        print(f"Summary saved: {summary_file}")
        print(f"{'='*70}\n")

if __name__ == "__main__":
    main()


# Dry run - see what will be processed
# python run_all_tissues.py --dry-run

# Run all tissues
# python run_all_tissues.py

# Run specific tissues only
# python run_all_tissues.py --tissues breast brain colorectal

# Skip files that already have outputs
# python run_all_tissues.py --skip-existing