#!/usr/bin/env python3
"""
Script to recursively find and delete all .log, .sbatch, and .batch files in a specified directory
Usage: python clean_logs.py <folder_path> [--dry-run]
"""

import os
import sys
import argparse
import glob
from pathlib import Path

def find_log_files(folder_path):
    """Find all .log, .sbatch, and .batch files recursively in the given folder"""
    log_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.log') or file.endswith('.sbatch') or file.endswith('.batch'):
                log_files.append(os.path.join(root, file))
    return log_files

def delete_log_files(log_files, dry_run=False):
    """Delete the log files (or just print them if dry_run=True)"""
    deleted_count = 0
    total_size = 0
    error_count = 0
    
    for i, log_file in enumerate(log_files, 1):
        try:
            if os.path.exists(log_file):
                file_size = os.path.getsize(log_file)
                total_size += file_size
                
                if dry_run:
                    if i <= 10:  # Only show first 10 in dry run
                        print(f"[DRY RUN] Would delete: {log_file} ({file_size} bytes)")
                else:
                    os.remove(log_file)
                    deleted_count += 1
                    if i % 100 == 0:  # Show progress every 100 files
                        print(f"Progress: {i}/{len(log_files)} files processed...")
        except Exception as e:
            error_count += 1
            if error_count <= 5:  # Only show first 5 errors
                print(f"Error deleting {log_file}: {e}")
    
    if error_count > 5:
        print(f"... and {error_count - 5} more errors")
    
    return deleted_count, total_size

def main():
    parser = argparse.ArgumentParser(description='Recursively find and delete all .log files')
    parser.add_argument('folder_path', help='Path to the folder to clean')
    parser.add_argument('--dry-run', action='store_true', 
                       help='Show what would be deleted without actually deleting')
    parser.add_argument('--confirm', action='store_true',
                       help='Skip confirmation prompt')
    
    args = parser.parse_args()
    
    # Check if folder exists
    if not os.path.exists(args.folder_path):
        print(f"Error: Folder '{args.folder_path}' does not exist")
        sys.exit(1)
    
    if not os.path.isdir(args.folder_path):
        print(f"Error: '{args.folder_path}' is not a directory")
        sys.exit(1)
    
    # Find all log, sbatch, and batch files
    print(f"Searching for .log, .sbatch, and .batch files in: {args.folder_path}")
    log_files = find_log_files(args.folder_path)
    
    if not log_files:
        print("No .log, .sbatch, or .batch files found")
        return
    
    print(f"Found {len(log_files)} .log, .sbatch, and .batch files")
    
    # Show first 10 files as examples
    if len(log_files) <= 10:
        for log_file in log_files:
            print(f"  {log_file}")
    else:
        print("Sample files:")
        for log_file in log_files[:5]:
            print(f"  {log_file}")
        print(f"  ... and {len(log_files) - 5} more files")
    
    
    # Confirmation
    if not args.dry_run and not args.confirm:
        response = input(f"\nAre you sure you want to delete {len(log_files)} log, sbatch, and batch files? (y/N): ")
        if response.lower() not in ['y', 'yes']:
            print("Operation cancelled")
            return
    
    # Delete files
    if args.dry_run:
        print("\n[DRY RUN MODE - No files will be deleted]")
        delete_log_files(log_files, dry_run=True)
    else:
        print(f"\nDeleting {len(log_files)} log, sbatch, and batch files...")
        deleted_count, deleted_size = delete_log_files(log_files, dry_run=False)
        print(f"\nCompleted: {deleted_count}/{len(log_files)} files deleted")
        print(f"Freed up: {deleted_size:,} bytes ({deleted_size/1024/1024:.2f} MB)")

if __name__ == "__main__":
    main()
