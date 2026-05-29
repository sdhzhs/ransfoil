#!/usr/bin/env python
"""
Compare two files with 'Scaling' markers.
Extracts 3 data lines after each 'Scaling' marker and compares items.
Usage: python comlibcall.py <file1> <file2> [tolerance]
"""

import sys
import re

def extract_scaling_data(filename):
    """
    Extract data from 3 lines following 'Scaling' marker.
    Returns list of data sets, where each set contains 3 lines of data after 'Scaling'.
    """
    scaling_blocks = []
    
    try:
        lines = open(filename, 'r').readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Check if this line contains 'Scaling' (case-sensitive)
            if 'Scaling' in line:
                # Extract data from the next 3 lines
                if i + 3 < len(lines):
                    block_data = []
                    for j in range(1, 4):  # Next 3 lines
                        data_line = lines[i + j].strip()
                        # Extract all numbers from the line (use non-capturing group for exponent)
                        numbers = re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?', data_line)
                        if numbers:
                            try:
                                data_values = [float(x) for x in numbers]
                                block_data.extend(data_values)
                            except ValueError:
                                pass
                    
                    if block_data:
                        scaling_blocks.append(block_data)
                    i += 3
            
            i += 1
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        sys.exit(1)
    
    return scaling_blocks


def compare_files(file1, file2, tolerance=1e-2):
    """
    Compare data from 3 lines after 'Scaling' markers between two files.
    """
    data1 = extract_scaling_data(file1)
    data2 = extract_scaling_data(file2)
    
    print(f"File 1: {file1}")
    print(f"File 2: {file2}")
    print(f"Tolerance: {tolerance}")
    print(f"Number of 'Scaling' blocks - File 1: {len(data1)}, File 2: {len(data2)}")
    print("-" * 70)
    
    if len(data1) != len(data2):
        print(f"Warning: Different number of Scaling blocks ({len(data1)} vs {len(data2)})")
    
    all_pass = True
    num_blocks = min(len(data1), len(data2))
    
    for block_idx in range(num_blocks):
        block1 = data1[block_idx]
        block2 = data2[block_idx]
        
        if len(block1) != len(block2):
            print(f"Scaling block {block_idx + 1}: Warning - different data counts ({len(block1)} vs {len(block2)})")
            all_pass = False
        
        num_values = min(len(block1), len(block2))
        block_pass = True
        
        for val_idx in range(num_values):
            val1 = block1[val_idx]
            val2 = block2[val_idx]
            
            # Calculate relative error
            if val2 != 0.0:
                abs_rel_err = abs((val1 - val2) / val2)
                rel_err = (val1 - val2) * 100.0 / val2
            else:
                abs_rel_err = abs(val1 - val2)
                rel_err = (val1 - val2) * 100.0
            
            if abs_rel_err > tolerance:
                if block_pass:
                    print(f"Scaling block {block_idx + 1}:")
                    block_pass = False
                    all_pass = False
                print(f"  Item {val_idx + 1}: File1={val1:.6E}, File2={val2:.6E}, RelErr={rel_err:.4f}%")
        
        if block_pass:
            print(f"Scaling block {block_idx + 1}: PASS (all items within tolerance)")
    
    print("-" * 70)
    if all_pass:
        print(f"PASS: All data items within tolerance {tolerance}")
        return 0
    else:
        print(f"FAIL: Some items exceed tolerance {tolerance}")
        return 1


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python comlibcall.py <file1> <file2> [tolerance]")
        print("  file1, file2: Files to compare")
        print("  tolerance: Relative error tolerance (default: 1e-2)")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    if len(sys.argv) > 3:
        try:
            tol = float(sys.argv[3])
        except ValueError:
            print(f"Error: Invalid tolerance value '{sys.argv[3]}'")
            sys.exit(1)
    else:
        tol = 1e-2
    
    result = compare_files(file1, file2, tol)
    sys.exit(result)
