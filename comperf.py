def readperfdata(filename):
    """Read iteration elapsed time from output file
    Format: 'Iteration elapsed time:           0 m          17 s'
    Returns: time in seconds (last occurrence)
    """
    last_time = None
    for line in open(filename):
        if 'Iteration elapsed time:' in line:
            # Extract minutes and seconds
            parts = [str for str in line.strip().split() if str != '']
            # Find 'm' and 's' markers
            minutes = 0
            seconds = 0
            for i, part in enumerate(parts):
                if part == 'm' and i > 0:
                    try:
                        minutes = int(parts[i-1])
                    except ValueError:
                        pass
                elif part == 's' and i > 0:
                    try:
                        seconds = int(parts[i-1])
                    except ValueError:
                        pass
            last_time = minutes * 60 + seconds
    return last_time

import sys

if len(sys.argv) < 3:
    print("Usage: python comperf.py <file1> <file2> [tolerance]")
    print("  file1: first output file path")
    print("  file2: second output file path (baseline/golden)")
    print("  tolerance: optional, acceptable performance difference (default: 0.1 = 10%)")
    sys.exit(1)

filename1 = sys.argv[1]
filename2 = sys.argv[2]

time1 = readperfdata(filename1)
time2 = readperfdata(filename2)

if time1 is None:
    print(f"Error: No 'Iteration elapsed time' found in {filename1}")
    sys.exit(1)

if time2 is None:
    print(f"Error: No 'Iteration elapsed time' found in {filename2}")
    sys.exit(1)

# Set tolerance (default 10% performance difference allowed)
if len(sys.argv) > 3:
    tol = float(sys.argv[3])
else:
    tol = 0.1

print(f"\nPerformance Comparison:")
print(f"{'Metric':<25} {'CurrentRun':<15} {'Golden':<15} {'Speedup':<12} {'Diff %':<12}")
print("-" * 79)

# Calculate speedup and difference
if time1 != 0.0:
    speedup = time2 / time1
else:
    speedup = 0.0

if time2 != 0.0:
    diff_percent = (time1 - time2) * 100 / time2
else:
    diff_percent = 0.0

# Format times as m:ss
def format_time(seconds):
    mins = seconds // 60
    secs = seconds % 60
    return f"{mins}m {secs}s"

time1_str = format_time(time1)
time2_str = format_time(time2)

# Determine pass/fail status
# Fail only if time1 is more than tol (10%) slower than golden (time2)
# Pass if time1 is faster or within tolerance
status = ""
flag = True
if time1 > time2 * (1 + tol):
    flag = False
    status = " * SLOW"

print(f"{'Iteration elapsed time':<25} {time1_str:<15} {time2_str:<15} {speedup:<12.4f} {diff_percent:<12.2f}{status}")
print("-" * 79)

# Check for significant performance improvement (10% or more faster)
if time2 != 0.0 and time1 < time2 * (1 - tol):
    improvement = (time2 - time1) * 100 / time2
    print(f"\n[IMPROVED] Performance improved by {improvement:.2f}% (more than {tol*100}% faster than golden)")

if flag:
    if time1 <= time2:
        print(f"\n[PASS] Performance is equal to or better than golden")
    else:
        print(f"\n[PASS] Performance is within acceptable range (not more than {tol*100}% slower than golden)")
    sys.exit(0)
else:
    print(f"\n[FAIL] Performance is too slow (more than {tol*100}% slower than golden)")
    sys.exit(1)
