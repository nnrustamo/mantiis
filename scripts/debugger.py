import re
from decimal import Decimal, getcontext
from pathlib import Path

# Set precision to at least 20 digits for safety
getcontext().prec = 20

# Regex to parse lines like: rank = 0, f[0][0] = 0.104690908478217
line_re = re.compile(r"rank\s*=\s*(\d+),\s*f\[(\d+)\]\[(\d+)\]\s*=\s*([0-9.eE+-]+)")

def parse_file(filename):
    data = {}
    with open(filename, "r") as f:
        for line in f:
            m = line_re.search(line)
            if m:
                rank = int(m.group(1))
                i, j = m.group(2), m.group(3)
                key = f"f[{i}][{j}]"
                val = Decimal(m.group(4))
                data.setdefault(rank, {})[key] = val
    return data

def compare_with_ground_truth(basefile, comparefiles):
    ground = parse_file(basefile)
    ground_rank0 = ground[0]  # ground truth is rank 0 only

    for f in comparefiles:
        test = parse_file(f)
        print(f"\nComparing {f} against {basefile}:")

        for rank, entries in test.items():
            for key, val in entries.items():
                if key in ground_rank0:
                    diff = val - ground_rank0[key]
                    if diff != 0:
                        print(f"Mismatch: rank {rank}, {key}: {val} vs {ground_rank0[key]} (diff={diff})")

if __name__ == "__main__":
    basefile = "1_procs.txt"
    comparefiles = sorted(Path(".").glob("*_procs.txt"))
    comparefiles = [str(f) for f in comparefiles if not str(f).startswith("1_")]

    compare_with_ground_truth(basefile, comparefiles)
