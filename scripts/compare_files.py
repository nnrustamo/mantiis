import sys

def compare_files(file1, file2, tol=1e-10):
    with open(file1, 'r') as f:
        nums1 = [float(line.strip()) for line in f if line.strip()]

    with open(file2, 'r') as f:
        nums2 = [float(line.strip()) for line in f if line.strip()]

    if len(nums1) != len(nums2):
        print(f"FAIL: files have different lengths: {len(nums1)} vs {len(nums2)}")
        return False

    max_diff = 0.0
    for i, (a, b) in enumerate(zip(nums1, nums2)):
        diff = abs(a - b)
        max_diff = max(max_diff, diff)

    if max_diff > tol:
        print(f"FAIL: Maximum difference {max_diff:.2e} exceeds tolerance {tol:.2e}")
        return False
    else:
        print(f"PASS: Maximum difference {max_diff:.2e} within tolerance {tol:.2e}")
        return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_files.py file1.txt file2.txt [tolerance]")
        sys.exit(1)

    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    passed = compare_files(sys.argv[1], sys.argv[2], tol)
    sys.exit(0 if passed else 1)
