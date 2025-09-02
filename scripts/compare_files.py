import sys

def compare_files(file1, file2, tol=1e-20):
    with open(file1, 'r') as f:
        nums1 = [float(line.strip()) for line in f if line.strip()]

    with open(file2, 'r') as f:
        nums2 = [float(line.strip()) for line in f if line.strip()]

    if len(nums1) != len(nums2):
        print(f"Error: files have different lengths: {len(nums1)} vs {len(nums2)}")
        return

    max_diff = 0.0
    for i, (a, b) in enumerate(zip(nums1, nums2)):
        diff = abs(a - b)
        if diff > tol:
            print(f"Line {i+1}: {a} vs {b}, |diff| = {diff}")
        max_diff = max(max_diff, diff)

    print(f"Maximum absolute difference: {max_diff}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_files.py file1.txt file2.txt")
        sys.exit(1)

    compare_files(sys.argv[1], sys.argv[2])
