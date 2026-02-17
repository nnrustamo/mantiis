#!/bin/bash

# Test runner for mantiis LBM code
# Compiles tests, runs them, and compares outputs to gold files

set -e  # Exit on first error during compilation

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Set library path for OpenCV
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH

# Configuration
NP=2                                    # Number of MPI processes
THREADS=1                               # Number of OpenMP threads
MULTIBLOCK="False"                      # Multiblock setting
GRID_SIZE=64                            # Grid size
ITERS=100                               # Number of iterations
TEST_DIR="test_input_output"            # Test input/output directory
GOLD_DIR="$TEST_DIR/gold"               # Gold files directory
COMPARE_SCRIPT="../scripts/compare_files.py"
TOLERANCE=1e-10

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track overall test status
TESTS_PASSED=0
TESTS_FAILED=0

print_header() {
    echo ""
    echo "=============================================="
    echo "$1"
    echo "=============================================="
}

print_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
    TESTS_PASSED=$((TESTS_PASSED + 1))
}

print_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    TESTS_FAILED=$((TESTS_FAILED + 1))
}

print_info() {
    echo -e "${YELLOW}[INFO]${NC} $1"
}

# Compare a single file pair
compare_file() {
    local test_file="$1"
    local gold_file="$2"
    local description="$3"
    
    if [[ ! -f "$test_file" ]]; then
        print_fail "$description: Test output file not found ($test_file)"
        return 1
    fi
    
    if [[ ! -f "$gold_file" ]]; then
        print_fail "$description: Gold file not found ($gold_file)"
        return 1
    fi
    
    if python3 "$COMPARE_SCRIPT" "$test_file" "$gold_file" "$TOLERANCE"; then
        print_pass "$description"
        return 0
    else
        print_fail "$description"
        return 1
    fi
}

# Run a test suite
run_test_suite() {
    local test_name="$1"
    local executable="$2"
    local gold_subdir="$3"
    
    print_header "Running $test_name test"
    
    # Run the test
    print_info "Executing: mpirun -np $NP ./$executable $THREADS $MULTIBLOCK $GRID_SIZE $TEST_DIR/ $ITERS"
    
    set +e  # Don't exit on test failure
    mpirun -np $NP ./$executable $THREADS $MULTIBLOCK $GRID_SIZE "$TEST_DIR/" $ITERS > /dev/null 2>&1
    local run_status=$?
    set -e
    
    if [[ $run_status -ne 0 ]]; then
        print_fail "$test_name: Execution failed with exit code $run_status"
        return 0  # Continue to next test
    fi
    
    print_info "Comparing results to gold files..."
    
    # Compare ux.txt (|| true prevents exit on failure)
    compare_file "$TEST_DIR/ux.txt" "$GOLD_DIR/$gold_subdir/ux.txt" "$test_name: ux.txt" || true
    
    # Compare uy.txt
    compare_file "$TEST_DIR/uy.txt" "$GOLD_DIR/$gold_subdir/uy.txt" "$test_name: uy.txt" || true
    
    # Compare rho.txt
    compare_file "$TEST_DIR/rho.txt" "$GOLD_DIR/$gold_subdir/rho.txt" "$test_name: rho.txt" || true
}

# Main execution
print_header "Building tests"
# print_info "Running make clean && make..."
# make clean
make

print_header "Test Configuration"
echo "  MPI processes:  $NP"
echo "  OpenMP threads: $THREADS"
echo "  Grid size:      $GRID_SIZE x $GRID_SIZE"
echo "  Iterations:     $ITERS"
echo "  Tolerance:      $TOLERANCE"
echo "  Test directory: $TEST_DIR"

# Run transport test
run_test_suite "Transport" "transport_test" "transport"

# Run adsorption test
run_test_suite "Adsorption" "adsorption_test" "adsorption"

# Summary
print_header "Test Summary"
echo -e "Passed: ${GREEN}$TESTS_PASSED${NC}"
echo -e "Failed: ${RED}$TESTS_FAILED${NC}"

if [[ $TESTS_FAILED -eq 0 ]]; then
    echo ""
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo ""
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
