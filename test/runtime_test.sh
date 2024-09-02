#!/bin/bash

# Parameters
omp_threads=(1 5 10)
grid_options=("true" "false")
nx_sizes=(1024 512 256 128 64)

# Path to the executable
executable="../main"

# Log file
log_file="run_log.txt"

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

#
> "$log_file"
original_dir=$(pwd)
cd "$(dirname "$executable")" || exit
if [ ! -x "$(basename "$executable")" ]; then
    echo -e "${RED}Error: Executable $(basename "$executable") not found or not executable${NC}" | tee -a "$original_dir/$log_file"
    exit 1
fi

# Loop over each combination of parameters
for threads in "${omp_threads[@]}"; do
    for grid in "${grid_options[@]}"; do
        for nx in "${nx_sizes[@]}"; do
            echo "Running with OpenMP threads=$threads, MultiGrid=$grid, Nx=$nx" | tee -a "$original_dir/$log_file"
            
            echo "Current directory: $(pwd)" | tee -a "$original_dir/$log_file"
            
            # Run the executable with the parameters and log output
            {
                echo "Running command: $(basename "$executable") $threads $grid $nx"
                ./$(basename "$executable") $threads $grid $nx
                exit_status=$?
                echo "Exit status: $exit_status"
                
                if [ $exit_status -eq 0 ]; then
                    echo -e "${GREEN}[SUCCESS]${NC}" | tee -a "$original_dir/$log_file"
                else
                    echo -e "${RED}[ERROR]${NC}" | tee -a "$original_dir/$log_file"
                fi
            } >> "$original_dir/$log_file" 2>&1
            
            if [ $? -ne 0 ]; then
                echo -e "${RED}Error: The executable did not run successfully with OpenMP threads=$threads, Grid Initialization=$grid, Nx=$nx${NC}" | tee -a "$original_dir/$log_file"
                exit 1
            fi

            echo -e "${GREEN}Completed run with OpenMP threads=$threads, Grid Initialization=$grid, Nx=$nx [SUCCESS]${NC}" | tee -a "$original_dir/$log_file"
        done
    done
done

cd "$original_dir" || exit
