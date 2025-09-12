#!/bin/bash

# Define base directory
BASE_DIR="/home/nrustamo/projects/HNO/image_generation/images"

# Define the list of directories to loop through
directories=(
    "nx1024_ny1024_scaling0.4_sigma20_octave20"
    "nx1024_ny1024_scaling0.4_sigma20_octave40"
    "nx1024_ny1024_scaling0.4_sigma50_octave20"
    "nx1024_ny1024_scaling0.4_sigma50_octave40"
    "nx1024_ny1024_scaling0.8_sigma20_octave20"
    "nx1024_ny1024_scaling0.8_sigma20_octave40"
    "nx1024_ny1024_scaling0.8_sigma50_octave20"
    "nx1024_ny1024_scaling0.8_sigma50_octave40"
)

# Loop through each directory and run the command
for dir in "${directories[@]}"; do
    # Construct the full path to the directory and ensure it ends with a slash
    full_path="$BASE_DIR/$dir/"
    
    # Run the command with the current directory
    echo "Running command for directory: $full_path"
    mpirun -np 12 ./mantiis 1 false 1024 "$full_path" 200000
    
    # Optional: Add a small sleep if you need to throttle execution
    # sleep 1
done
