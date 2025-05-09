import re

def process_run_log(file_path):
    # Regular expressions to match relevant lines
    re_threads = re.compile(r'Running with OpenMP threads=(\d+)')
    re_multigrid = re.compile(r'MultiGrid=(true|false)')
    re_domain_size = re.compile(r'Nx=(\d+)')
    re_active_cells = re.compile(r'The number of active cells: (\d+)')
    re_simulation_time = re.compile(r'Simulation run time \(milliseconds\) (\d+)')
    re_total_run_time = re.compile(r'Total run time \(milliseconds\) (\d+)')
    
    results = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i]
        print(line)
        if 'Running with OpenMP threads=' in line:
            # Extract thread count, multigrid setting, and domain size
            threads = re_threads.search(line).group(1)
            multigrid = re_multigrid.search(line).group(1)
            domain_size = re_domain_size.search(line).group(1)
            
            active_cells = None
            simulation_time = None
            total_run_time = None

        i += 1
        line = lines[i]
            
        # Extract relevant details
        while i < len(lines) and 'Running with OpenMP threads=' not in line:
            if 'The number of active cells:' in line:
                active_cells = re_active_cells.search(line).group(1)
                i += 1
                line = lines[i]
            if 'Simulation run time (milliseconds)' in line:
                simulation_time = re_simulation_time.search(line).group(1)
                i += 1
                line = lines[i]
            if 'Total run time (milliseconds)' in line:
                total_run_time = re_total_run_time.search(line).group(1)
                i += 1
                line = lines[i]
            i += 1
            if i >= len(lines):
                break
            else:
                line = lines[i]
            
        if active_cells and simulation_time and total_run_time:
            time_per_iteration = float(simulation_time) / 100
            
            results.append({
                'threads': int(threads),
                'MultiGrid': multigrid,
                'Nx': int(domain_size),
                'Active Cells': int(active_cells),
                'Simulation Time (s)': float(simulation_time) / 1000,
                'Total Run Time (s)': float(total_run_time)/ 1000,
                'Time Per Iteration (s)': float(time_per_iteration) / 1000
            })
        else:
            i += 1
    
    return results

file_path = 'run_log.txt'
data = process_run_log(file_path)

import pandas as pd
excel_file_path = 'run_log_results.xlsx'
df = pd.DataFrame(data)
df.to_excel(excel_file_path, index=False, engine='openpyxl')
print(f"Data successfully saved to {excel_file_path}")

