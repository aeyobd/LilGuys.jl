#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile
import subprocess

import re



def main():
    args = parse_args()
    params = retrieve_parameters(args.param)
    if not args.resubmit:
        clean_dir(params['OutputDir'])

    set_default_mem_time(args, params)
    set_default_name(args)

    exe = os.path.join(os.getenv("GADGET_PATH"), args.executable)
    if args.resubmit == -1:
        script = create_sbatch_script(args, exe, 1)
    elif args.resubmit is not None:
        script = create_sbatch_script(args, exe, f"2 {args.resubmit}")
    else:
        script = create_sbatch_script(args, exe, "")

    print(script)
    tmpfile = _create_temp_script(script)

    # Submit the job
    submit_job(tmpfile)
    os.remove(tmpfile)



def clean_dir(directory):
    """
    check if a directory exists and if it does, ask the user if they want to overwrite it
    """
    if os.path.exists(directory):
        ans = input(f"Overwrite {directory}? (y/n)")
        if ans.lower() == 'y':
            os.removedirs(directory)
        else:
            sys.exit(1)
    os.makedirs(directory)


def _create_temp_script(script):
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.sh', prefix='slurm_', dir='.') as tmpfile:
        tmpfile.write(script)
        print(f"Created temporary script at {tmpfile.name}")

    return tmpfile.name



def parse_args():
    parser = argparse.ArgumentParser(description='Submit a job to SLURM with default or specified parameters.')
    parser.add_argument('executable', type=str, help='the executable to run (see $GADGET_PATH)')
    parser.add_argument('--ntasks', type=int, default=2, help='Number of CPUs per task (default: 4)')
    parser.add_argument('--mem', type=str, default=None, help='Memory per job (e.g., 4G, 16G) (default: 4G)')
    parser.add_argument('--time', type=str, default=None, help='Wall-clock time limit (HH:MM:SS) (default: 01:00:00)')
    parser.add_argument('--name', type=str, default=None, help='SLURM job name (default: pwd)')
    parser.add_argument('--partition', type=str, default='cosma', help='SLURM partition (default: short)')
    parser.add_argument('--account', type=str, default=os.getenv("SLURM_ACCOUNT"), help='slurm account')
    parser.add_argument('-p', '--param', type=str, default="param.txt", help='parameter file')
    parser.add_argument('--resubmit', default=None, type=int, help='resubmit the job?')
    args = parser.parse_args()

    return args



def submit_job(script_path):
    print("submitting")
    result = subprocess.run(['sbatch', script_path], capture_output=True, text=True)
    print(result.stderr)
    print(result.stdout)
    print("submitted")



def create_sbatch_script(args, executable, resubmit):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    return f"""#!/bin/bash
#SBATCH --job-name  {args.name}
#SBATCH --output    %J.out
#SBATCH --time      {args.time}
#SBATCH --ntasks    {args.ntasks}
#SBATCH --mem       {args.mem}
#SBATCH --partition {args.partition}
#SBATCH --account   {args.account}

source {script_dir}/slurm_header.sh
mpirun -np $SLURM_NTASKS {executable} param.txt {resubmit}

scontrol show job $SLURM_JOB_ID
"""




def set_default_name(args):
    if args.name is None:
        pwd = os.getcwd()
        args.name = extract_path_after_substring(pwd, 'models/')

    return args


def extract_path_after_substring(full_path, substring):
    # Find the index of the substring in the path
    index = full_path.find(substring)
    if index != -1:
        # Add the length of substring to index to start after the substring
        return full_path[index + len(substring):]
    else:
        # Return None or an appropriate message if the substring is not found
        return None


def set_default_mem_time(args, params):
    if args.mem is None:
        args.mem = params['MaxMemSize'] + 'MB'
    if args.time is None:
        time = int(params['TimeLimitCPU'])
        hours = time // 3600
        minutes = (time % 3600) // 60
        seconds = time % 60
        args.time = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

    return args


def retrieve_parameters(filename):
    """
    Retrieve parameters from a file.

    Parameters
    ----------
    filename : str

    Returns
    -------
    parameters : dict
        A dictionary of parameters
    """

    parameters = {}
    with open(filename, 'r') as file:
        for line in file:
            # Strip comments
            key, value = _read_param_line(line)
            if key:
                parameters[key] = value

    return parameters



def _read_param_line(line):
    """
    Read a line from a parameter file.

    Parameters
    ----------
    line : str
        The line to read

    Returns
    -------
    key : str
        The key of the parameter

    value : str
        The value of the parameter
    """

    key = None
    value = None

    # Strip comments
    line = line.split('#')[0].split('%')[0].strip()
    if line:  # If there's anything left after stripping comments and whitespace
        # Split on whitespace and take the first two elements as key and value
        parts = line.split()
        if len(parts) >= 2:
            key = parts[0].strip()
            value = parts[1].strip()

    return key, value



if __name__ == '__main__':
    main()
