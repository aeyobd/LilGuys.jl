#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile
import subprocess
import shutil

import re



def main():
    args = parse_args()
    params = retrieve_parameters(args.param)
    if not args.resubmit:
        clean_dir(params['OutputDir'])

    set_defaults(args, params)

    exe = os.path.join(os.getenv("GADGET_PATH"), args.executable)
    if args.resubmit == -1:
        script = create_sbatch_script(args, exe, 1)
    elif args.resubmit is not None:
        script = create_sbatch_script(args, exe, f"2 {args.resubmit}")
    else:
        script = create_sbatch_script(args, exe, "")

    submit_job(script)


def parse_args():
    parser = argparse.ArgumentParser(description='Submit a job to SLURM with default or specified parameters.')
    parser.add_argument('executable', type=str, 
                        help='the executable to run  in $GADGET_PATH')
    parser.add_argument('--ntasks', type=int, default=4, 
                        help='Number of CPUs per task (default: 4)')
    parser.add_argument('--mem', type=str, default=None, 
                        help='Memory per job (default: from param.txt)')
    parser.add_argument('--time', type=str, default=None, 
                        help='(wall) time limit (HH:MM:SS) (default: from param.txt)')
    parser.add_argument('--name', type=str, default=None, 
                        help='SLURM job name (default: based on directory name)')
    parser.add_argument('--partition', type=str, default='cosma', 
                        help='SLURM partition (default: cosma)')
    parser.add_argument('--account', type=str, default=os.getenv("SLURM_ACCOUNT"), 
                        help='SLURM account (default: $SLURM_ACCOUNT)')
    parser.add_argument('-p', '--param', type=str, default="param.txt", 
                        help='Gadget parameter file (default: param.txt)')
    parser.add_argument('--resubmit', default=None, type=int, 
                        help='if specified, then we are resubmitting a job. Options are -1 (for last restartfile) or a nonnegative integer for the snapshot to restart from.')
    args = parser.parse_args()

    return args


def clean_dir(directory):
    """
    check if a directory exists and if it does, ask the user if they want to overwrite it
    """
    if os.path.exists(directory):
        ans = input(f"Overwrite {directory}? (y/n) ")
        if ans.lower() == 'y':
            shutil.rmtree(directory)
        else:
            sys.exit(1)
    os.makedirs(directory)


def create_sbatch_script(args, executable, resubmit):
    """
    Create a SLURM script for submitting a job.
    """

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



def submit_job(script):
    print("submitting")

    result = subprocess.run(
        ['sbatch'], input=script,
        capture_output=True, text=True
    )

    print(result.stderr)
    print(result.stdout)



def set_defaults(args, params):
    _set_default_mem(args, params)
    _set_default_time(args, params)
    _set_default_name(args)


def _set_default_name(args):
    if args.name is None:
        pwd = os.getcwd()
        args.name = _extract_path_after_substring(pwd, 'models/')

    return args


def _extract_path_after_substring(full_path, substring):
    # Find the index of the substring in the path
    index = full_path.find(substring)
    if index != -1:
        # Add the length of substring to index to start after the substring
        return full_path[index + len(substring):]
    else:
        # Return None or an appropriate message if the substring is not found
        return full_path


def _set_default_mem(args, params):
    """
    If the args do not specify `mem`, set it to the value in the parameter file.
    """
    if args.mem is None:
        args.mem = params['MaxMemSize'] + 'MB'

    return args


def _set_default_time(args, params):
    """
    If the args do not specify `time`, set it to the value in the parameter file.
    """
    if args.time is None:
        time = int(params['TimeLimitCPU'])
        hours = time // 3600
        minutes = (time % 3600) // 60
        seconds = time % 60
        args.time = f"{hours:02d}:{minutes:02d}:{seconds:02d}"



def retrieve_parameters(filename):
    """
    Retrieve parameters from a Gadget parameter file.
    The first word is the key, and the value can be seperated by arbitrary whitespace.
    Comments are denoted by '#' or '%' anywhere in the line.

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

            key, value = _read_param_line(line)

            if key:
                parameters[key] = value

    return parameters



def _read_param_line(line):
    """
    Read a line from a Gadget parameter file.

    Parameters
    ----------
    line : str
        The line to read

    Returns
    -------
    key, value: str, str
        The name and value of the parameter. Returns None, None if the line is empty.
    """

    key = None
    value = None

    line_no_comments = line.split('#')[0].split('%')[0].strip()

    if line_no_comments:
        parts = line_no_comments.split()
        if len(parts) == 2:
            key = parts[0].strip()
            value = parts[1].strip()
        else:
            print(f"Warning: Ignoring line '{line_no_comments}'")

    return key, value



def _create_temp_script(script):
    """writes the text to a temporary file and returns the filename. Not used, can remove..."""
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.sh', prefix='slurm_', dir='.') as tmpfile:
        tmpfile.write(script)
        print(f"Created temporary script at {tmpfile.name}")

    return tmpfile.name




if __name__ == '__main__':
    main()
