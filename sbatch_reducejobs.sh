#!/bin/bash

#SBATCH -J reducejobs
#SBATCH -n 2                    # Number of cores
#SBATCH -N 1                    # Core on machine
#SBATCH -t 12:00:00             # H:M:S
#SBATCH --mem 16000              # 8Gb Memory request
#SBATCH -p itc_cluster,conroy
#SBATCH -o reducejobs_o/reducejobs.o
#SBATCH -e reducejobs_e/reducejobs.e

# ${1} is e.g. p0.30 (desired [Fe/H]) & ${2} is v/vcrit e.g. 0.1. ${3} would be convective core ovsh, e.g. 0.016
#
# Usage:
#
#    >> ./run_reducejobs.sh p0.30 0.1
#
# to process ~ feh_p0.30_afe_p0.0_vvcrit0.1.
#

if [[ -z "${3}" ]]; then
    ./reduce_jobs.py MIST_v1.0/feh_${1}_afe_p0.0_vvcrit${2}_TP
else
    ./reduce_jobs.py MIST_v1.0/feh_${1}_afe_p0.0_vvcrit${2}_cov${3}_TP
fi

