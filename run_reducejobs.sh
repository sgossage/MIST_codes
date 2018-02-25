#!/bin/bash

# ${1} is e.g. p0.30 (desired [Fe/H]) & ${2} is v/vcrit e.g. 0.1. ${3} would be convective core ovsh, e.g. 0.016
#
# Usage:
#
#    >> ./run_reducejobs.sh p0.30 0.1
#
# to process ~ feh_p0.30_afe_p0.0_vvcrit0.1.
#

if [[ -z "${3}" ]]; then
    srund reduce_jobs.py MIST_v1.0/feh_${1}_afe_p0.0_vvcrit${2}_TP
else
    srund reduce_jobs.py MIST_v1.0/feh_${1}_afe_p0.0_vvcrit${2}_cov${3}_TP
fi

